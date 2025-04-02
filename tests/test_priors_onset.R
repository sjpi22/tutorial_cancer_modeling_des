###########################  Unit test: Estimating simple prior distributions   ################
#
#  Objective: Estimate plausible parameter priors for simple time-to-event  
#  distributions between cancer states given data for each state
########################### <<<<<>>>>> ##############################################

rm(list = ls()) # Clean environment
options(scipen = 999) # View data without scientific notation

#### 1.Libraries and functions  ==================================================

###### 1.1 Load packages
library(readxl)
library(tidyverse)
library(lhs)
library(assertthat)
library(data.table)
library(cobs) # For fitting constrained B-splines
library(sandwich) # For heteroskedasticity-robust linear model standard errors
library(lhs) # For Latin hypercube sampling
library(tools) # For title case
library(survival)

###### 1.2 Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)

# Define nonnegative smooth function for clinical cancer incidence (note: fit_spline must be a global variable)
ir_clinical <- function(x) {
  return(pmax(0, predict(fit_spline, x)[, "fit"]))
}


#### 2. General parameters ========================================================

###### 2.1 Configurations
# Run file to process configurations
source("configs/process_configs.R")
list2env(configs, envir = .GlobalEnv)

# Get list of relevant output file paths and load to global environment
l_filepaths <- update_config_paths("files_tests", configs$paths)
list2env(l_filepaths, envir = .GlobalEnv)

###### 2.2 File paths
path_truth <- "_ground_truth"
file_distr <- file.path(path_truth, "true_params.xlsx")
file.surv <- file.path(paths$data, "relative_survival_cancer.csv")

###### 2.3 Other parameters
conf_level <- 0.95 # For generating bounds
multiplier_bounds <- 0.2
age_interval <- 0.25
v_cols <- c("targets", "ci_lb", "ci_ub")
v_colors <- c("green", "red", "orange")
var_index <- "target_index"
n_samp <- 50 # Number of cohorts to simulate for checking coverage
seed <- 123
eps <- 1e-8


#### 3. Load data  ===========================================

# Set seed
set.seed(seed) 

# Load ground truth model parameters
l_params_model <- do.call(load_model_params, c(
  modifyList(params_model,
             list(file.surv = NULL),
             keep.null = T),
  list(seed = NULL,
       file.distr = file_distr)
))

# Load survival data
df_surv <- load_surv_data(params_model$file.surv)

# Load targets
l_targets <- load_calibration_targets(params_calib$l_params_outcome)

# Process targets
for (target in names(params_calib$l_outcome_params)) {
  # Process incidence data
  if (params_calib$l_outcome_params[[target]][["outcome_type"]] == "incidence") {
    # Calculate approximate confidence intervals for incidence
    l_targets[[target]] <- l_targets[[target]] %>%    
      mutate(ci_lb = pmax(0, targets - qnorm((1 + conf_level)/2)*se),
             ci_ub = targets + qnorm((1 + conf_level)/2)*se)
    
    # Rescale incidence values by unit
    for (val in c(v_cols, "se")) {
      l_targets[[target]][[val]] <- l_targets[[target]][[val]] / l_targets[[target]]$unit
    }
  }
}

# Set variables dependent on parameters
max_age <- l_params_model$max_age
v_ages <- seq(0, max_age, 0.25)

# Ordered list of prevalence targets
if (params_model$lesion_state == T) { 
  l_prev_targets <- "prevalence_lesion"
  l_prev_targets <- c(l_prev_targets, "prevalence")
  idx_preclinical <- -1
} else {
  l_prev_targets <- "prevalence"
  idx_preclinical <- 1
}

# Get variable for disease onset
var_onset <- paste0("time_H_", l_params_model$v_states[2])

# Get variable for censoring incidence
var_censor <- params_calib$l_params_outcome$incidence$lit_params$censor_var

# Simulate ground truth data for comparisons
m_cohort <- run_base_model(l_params_model)

# Separate patient and lesion data as necessary
if (is.data.table(m_cohort)) {
  m_patients <- m_cohort
} else {
  m_patients <- m_cohort$patient_level
  m_lesions <- m_cohort$lesion_level
}

# Save true survival distribution
d_time_C_Dc_true <- l_params_model[paste0("d_time_C", l_params_model$v_cancer, "_Dc")]

# Swap empirical survival distribution
for (i in unique(df_surv$stage)) {
  # Filter survival data to stage at diagnosis
  temp_df_surv <- df_surv[df_surv$stage == i, ]
  
  # Calculate probability mass function from CDF
  probs <- diff(temp_df_surv$pct_died)
  probs <- c(probs, 1 - sum(probs))
  
  # Create distribution data
  l_params_model[[paste0("d_time_C", i, "_Dc")]] <- list(
    distr = "empirical", 
    params = list(xs = temp_df_surv$years_from_dx, 
                  probs = probs, 
                  max_x = max_age), 
    src = "known")
}


#### 4. Analysis  ===========================================

##### 4.0 Plot calibration targets
df_prevalence <- l_targets[l_prev_targets]
df_incidence <- l_targets$incidence %>%
  mutate(targets = targets / unit,
         se = se / unit) %>%
  mutate(ci_lb = targets - qnorm((1 + conf_level)/2)*se,
         ci_ub = targets + qnorm((1 + conf_level)/2)*se)

# Plot prevalence
plot(df_prevalence[[1]][[var_index]], df_prevalence[[1]][[v_cols[1]]], 
     ylim = c(0, max(df_prevalence[[1]][[v_cols[3]]])),
     xlab = "Age", ylab = "Prevalence")
arrows(df_prevalence[[1]][[var_index]], df_prevalence[[1]][[v_cols[2]]], 
       df_prevalence[[1]][[var_index]], df_prevalence[[1]][[v_cols[3]]], 
       length = 0.05, angle = 90, code = 3)

if (length(l_prev_targets) > 1) {
  # Plot prevalence of 2nd state
  points(df_prevalence[[2]][[var_index]], df_prevalence[[2]][[v_cols[1]]], col = "red")
  arrows(df_prevalence[[2]][[var_index]], df_prevalence[[2]][[v_cols[2]]], 
         df_prevalence[[2]][[var_index]], df_prevalence[[2]][[v_cols[3]]], 
         col = "red", length = 0.05, angle = 90, code = 3)
}

# Plot incidence
plot(df_incidence[[var_index]], df_incidence[[v_cols[1]]], 
     ylim = c(0, max(df_incidence[[v_cols[3]]])),
     xlab = "Age", ylab = "Incidence")
arrows(df_incidence[[var_index]], df_incidence[[v_cols[2]]], 
       df_incidence[[var_index]], df_incidence[[v_cols[3]]], length = 0.05, angle = 90, code = 3)

##### 4.1 Derive average time from diagnosis to death by fitting splines
# Initialize simulation population 
m_death <- data.table()

# Sample every three values of ages for spline knots
v_knots_Dc <- l_params_model[["d_time_C1_Dc"]]$params$xs
v_knots_Dc <- c(v_knots_Dc[seq(1, length(v_knots_Dc), 3)], max_age)

# Fit splines to survival
par(mfrow = c(length(l_params_model$v_cancer)/2, 2))
d_time_C_Dc_spline <- list()
for (i in l_params_model$v_cancer) {
  # Calculate cumulative percentage dead from survival distribution
  pct_Dc <- cumsum(l_params_model[[paste0("d_time_C", i, "_Dc")]]$params$probs)
  
  # Convert to cumulative hazard
  chaz_Dc <- -log(1 - pct_Dc)
  
  # Fit spline to survival data
  fit_spline_Dc <- cobs(
    x = l_params_model[[paste0("d_time_C", i, "_Dc")]]$params$xs[-1],
    y = head(chaz_Dc, -1), 
    constraint = c("increase"),
    knots = v_knots_Dc,
    pointwise = matrix(c(0, 0, 0), ncol = 3))
  
  # Plot for verification
  plot(fit_spline_Dc, xlim = c(0, 20),
       xlab = "Proportion dead", ylab = "Time to death",
       main = paste0("Stage ", i))
  
  # Calculate CDF from 0 to max_age bounding at 1
  cdf_Dc <- pmax(0, pmin(1 - exp(-predict(fit_spline_Dc, v_ages)[, "fit"]), 1))
  
  # Calculate probability mass function from CDF
  probs <- diff(cdf_Dc)
  probs <- c(probs, 1 - sum(probs))
  
  # Create spline distribution data
  d_time_C_Dc_spline[[i]] <- list(distr = "empirical", 
                                  params = list(xs = v_ages, 
                                                probs = probs, 
                                                max_x = max_age))
  
  # Simulate time to death with spline
  m_death[, paste0("time_C", i, "_Dc") := query_distr(
    "r", l_params_model$n_cohort, 
    d_time_C_Dc_spline[[i]]$distr, 
    d_time_C_Dc_spline[[i]]$params
  )]
}

# Get expectation of time to death
mean_time_C_Dc <- colMeans(m_death)
mean_Dc <- sum(l_targets$stage_distr$targets * mean_time_C_Dc)

# Create rough estimate exponential distribution for time from diagnosis to death from cancer
d_time_C_Dc_est <- list(distr = "exp", params = list(rate = 1/mean_Dc))

# Compare to true expectation
mean_time_C_Dc_true <- sapply(d_time_C_Dc_true, function(distr) integrate(function(x) x*query_distr("d", x, distr$distr, distr$params), lower = 0, upper = Inf)$value)
mean_Dc_true <- sum(l_targets$stage_distr$targets * mean_time_C_Dc_true)
print("Estimated vs. true mean time to death from cancer:")
print(paste(round(mean_Dc, 3), "vs", round(mean_Dc_true, 3)))

##### 4.2 Derive time-to-event distribution from targets using splines
# Sample every three values of ages for spline knots
v_knots <- df_incidence$age_start[seq(2, length(df_incidence$age_start), 3)]

# Calculate clinical cancer probability corresponding to confidence intervals of cancer incidence
l_p_clinical <- list()
par(mfrow = c(1, 3))
for (val in v_cols) {
  # Fit spline for clinical cancer incidence
  fit_spline <- cobs(
    x = df_incidence[[var_index]],
    y = df_incidence[[val]], 
    constraint = "increase",
    w = 1/(df_incidence$se^2),
    knots = c(0, v_knots, max(df_incidence$age_end)))
  
  # Labels for plotting spline fit
  if (val == v_cols[1]) {
    label <- "mean"
  } else if (val == v_cols[2]) {
    label <- "CI LB"
  } else {
    label <- "CI UB"
  }
  
  # If cancer patients were considered not at risk after diagnosis in incidence calculate,
  # convert clinical cancer cumulative hazard to cumulative probability,
  # but if they were considered still at risk, the cumulative hazard is actually the cumulative probability
  if (var_censor == "time_H_C") {
    # Plot spline fit
    plot(fit_spline, xlab = "Age", ylab = "Incidence", main = paste("Spline fit:", label),
         xlim = c(0, max(df_incidence$age_end) + 1))
    arrows(df_incidence[[var_index]], df_incidence$ci_lb, 
           df_incidence[[var_index]], df_incidence$ci_ub, length = 0.05, angle = 90, code = 3)
    
    # Integrate spline to estimate cumulative hazard of clinical cancer at ages in preclinical cancer data
    chaz_clinical <- sapply(df_prevalence[[-1]][[var_index]],
                            function(t) integrate(ir_clinical, lower = 0, upper = t)[["value"]])
    
    # Calculate cumulative probability of clinical cancer
    l_p_clinical[[val]] <- 1 - exp(-chaz_clinical)
  } else if (var_censor == "time_H_D") {
    # Estimate percent dead from cancer
    pct_dead <- sapply(df_incidence[[var_index]],
                       function(t) integrate(function(u) ir_clinical(u)*query_distr("p", t - u, d_time_C_Dc_est$distr, d_time_C_Dc_est$params), lower = 0, upper = t)[["value"]])
    
    # Recalculate probability density of clinical cancer, scaling by percent dead from cancer
    pdf_clinical <- df_incidence[[val]] * (1 - pct_dead)
    
    # Refit spline to clinical cancer probability density
    fit_spline <- cobs(
      x = df_incidence[[var_index]],
      y = pdf_clinical, 
      constraint = "increase",
      w = 1/(df_incidence$se^2),
      knots = c(0, v_knots, max(df_incidence$age_end) + 1))
    
    # Plot spline fit
    plot(fit_spline, xlab = "Age", ylab = "Probability", main = paste("Spline fit:", label),
         xlim = c(0, max(df_incidence$age_end) + 1))
    arrows(df_incidence[[var_index]], df_incidence[[v_cols[2]]], 
           df_incidence[[var_index]], df_incidence[[v_cols[3]]], length = 0.05, angle = 90, code = 3)
    
    # Calculate cumulative probability of clinical cancer by integrating over density
    l_p_clinical[[val]] <- sapply(df_prevalence[[-1]][[var_index]],
                                  function(t) integrate(ir_clinical, lower = 0, upper = t)[["value"]])
  }
  
  # Calculate probabilities for time to preclinical cancer by adding and scaling preclinical prevalence and clinical probability
  df_prevalence[[idx_preclinical]][[paste("p", val, sep = "_")]] <- df_prevalence[[idx_preclinical]][[val]] * (1 - l_p_clinical[[val]]) + l_p_clinical[[val]]
  
  # Calculate probabilities for time to precancerous lesion by adding and scaling lesion prevalence and preclinical probability
  # Note - multiplied by p_clinical rather than p_preclinical because screening study includes preclinical cases in the denominator
  if (params_model$lesion_state == T) {
    df_prevalence[[1]][[paste("p", val, sep = "_")]] <- df_prevalence[[1]][[val]] * (1 - l_p_clinical[[val]]) + df_prevalence[[-1]][[paste("p", val, sep = "_")]]
  }
}

# Plot estimates against truth
pdf(file_fig_prior_cdf, width = 6, height = 6)
par(mfrow = c(1, 1))
for (val in v_cols) {
  # Plot targets
  if (val == v_cols[1]) {
    # Plot estimated disease onset CDF
    plot(df_prevalence[[1]]$target_index, df_prevalence[[1]]$p_targets,
         col = v_colors[1], type = "p", ylim = c(0, max(df_prevalence[[1]]$p_ci_ub)),
         xlab = "Age", ylab = "Probability")
    
    # Plot prevalence for comparison
    arrows(df_prevalence[[1]]$target_index, df_prevalence[[1]]$ci_lb, 
           df_prevalence[[1]]$target_index, df_prevalence[[1]]$ci_ub, length = 0.05, angle = 90, code = 3,
           col = v_colors[1], lty = 3)
    
    # Plot preclinical cancer CDF if model includes lesion state
    if (params_model$lesion_state == T) {
      points(df_prevalence[[-1]]$target_index, df_prevalence[[-1]]$p_targets, col = v_colors[3])
      
      # Plot prevalence for comparison
      arrows(df_prevalence[[-1]]$target_index, df_prevalence[[-1]]$ci_lb, 
             df_prevalence[[-1]]$target_index, df_prevalence[[-1]]$ci_ub, length = 0.05, angle = 90, code = 3,
             col = v_colors[3], lty = 3)
    }
    
    # Plot estimated clinical cancer CDF
    points(df_prevalence[[1]]$target_index, l_p_clinical[[val]], col = v_colors[2])
    
  } else {
    # Plot error bounds
    lines(df_prevalence[[1]]$target_index, df_prevalence[[1]][[paste("p", val, sep = "_")]], col = v_colors[1], type = "l", lty = 2)
    if (params_model$lesion_state == T) {
      lines(df_prevalence[[-1]]$target_index, df_prevalence[[-1]][[paste("p", val, sep = "_")]], col = v_colors[3], type = "l", lty = 2)
    }
    lines(df_prevalence[[1]]$target_index, l_p_clinical[[val]], col = v_colors[2], type = "l", lty = 2)
  }
}

# Calculate true disease distributions for plotting
p_onset_true <- query_distr("p", v_ages, l_params_model[[(paste0("d_", var_onset))]]$distr, l_params_model[[(paste0("d_", var_onset))]]$params) 
lines(v_ages, p_onset_true, col = v_colors[1])

# Plot true state 2 (if applicable) and clinical cancer CDF
for (i in 1:length(l_prev_targets)) {
  # Get variable name
  var_dx <- paste0("time_H_", l_params_model$v_states[i + 2])
  
  # Fit KM curve to disease state data, censoring at death time
  fit_surv <- survfit(Surv(m_patients[, pmin(get(var_dx), time_H_D, na.rm = T)], 
                           m_patients[, pmin(get(var_dx), Inf, na.rm = T) < time_H_D]) ~ 1)
  
  # Calculate CDF at each time point
  cdf_surv <- summary(fit_surv)
  df_cdf_surv = data.frame(
    time = cdf_surv$time,
    pct = 1 - cdf_surv$surv
  )
  
  # Plot
  lines(df_cdf_surv$time, df_cdf_surv$pct, col = v_colors[4 - i])
}

# Add legend
v_labs <- c("Preclinical", "Clinical")
if (params_model$lesion_state == T) {
  v_labs <- c("Lesion", v_labs)
  v_colors_plot <- v_colors
} else {
  v_colors_plot <- head(v_colors, -1)
}

legend("topleft",
       legend = c(v_labs, "Estimated CDF", "CI of estimated CDF", "True CDF", "Prevalence"), 
       col = c(v_colors_plot, rep("black", 4)), pch = c(rep(15, length(v_labs)), 1, rep(NA, 3)),
       lty = c(rep(NA, length(v_labs)), NA, 2, 1, 3),
       bty = "n", border = F)
dev.off()

##### 4.3 Fit Weibull distribution to time to disease onset
# Calculate Weibull x transformation
df_prevalence[[1]] <- df_prevalence[[1]] %>%
  mutate(x_transformed = log(get(var_index)))

# Calculate Weibull y transformation
for (val in v_cols) {
  df_prevalence[[1]][[paste("y", val, sep = "_")]] <- log(-log(1 - df_prevalence[[1]][[paste("p", val, sep = "_")]]))
}

# Get Weibull estimates using weighted linear regression
fit_lm_mean <- lm(y_targets ~ x_transformed, data = df_prevalence[[1]], 
                  weights = 1/(y_ci_ub - y_ci_lb)^2) # Weight scales CI to log level

# Get shape and scale estimates for time to disease onset distribution
coefs_onset <- fit_lm_mean$coefficients
shape_onset <- coefs_onset[2]
scale_onset <- exp(-coefs_onset[1]/shape_onset)

# Set estimated distribution for time to disease onset
d_time_onset_est <- list()
d_time_onset_est[[v_cols[1]]] <- list(distr = "weibull", params = list(shape = shape_onset, scale = scale_onset))

# Calculate heteroskedasticity-robust standard errors
vcov <- vcovHC(fit_lm_mean)
stderrorHC <- sqrt(diag(vcov))

# Calculate confidence interval of estimates
coefs_onset_lb <- coefs_onset - qnorm((1 + conf_level)/2) * stderrorHC
coefs_onset_ub <- coefs_onset + qnorm((1 + conf_level)/2) * stderrorHC
coefs_onset_ci <- cbind(coefs_onset_lb, coefs_onset_ub)

# Plot transformations to validate Weibull fit - should be linear
for (val in v_cols) {
  if (val == v_cols[1]) {
    plot(df_prevalence[[1]]$x_transformed, df_prevalence[[1]][[paste("y", val, sep = "_")]], 
         xlab = "log(age)", ylab = "log(-log(1 - CDF))",
         main = "Transformations for Weibull regression")
  } else {
    points(df_prevalence[[1]]$x_transformed, df_prevalence[[1]][[paste("y", val, sep = "_")]])
  }
}
abline(coefs_onset[1], coefs_onset[2], col = "red")
abline(coefs_onset_lb[1], coefs_onset_ub[2], col = "red", lty = 2)
abline(coefs_onset_ub[1], coefs_onset_lb[2], col = "red", lty = 2)

# Convert CIs to shape and scale - note, dividing min intercept by max slope 
# and vice versa for more accurate coverage of line
shape_onset_ci <- coefs_onset_ci[2, ]
scale_onset_ci <- rev(exp(-rev(coefs_onset_ci[1,])/shape_onset_ci))
d_time_onset_est[[v_cols[3]]] <- list(distr = "weibull", params = list(shape = shape_onset_ci[1], scale = scale_onset_ci[1])) # Reversed order since higher Weibull scale --> lower CDF
d_time_onset_est[[v_cols[2]]] <- list(distr = "weibull", params = list(shape = shape_onset_ci[2], scale = scale_onset_ci[2]))

# Plot true vs. fitted Weibull parameters
pdf(file_fig_prior_onset, width = 6, height = 6)
curve(pweibull(x, l_params_model[[paste0("d_", var_onset)]]$params$shape, l_params_model[[paste0("d_", var_onset)]]$params$scale), from = 0, to = 10,
      xlab = "Time to disease onset", ylab = "Density", main = "Weibull fit to disease onset data")
lines(0:10, pweibull(0:10, shape_onset, scale_onset), col = "red")
lines(0:10, pweibull(0:10, shape_onset_ci[1], scale_onset_ci[1]), col = "red", lty = 2)
lines(0:10, pweibull(0:10, shape_onset_ci[2], scale_onset_ci[2]), col = "red", lty = 2)
legend("topleft",
       legend = c("True CDF", "Initial estimate of CDF", "Prior bound CDFs"), 
       col = c("black", "red", "red"),
       lty = c(1, 1, 2),
       bty = "n", border = F)
dev.off()

##### 4.4 Calculate time from lesion to cancer onset
if (params_model$lesion_state == T) {
  # Calculate SE of preclinical cancer CDF
  se_preclin <- (df_prevalence[[2]]$ci_ub - df_prevalence[[2]]$ci_lb)/2/qnorm((1 + conf_level)/2)
  
  # Set initial distribution for time from lesion to cancer onset
  d_time_onset_cancer <- list(distr = "weibull", params = list(shape = 1, scale = 1))
  
  # Initialize matrix to hold parameters
  df_params_L_P <- matrix(nrow = 0, ncol = 2)
  
  for (val in v_cols) {
    # Find parameters that optimize fit
    fit_cdf <- optim(par = c(3, 10), fn = function(x) {
      x <- as.list(x)
      names(x) <- names(d_time_onset_cancer$params)
      obj_fn_cdf(x, d_time_onset_est[[val]], d_time_onset_cancer, 
                 df_prevalence[[2]][[var_index]], 
                 df_prevalence[[2]][[paste0("p_", val)]],
                 se_preclin,
                 method = "ll")
    }, control = list(fnscale = -1))
    
    # Output warning if optimization did not converge
    if (fit_cdf$convergence != 0) {
      warning("Optimization did not converge")
    }
    
    # Append to matrix
    df_params_L_P <- rbind(df_params_L_P, fit_cdf$par)
  }
  
  # Convert matrix to data frame
  df_params_L_P <- data.frame(v_cols, df_params_L_P)
  names(df_params_L_P) <- c("value", names(d_time_onset_cancer$params))
  
  #### Check recovery of distribution ####
  # Get variable name
  var_dx <- "time_H_P"
  
  # Fit KM curve to disease state data, censoring at death time
  fit_surv <- survfit(Surv(m_patients[, pmin(get(var_dx), time_H_D, na.rm = T)], 
                           m_patients[, pmin(get(var_dx), Inf, na.rm = T) < time_H_D]) ~ 1)
  
  # Calculate CDF at each time point
  cdf_surv <- summary(fit_surv)
  df_cdf_surv = data.frame(
    time = cdf_surv$time,
    pct = 1 - cdf_surv$surv
  )
  
  # Plot true CDF
  plot(df_cdf_surv$time, df_cdf_surv$pct, col = v_colors[4 - i], type = "l",
       xlab = "Age", ylab = "Probability", main = "Time to preclinical cancer")
  
  # Plot estimated CDFs
  for (val in v_cols) {
    # Update distribution for time from lesion to cancer onset
    params_L_P <- df_params_L_P[df_params_L_P$value == val, -1]
    names(params_L_P) <- names(d_time_onset_cancer$params)
    d_time_onset_cancer$params <- modifyList(d_time_onset_cancer$params, params_L_P)
    
    # Calculate CDF of sum at indices
    p_pred_P <- calc_cdf_sum(df_prevalence[[2]][[var_index]],
                             d_time_onset_est[[val]],
                             d_time_onset_cancer)
    if (val == "targets") {
      points(df_prevalence[[2]][[var_index]], p_pred_P)
    } else {
      lines(df_prevalence[[2]][[var_index]], p_pred_P, lty = 2)
    }
  }
}
  
# The following sections have only been tested on the no-lesion model
if (params_model$lesion_state == F) {
  ##### 4.4 Fit exponential distribution to time from preclinical to clinical cancer
  # Create initial exponential distribution estimate for time from preclinical to clinical cancer
  # with rate = 1 / age difference between the clinical and preclinical distribution quantiles
  # at the maximum probability observed in the clinical distribution
  age_P_max_C <- query_distr("q", max(l_p_clinical[[v_cols[1]]]), d_time_onset_est$distr, d_time_onset_est$params)
  rate_P_C_init <- 1/(max(df_prevalence[[idx_preclinical]][[var_index]]) - age_P_max_C)
  d_time_P_C_est <- list(distr = "exp", params = list(rate = rate_P_C_init))
  
  # Create objective function for difference between predicted and observed clinical cancer CDF
  obj_fn_clinical <- function(rate_P_C, val) {
    # Update distribution for time from preclinical to clinical
    d_time_P_C_est$params$rate <- rate_P_C
    
    # Calculate difference between observed and estimated clinical cancer CDF
    diff_cdf <- l_p_clinical[[val]] - sapply(df_prevalence[[idx_preclinical]][[var_index]], function(t) {
      # Calculates sum of the two distributions
      integrate(function(u)
        query_distr("p", t - u, d_time_onset_est$distr, d_time_onset_est$params) * 
          query_distr("d", u, d_time_P_C_est$distr, d_time_P_C_est$params),
        lower = 0, upper = t)[["value"]]
    })
    
    # Return sum of squared differences
    return(sum(diff_cdf^2))
  }
  
  # Perform quick optimization to refine parameter estimate for time from preclinical to clinical cancer
  l_rate_P_C_est <- list()
  for (val in v_cols) {
    l_rate_P_C_est[[val]] <- optimize(function(x) obj_fn_clinical(x, val = val), 
                                      interval = c(d_time_P_C_est$params$rate / 4, d_time_P_C_est$params$rate * 4))$minimum
  }
  
  ##### 4.5 Check coverage of targets
  # Create dataframe of parameter priors
  df_prior_vals <- data.frame(rbind(
    unname(c(shape_onset, shape_onset_ci, d_time_onset$params$shape)),
    unname(c(scale_onset, scale_onset_ci, d_time_onset$params$scale)),
    unname(c(unlist(l_rate_P_C_est), d_time_P_C$params$rate))
  ))
  colnames(df_prior_vals) <- c("init", "min", "max", "truth")
  
  # Create map of prior distributions
  df_priors <- data.frame(
    var_name = c("d_time_H_P", "d_time_H_P", "d_time_P_C"),
    param_name = c("shape", "scale", "rate"),
    df_prior_vals
  )
  df_priors$var_id <- paste(df_priors$var_name, df_priors$param_name, sep = ".")
  
  # Get number of params to calibrate and number of samples
  n_param <- nrow(df_priors)
  
  # Sample unit Latin Hypercube
  m_lhs_unit <- randomLHS(n_samp, n_param)
  
  # Rescale to min/max of each parameter
  m_param_samp <- matrix(nrow = n_samp, ncol = n_param)
  for (i in 1:n_param) {
    m_param_samp[, i] <- qunif(
      m_lhs_unit[, i],
      min = df_priors$min[i],
      max = df_priors$max[i])
  }
  colnames(m_param_samp) <- df_priors$var_id
  
  # Set estimated model parameters
  l_params_est <- update_param_from_map(l_params_model, df_prior_vals$init, df_priors)
  
  # Simulate model with updated parameters and check coverage of targets
  l_outputs_est <- list(params_to_targets(l_params_est, NULL, df_priors, var_censor, v_ages, alpha, ir_unit = 1))
  v_target_max <- c(prevalence = 0, incidence = 0)
  for (i in 1:n_samp) {
    # Calculate targets for parameter set
    l_targets_temp <- params_to_targets(l_params_est, m_param_samp[i, ], df_priors, var_censor, v_ages, alpha, ir_unit = 1)
    
    # Add targets to list
    l_outputs_est <- c(l_outputs_est,
                       list(l_targets_temp))
    
    # Keep track of maximum value for plotting
    v_target_max <- pmax(v_target_max, c(
      max(l_targets_temp$prevalence[[v_cols[1]]]),
      max(l_targets_temp$incidence[[v_cols[1]]])))
  }
  
  # Plot coverage of targets and visually verify
  par(mfrow = c(1, 2))
  for (target in c("prevalence", "incidence")) {
    # Plot main estimate
    plot(l_outputs_est[[1]][[target]][[var_index]], l_outputs_est[[1]][[target]][[v_cols[1]]], type = "l", col = "blue", lwd = 2,
         ylim = c(0, v_target_max[target]),
         xlab = "Age", ylab = toTitleCase(target))
    
    # Plot LHS estimates
    for (i in 1:n_samp) {
      lines(l_outputs_est[[i + 1]][[target]][[var_index]], l_outputs_est[[i + 1]][[target]][[v_cols[1]]], col = "gray", lty = 2)
    }
    
    # Plot truth on top
    points(l_targets[[target]][[var_index]], l_targets[[target]][[v_cols[1]]])
    arrows(l_targets[[target]][[var_index]], l_targets[[target]]$ci_lb, 
           l_targets[[target]][[var_index]], l_targets[[target]]$ci_ub, length = 0.05, angle = 90, code = 3)
  }
  
  # Verify that prior bounds cover the truth
  df_priors <- df_priors %>%
    mutate(coverage = (truth >= min & truth <= max))
  if (sum(df_priors$coverage == FALSE) > 0) {
    stop("Prior bounds do not cover the truth")
  } else {
    print("Prior bounds cover the truth")
  }
}