###########################  Create prior distributions  ##########################
#
#  Objective: Refine data-driven prior distributions for model parameters
########################### <<<<<>>>>> #########################################

rm(list = ls()) # Clean environment
options(scipen = 999) # View data without scientific notation

#### 1.Libraries and functions  ==================================================

###### 1.1 Load packages
library(readxl)
library(tidyverse)
library(assertthat)
library(data.table)
library(cobs) # For fitting constrained B-splines
library(sandwich) # For heteroskedasticity-robust linear model standard errors

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

# Extract relevant parameters from configs
params_model <- configs$params_model
params_calib <- configs$params_calib

###### 2.2 Other parameters
conf_level <- 0.95 # For generating bounds
multiplier_bounds <- 0.2
age_interval <- 0.25
v_cols <- c("targets", "ci_lb", "ci_ub")
v_colors <- c("green", "red", "orange")
var_index <- "target_index"


#### 3. Load data  ===========================================

# Set seed
set.seed(params_calib$seed_calib)

# Load model parameters
l_params_model <- do.call(load_model_params, c(
  params_model,
  list(seed = params_calib$seed_calib)
))

# Load targets
l_targets <- load_calibration_targets(params_calib$l_outcome_params)

# Process incidence data
for (target in names(params_calib$l_outcome_params)) {
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
var_censor <- params_calib$l_outcome_params$incidence$lit_params$censor_var


#### 4. Derive prior distribution for time to disease onset  ===========================================

##### 4.0 Plot calibration targets
df_prevalence <- l_targets[l_prev_targets]
df_incidence <- l_targets$incidence

# Plot prevalence
par(mfrow = c(1, 2))
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
v_knots_Dc <- l_params_model$d_time_C1_Dc$params$xs
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

# Plots
par(mfrow = c(1, 1))
for (val in v_cols) {
  # Plot targets
  if (val == v_cols[1]) {
    # Plot estimated disease onset CDF
    plot(df_prevalence[[1]][[var_index]], df_prevalence[[1]]$p_targets,
         col = v_colors[1], type = "p", ylim = c(0, max(df_prevalence[[1]]$p_ci_ub)),
         xlab = "Age", ylab = "Probability")
    
    # Plot prevalence for comparison
    arrows(df_prevalence[[1]][[var_index]], df_prevalence[[1]]$ci_lb, 
           df_prevalence[[1]][[var_index]], df_prevalence[[1]]$ci_ub, length = 0.05, angle = 90, code = 3,
           col = v_colors[1], lty = 3)
    
    # Plot preclinical cancer CDF if model includes lesion state
    if (params_model$lesion_state == T) {
      points(df_prevalence[[-1]][[var_index]], df_prevalence[[-1]]$p_targets, col = v_colors[3])
      
      # Plot prevalence for comparison
      arrows(df_prevalence[[-1]][[var_index]], df_prevalence[[-1]]$ci_lb, 
             df_prevalence[[-1]][[var_index]], df_prevalence[[-1]]$ci_ub, length = 0.05, angle = 90, code = 3,
             col = v_colors[3], lty = 3)
    }
    
    # Plot estimated clinical cancer CDF
    points(df_prevalence[[1]][[var_index]], l_p_clinical[[val]], col = v_colors[2])
    
  } else {
    # Plot error bounds
    lines(df_prevalence[[1]][[var_index]], df_prevalence[[1]][[paste("p", val, sep = "_")]], col = v_colors[1], type = "l", lty = 2)
    if (params_model$lesion_state == T) {
      lines(df_prevalence[[-1]][[var_index]], df_prevalence[[-1]][[paste("p", val, sep = "_")]], col = v_colors[3], type = "l", lty = 2)
    }
    lines(df_prevalence[[1]][[var_index]], l_p_clinical[[val]], col = v_colors[2], type = "l", lty = 2)
  }
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
       legend = c(v_labs, "Estimated CDF", "Estimate CI", "Prevalence"), 
       col = c(v_colors_plot, rep("black", 3)), pch = c(rep(0, length(v_labs)), 1, rep(NA, 2)),
       lty = c(rep(NA, length(v_labs)), NA, 2, 3),
       bty = "n", border = F, ncol = 2)

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
d_time_onset_est <- list(distr = "weibull", params = list(shape = shape_onset, scale = scale_onset))

# Calculate heteroskedasticity-robust standard errors
vcov <- vcovHC(fit_lm_mean)
stderrorHC <- sqrt(diag(vcov))

# Calculate confidence interval of estimates
coefs_onset_lb <- coefs_onset - qnorm((1 + conf_level)/2) * stderrorHC
coefs_onset_ub <- coefs_onset + qnorm((1 + conf_level)/2) * stderrorHC
coefs_onset_ci <- cbind(coefs_onset_lb, coefs_onset_ub)

# Convert CIs to shape and scale - note, dividing min intercept by max slope 
# and vice versa for more accurate coverage of line
shape_onset_ci <- coefs_onset_ci[2, ]
scale_onset_ci <- rev(exp(-rev(coefs_onset_ci[1,])/shape_onset_ci))

# Expand bounds by multiplier
v_multipliers <- c(1 - multiplier_bounds, 1 + multiplier_bounds)
shape_onset_bounds <- shape_onset_ci * v_multipliers
scale_onset_bounds <- scale_onset_ci * v_multipliers

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

# Plot fitted Weibull distribution of disease onset
plot(v_ages, pweibull(v_ages, shape_onset, scale_onset), type = "l",
     xlab = "Age", ylab = "Probability", main = "Fitted Weibull distribution")
lines(v_ages, pweibull(v_ages, shape_onset_ci[1], scale_onset_ci[1]), lty = 2)
lines(v_ages, pweibull(v_ages, shape_onset_ci[2], scale_onset_ci[2]), lty = 2)

# Plot against estimated distribution of disease onset
points(df_prevalence[[1]][[var_index]], df_prevalence[[1]][[paste0("p_", v_cols[1])]])
arrows(df_prevalence[[1]][[var_index]], df_prevalence[[1]][[paste0("p_", v_cols[2])]], 
       df_prevalence[[1]][[var_index]], df_prevalence[[1]][[paste0("p_", v_cols[3])]], 
       length = 0.05, angle = 90, code = 3)

# Plot against prevalence targets
points(df_prevalence[[1]][[var_index]], df_prevalence[[1]][[v_cols[1]]], col = "blue")
arrows(df_prevalence[[1]][[var_index]], df_prevalence[[1]][[v_cols[2]]], 
       df_prevalence[[1]][[var_index]], df_prevalence[[1]][[v_cols[3]]], 
       length = 0.05, angle = 90, code = 3, col = "blue")

##### 4.4 Update prior distribution
# Update prior dataframe
df_priors <- read.csv(params_calib$file_prior)
df_priors[df_priors$var_id == paste0("d_time_H_", l_params_model$v_states[2], ".shape"), c("min", "max")] <- as.list(shape_onset_bounds)
df_priors[df_priors$var_id == paste0("d_time_H_", l_params_model$v_states[2], ".scale"), c("min", "max")] <- as.list(scale_onset_bounds)

# Overwrite prior file
write.csv(df_priors, file = params_calib$file_priors, row.names = FALSE)
