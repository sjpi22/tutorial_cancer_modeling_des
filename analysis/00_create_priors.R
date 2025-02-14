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

###### 2.2 Other parameters
alpha <- 0.01
multiplier_bounds <- 0.2
age_interval <- 0.25
v_cols <- c("targets", "ci_lb", "ci_ub")


#### 3. Load data  ===========================================

# Load targets
l_true_targets <- load_calibration_targets(params_calib$l_outcome_params)

# Process incidence data
for (target in names(params_calib$l_outcome_params)) {
  if (params_calib$l_outcome_params[[target]][["outcome_type"]] == "incidence") {
    # Calculate approximate confidence intervals for incidence
    l_true_targets[[target]] <- l_true_targets[[target]] %>%    
      mutate(ci_lb = pmax(0, targets - qnorm(1 - alpha/2)*se),
             ci_ub = targets + qnorm(1 - alpha/2)*se)
    
    # Rescale incidence values by unit
    for (val in c(v_cols, "se")) {
      l_true_targets[[target]][[val]] <- l_true_targets[[target]][[val]] / l_true_targets[[target]]$unit
    }
  }
}

# Load ground truth model parameters
l_params_model <- do.call(load_model_params, c(
  params_model,
  list(seed = params_calib$seed_calib)
))

# Set variables dependent on parameters
max_age <- l_params_model$max_age
v_ages <- seq(0, max_age, 0.25)
if (params_model$lesion_state == T) {
  target1 <- "prevalence_lesion"
  target2 <- "prevalence"
} else {
  target1 <- "prevalence"
}

# Set seed
set.seed(params_calib$seed_calib)


#### 4. Derive prior distribution for time to disease onset  ===========================================

##### 4.1 Derive average time from diagnosis to death by fitting splines

# Initialize simulation population 
m_death <- data.table()

# Sample every three values of ages for spline knots
v_knots_Dc <- l_params_model$d_time_C_Dc[[1]]$params$xs
v_knots_Dc <- c(v_knots_Dc[seq(1, length(v_knots_Dc), 4)], max_age)

# Fit splines to survival
par(mfrow = c(length(l_params_model$v_cancer)/2, 2))
d_time_C_Dc <- list()
for (i in l_params_model$v_cancer) {
  # Calculate cumulative percentage dead from survival distribution
  pct_died <- cumsum(l_params_model$d_time_C_Dc[[i]]$params$probs)
  
  # Fit spline to survival data
  fit_spline_Dc <- cobs(
    x = l_params_model$d_time_C_Dc[[i]]$params$xs[-1],
    y = head(pct_died, -1), 
    constraint = c("increase", "concave"),
    knots = v_knots_Dc,
    pointwise = matrix(c(0, 0, 0), ncol = 3))
  
  # Plot for verification
  plot(fit_spline_Dc, xlim = c(0, 20), ylim = c(0, 1),
       xlab = "Proportion dead", ylab = "Time to death",
       main = paste0("Stage ", i))
  
  # Calculate CDF from 0 to max_age bounding at 1
  cdf_Dc <- pmax(0, pmin(predict(fit_spline_Dc, v_ages)[, "fit"], 1))
  
  # Calculate probability mass function from CDF
  probs <- diff(cdf_Dc)
  probs <- c(probs, 1 - sum(probs))
  
  # Create spline distribution data
  d_time_C_Dc[[i]] <- list(distr = "empirical", 
                           params = list(xs = v_ages, 
                                         probs = probs, 
                                         max_x = max_age))
  
  # Simulate time to death with spline
  m_death[, paste0("time_C", i, "_Dc") := query_distr(
    "r", l_params_model$n_cohort, 
    d_time_C_Dc[[i]]$distr, 
    d_time_C_Dc[[i]]$params
  )]
}

# Get expectation of time to death
mean_Dc <- sum(l_true_targets$stage_distr$targets * colMeans(m_death))

##### 4.2 Derive time-to-event distribution from targets using splines

# Sample every three values of ages for spline knots
v_knots <- l_true_targets$incidence$age_start[seq(2, length(l_true_targets$incidence$age_start), 3)]

# Calculate clinical cancer probability corresponding to confidence intervals of cancer incidence
l_p_clinical <- list()
par(mfrow = c(1, 3))
for (val in v_cols) {
  # Fit spline for clinical cancer incidence
  fit_spline <- cobs(
    x = l_true_targets$incidence$target_index,
    y = l_true_targets$incidence[[val]], 
    constraint = "increase",
    w = 1/(l_true_targets$incidence$se^2),
    knots = c(0, v_knots, max(l_true_targets$incidence$age_end) + 1))
  
  # Plot spline fit
  if (val == "targets") {
    label <- "mean"
  } else if (val == "ci_lb") {
    label <- "CI LB"
  } else {
    label <- "CI UB"
  }
  
  # If cancer patients were considered not at risk after diagnosis in incidence calculate,
  # convert clinical cancer cumulative hazard to cumulative probability,
  # but if they were considered still at risk, the cumulative hazard is actually the cumulative probability
  if (params_calib$l_outcome_params$incidence$censor_var == "time_H_C") {
    # Plot spline fit
    plot(fit_spline, xlab = "Age", ylab = "Incidence", main = paste("Spline fit:", label),
         xlim = c(0, max(l_true_targets$incidence$age_end) + 1))
    
    # Integrate spline to estimate cumulative hazard of clinical cancer at ages in preclinical cancer data
    chaz_clinical <- sapply(l_true_targets$prevalence$target_index,
                            function(t) integrate(ir_clinical, lower = 0, upper = t)[["value"]])
    
    # Calculate cumulative probability of clinical cancer
    l_p_clinical[[val]] <- 1 - exp(-chaz_clinical)
  } else if (params_calib$l_outcome_params$incidence$censor_var == "time_H_D") {
    # Estimate percent dead from cancer
    pct_dead <- sapply(l_true_targets$incidence$target_index,
                       function(t) integrate(function(u) ir_clinical(u)*(t - pmax(u, t - mean_Dc))/mean_Dc, lower = 0, upper = t)[["value"]])
    
    # Recalculate probability density of clinical cancer, scaling by percent dead from cancer
    pdf_clinical <- l_true_targets$incidence[[val]] * (1 - pct_dead)
    
    # Refit spline to clinical cancer probability density
    fit_spline <- cobs(
      x = l_true_targets$incidence$target_index,
      y = pdf_clinical, 
      constraint = "increase",
      w = 1/(l_true_targets$incidence$se^2),
      knots = c(0, v_knots, max(l_true_targets$incidence$age_end) + 1))
    
    # Plot spline fit
    plot(fit_spline, xlab = "Age", ylab = "Probability", main = paste("Spline fit:", label),
         xlim = c(0, max(l_true_targets$incidence$age_end) + 1))
    arrows(l_true_targets$incidence$target_index, l_true_targets$incidence$ci_lb, 
           l_true_targets$incidence$target_index, l_true_targets$incidence$ci_ub, length = 0.05, angle = 90, code = 3)
    
    # Calculate cumulative probability of clinical cancer by integrating over density
    l_p_clinical[[val]] <- sapply(l_true_targets$prevalence$target_index,
                                  function(t) integrate(ir_clinical, lower = 0, upper = t)[["value"]])
  }
  
  # Calculate probabilities for time to preclinical cancer by adding and scaling preclinical prevalence and clinical probability
  l_true_targets$prevalence[[paste("p", val, sep = "_")]] <- l_true_targets$prevalence[[val]] * (1 - l_p_clinical[[val]]) + l_p_clinical[[val]]
  
  # Calculate probabilities for time to precancerous lesion by adding and scaling lesion prevalence and preclinical probability
  if (params_model$lesion_state == T) {
    l_true_targets[[target1]][[paste("p", val, sep = "_")]] <- l_true_targets[[target1]][[val]] * (1 - l_p_clinical[[val]]) + l_p_clinical[[val]]
  }
}

# Plots
par(mfrow = c(1, 1))
for (val in v_cols) {
  # Plot targets
  if (val == "targets") {
    # Plot estimated disease onset CDF
    plot(l_true_targets[[target1]]$target_index, l_true_targets[[target1]]$p_targets,
         type = "l", ylim = c(0, max(l_true_targets[[target1]]$p_ci_ub)),
         xlab = "Age", ylab = "Probability")
    
    # Plot preclinical cancer CDF if model includes lesion state
    if (params_model$lesion_state == T) {
      lines(l_true_targets[[target2]]$target_index, l_true_targets[[target2]]$p_targets, col = "purple", type = "l")
    }
    
    # Plot estimated clinical cancer CDF
    lines(l_true_targets$prevalence$target_index, l_p_clinical[[val]], col = "red", type = "l")
  } else {
    # Plot error bounds
    lines(l_true_targets[[target1]]$target_index, l_true_targets[[target1]][[paste("p", val, sep = "_")]], col = "black", type = "l", lty = 2)
    if (params_model$lesion_state == T) {
      lines(l_true_targets[[target2]]$target_index, l_true_targets[[target2]][[paste("p", val, sep = "_")]], col = "purple", type = "l", lty = 2)
    }
    lines(l_true_targets$prevalence$target_index, l_p_clinical[[val]], col = "red", type = "l", lty = 2)
  }
}

# Plot preclinical cancer prevalence for comparison
arrows(l_true_targets$prevalence$target_index, l_true_targets$prevalence$ci_lb, 
       l_true_targets$prevalence$target_index, l_true_targets$prevalence$ci_ub, length = 0.05, angle = 90, code = 3,
       col = "blue")

# Add legend
v_labs <- c("Preclinical CDF", "Clinical CDF", "Preclinical Prevalence")
if (params_model$lesion_state == T) {
  v_labs <- c("Lesion CDF", v_labs)
  v_colors <- c("black", "purple", "red", "blue")
} else {
  v_colors <- c("black", "red", "blue")
}
legend("topleft", legend = v_labs, col = v_colors, lty = rep(1, length(v_labs)))

##### 4.5 Fit Weibull distribution to time to disease onset

# Calculate Weibull x transformation
l_true_targets[[target1]] <- l_true_targets[[target1]] %>%
  mutate(x_transformed = log(target_index))

# Calculate Weibull y transformation
for (val in v_cols) {
  l_true_targets[[target1]][[paste("y", val, sep = "_")]] <- log(-log(1 - l_true_targets[[target1]][[paste("p", val, sep = "_")]]))
}

# Plot transformations to validate Weibull fit - should be linear
for (val in v_cols) {
  if (val == v_cols[1]) {
    plot(l_true_targets[[target1]]$x_transformed, l_true_targets[[target1]][[paste("y", val, sep = "_")]], 
         xlab = "log(age)", ylab = "log(-log(1 - CDF))",
         main = "Transformations for Weibull regression",
         ylim = range(l_true_targets[[target1]][paste("y", v_cols, sep = "_")]))
  } else {
    points(l_true_targets[[target1]]$x_transformed, l_true_targets[[target1]][[paste("y", val, sep = "_")]])
  }
}

# Get Weibull estimates using weighted linear regression
fit_lm_mean <- lm(y_targets ~ x_transformed, data = l_true_targets[[target1]], 
                  weights = 1/(y_ci_ub - y_ci_lb)^2) # Weight scales CI to log level

# Get shape and scale estimates for time to preclinical cancer distribution
coefs_H_P <- fit_lm_mean$coefficients
shape_H_P <- coefs_H_P[2]
scale_H_P <- exp(-coefs_H_P[1]/shape_H_P)

# Set estimated distribution for time to preclinical cancer
d_time_H_P_est <- list(distr = "weibull", params = list(shape = shape_H_P, scale = scale_H_P))

# Calculate heteroskedasticity-robust standard errors
vcov <- vcovHC(fit_lm_mean)
stderrorHC <- sqrt(diag(vcov))

# Calculate confidence interval of estimates
coefs_H_P_lb <- coefs_H_P - qnorm(1 - alpha/2) * stderrorHC
coefs_H_P_ub <- coefs_H_P + qnorm(1 - alpha/2) * stderrorHC
coefs_H_P_ci <- cbind(coefs_H_P_lb, coefs_H_P_ub)

# Convert CIs to shape and scale
shape_H_P_ci <- coefs_H_P_ci[2, ]
scale_H_P_ci <- rev(exp(-coefs_H_P_ci[1,]/shape_H_P_ci))

# Increase CIs by multiplier bounds
v_multiplier <- c(1 - multiplier_bounds, 1 + multiplier_bounds)
shape_H_P_bounds <- shape_H_P_ci * v_multiplier
scale_H_P_bounds <- scale_H_P_ci * v_multiplier

##### 4.4 Update prior distribution

# Update prior dataframe
df_priors <- read.csv(params_calib$file_prior)
df_priors[df_priors$var_id == paste0("d_time_H_", l_params_model$v_states[2], ".shape"), c("min", "max")] <- as.list(shape_H_P_bounds)
df_priors[df_priors$var_id == paste0("d_time_H_", l_params_model$v_states[2], ".scale"), c("min", "max")] <- as.list(scale_H_P_bounds)

# Overwrite prior file
write.csv(df_priors, file = params_calib$file_priors, row.names = FALSE)
