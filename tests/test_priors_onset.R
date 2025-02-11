###########################  Unit test: Estimating simple prior distributions   ################
#
#  Objective: Estimate plausible parameter priors for simple time-to-event  
#  distributions between cancer states given data for each state
########################### <<<<<>>>>> ##############################################


#### 1.Libraries and functions  ==================================================
# Clean environment
rm(list = ls())
options(scipen = 99)

library(readxl)
library(tidyverse)
library(lhs)
library(assertthat)
library(data.table)
library(cobs) # For fitting constrained B-splines
library(sandwich) # For heteroskedasticity-robust linear model standard errors
library(lhs) # For Latin hypercube sampling
library(tools) # For title case

###### 1.1 Load functions =================================================
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)

# Function for running simplified cancer model to sample patients
run_cancer_model <- function(l_params_model) {
  m_times <- with(l_params_model, {
    # Sample patients
    m_times <- data.table(
      id = 1:n_cohort
    )
    
    # Generate time to event data
    m_times[, `:=` (time_H_P = query_distr("r", .N, d_time_H_P$distr, d_time_H_P$params),
                    time_P_C = query_distr("r", .N, d_time_P_C$distr, d_time_P_C$params),
                    time_C_D = query_distr("r", .N, d_time_C_D$distr, d_time_C_D$params),
                    time_H_Do = query_distr("r", .N, d_time_H_Do$distr, d_time_H_Do$params))]
    m_times[, time_H_C := time_H_P + time_P_C]
    m_times[, time_H_Dc := time_H_C + time_C_D]
    
    # Summarize time to death from all causes
    if (include_Do) {
      m_times[, time_H_D := pmin(time_H_Do, time_H_Dc)]
    } else {
      m_times[, time_H_D := time_H_Dc]
    }
    
    # Calculate censor time for screening studies
    m_times[, time_screen_censor := pmin(time_H_D, time_H_C)]
  })
  
  return(m_times)
}

# Function for calculating calibration targets
calc_targets <- function(m_times, var_censor, v_ages, 
                         alpha = 0.05, ir_unit = 100000,
                         n_prevalence_sample = NULL) {
  # Calculate prevalence of preclinical cancer
  if (is.null(n_prevalence_sample)) n_prevalence_sample <- nrow(m_times)
  df_prevalence <- calc_prevalence(m_times, "time_H_P", "time_H_C", "time_screen_censor", v_ages = v_ages,
                                   n_screen_sample = n_prevalence_sample) %>%
    mutate(age_midpoint = (age_start + age_end)/2,
           ci_lb = pmax(0, value - qnorm(1 - alpha/2)*se),
           ci_ub = value + qnorm(1 - alpha/2)*se)
  
  # Calculate incidence of clinical cancer
  df_incidence <- calc_incidence(m_times, "time_H_C", var_censor, v_ages = v_ages,
                                 rate_unit = ir_unit) %>%
    mutate(age_midpoint = (age_start + age_end)/2,
           ci_lb = pmax(0, value - qnorm(1 - alpha/2)*se),
           ci_ub = value + qnorm(1 - alpha/2)*se)
  
  return(list(prevalence = df_prevalence, incidence = df_incidence))
}

# Function for converting vector of parameters to calibration outputs
params_to_targets <- function(
    l_params_model, 
    v_params, 
    df_priors, 
    var_censor, 
    v_ages, 
    alpha = 0.05, 
    ir_unit = 100000, 
    n_prevalence_sample = NULL
) {
  # Update parameters
  if (!is.null(v_params)) {
    l_params_updated <- update_param_from_map(l_params_model, v_params, df_priors)
  } else {
    l_params_updated <- l_params_model
  }
  
  # Run cancer model
  m_times <- run_cancer_model(l_params_updated)
  
  # Calculate calibration targets
  l_targets <- calc_targets(m_times, var_censor, v_ages, alpha, ir_unit, n_prevalence_sample)
  
  return(l_targets)
}

# Define nonnegative smooth function for clinical cancer incidence (note: fit_spline must be a global variable)
ir_clinical <- function(x) {
  return(pmax(0, predict(fit_spline, x)[, "fit"]))
}


#### 2. General parameters ========================================================

###### 2.1 model parameters
n_cohort <- 100000
include_Do <- TRUE # Flag to include death from other causes
v_ages <- 0:10
var_censor <- "time_H_D"
alpha <- 0.01
v_cols <- c("value", "ci_lb", "ci_ub")
n_samp <- 50 # Number of cohorts to simulate for checking coverage
set.seed(123) # Set seed for reproducibility

# Set distribution variables
d_time_H_P <- list(distr = "weibull", params = list(shape = 2, scale = 10))
d_time_P_C <- list(distr = "exp", params = list(rate = 0.16))
d_time_C_D <- list(distr = "exp", params = list(rate = 0.01))
d_time_H_Do <- list(distr = "weibull", params = list(shape = 5, scale = 8))


#### 3. Analysis ========================================================

##### 3.1 Generate and plot calibration targets

# Create model parameters list
l_params_model <- list(
  n_cohort = n_cohort, 
  d_time_H_P = d_time_H_P, 
  d_time_P_C = d_time_P_C, 
  d_time_C_D = d_time_C_D,
  d_time_H_Do = d_time_H_Do,
  include_Do = include_Do)

# Sample patients
m_times <- run_cancer_model(l_params_model)

# Calculate prevalence and incidence targets
l_targets <- calc_targets(m_times, var_censor, v_ages, alpha, ir_unit = 1)
df_prevalence <- l_targets$prevalence
df_incidence <- l_targets$incidence

# Plot prevalence
plot(df_prevalence$age_midpoint, df_prevalence$value, 
     ylim = range(df_prevalence$value) + c(-diff(range(df_prevalence$value)), diff(range(df_prevalence$value))))
arrows(df_prevalence$age_midpoint, df_prevalence$ci_lb, 
       df_prevalence$age_midpoint, df_prevalence$ci_ub, length = 0.05, angle = 90, code = 3)

# Plot incidence
plot(df_incidence$age_midpoint, df_incidence$value, 
     ylim = c(0, max(df_incidence$value)))
arrows(df_incidence$age_midpoint, df_incidence$ci_lb, 
       df_incidence$age_midpoint, df_incidence$ci_ub, length = 0.05, angle = 90, code = 3)


##### 3.2 Derive time-to-event distribution from targets using splines

# Sample every three values of ages for spline knots
v_knots <- df_incidence$age_start[seq(2, length(df_incidence$age_start), 3)]

# Calculate clinical cancer probability corresponding to confidence intervals of cancer incidence
l_p_clinical <- list()
par(mfrow = c(1, 3))
for (val in v_cols) {
  # Fit spline for clinical cancer incidence
  fit_spline <- cobs(
    x = df_incidence$age_midpoint,
    y = df_incidence[[val]], 
    constraint = "increase",
    w = 1/(df_incidence$se^2),
    knots = c(0, v_knots, max(df_incidence$age_end) + 1))
  
  # Plot spline fit
  if (val == "value") {
    label <- "mean"
  } else if (val == "ci_lb") {
    label <- "CI lower bound"
  } else {
    label <- "CI upper bound"
  }
  
  # If cancer patients were considered not at risk after diagnosis in incidence calculate,
  # convert clinical cancer cumulative hazard to cumulative probability,
  # but if they were considered still at risk, the cumulative hazard is actually the cumulative probability
  if (var_censor == "time_H_C") {
    # Plot spline fit
    plot(fit_spline, xlab = "Age", ylab = "Incidence", main = paste("Spline fit:", label),
         xlim = c(0, max(df_incidence$age_end) + 1))
    arrows(df_incidence$age_midpoint, df_incidence$ci_lb, 
           df_incidence$age_midpoint, df_incidence$ci_ub, length = 0.05, angle = 90, code = 3)
    
    # Integrate spline to estimate cumulative hazard of clinical cancer at ages in preclinical cancer data
    chaz_clinical <- sapply(df_prevalence$age_midpoint,
                            function(t) integrate(ir_clinical, lower = 0, upper = t)[["value"]])
    
    # Calculate cumulative probability of clinical cancer
    l_p_clinical[[val]] <- 1 - exp(-chaz_clinical)
  } else if (var_censor == "time_H_D") {
    # Estimate percent dead from cancer
    pct_dead <- sapply(df_prevalence$age_midpoint,
                       function(t) integrate(function(u) ir_clinical(u)*(t - pmax(u, t - 1/d_time_C_D$params$rate))*d_time_C_D$params$rate, lower = 0, upper = t)[["value"]])
    
    # Recalculate probability density of clinical cancer, scaling by percent dead from cancer
    pdf_clinical <- df_incidence[[val]] * (1 - pct_dead)
    
    # Refit spline to clinical cancer probability density
    fit_spline <- cobs(
      x = df_incidence$age_midpoint,
      y = pdf_clinical, 
      constraint = "increase",
      w = 1/(df_incidence$se^2),
      knots = c(0, v_knots, max(df_incidence$age_end) + 1))
    
    # Plot spline fit
    plot(fit_spline, xlab = "Age", ylab = "Probability", main = paste("Spline fit:", label),
         xlim = c(0, max(df_incidence$age_end) + 1))
    arrows(df_incidence$age_midpoint, df_incidence$ci_lb, 
           df_incidence$age_midpoint, df_incidence$ci_ub, length = 0.05, angle = 90, code = 3)
    
    # Calculate cumulative probability of clinical cancer by integrating over density
    l_p_clinical[[val]] <- sapply(df_prevalence$age_midpoint,
                                  function(t) integrate(ir_clinical, lower = 0, upper = t)[["value"]])
  }
  
  # Calculate probabilities for time to preclinical cancer by adding and scaling preclinical prevalence and clinical probability
  df_prevalence[[paste("p", val, sep = "_")]] <- df_prevalence[[val]] * (1 - l_p_clinical[[val]]) + l_p_clinical[[val]]
}

# Calculate true preclinical cancer distribution for plotting
p_preclinical_true <- query_distr("p", seq(0, max(v_ages), 0.5), d_time_H_P$distr, d_time_H_P$params) 

# Plot true clinical cancer CDF
par(mfrow = c(1, 1))
plot(sort(m_times$time_H_C), 1:nrow(m_times)/nrow(m_times), type = "l",  # True clinical cancer
     xlim = c(0, max(v_ages)), ylim = c(0, max(l_p_clinical$ci_ub)),
     xlab = "Age", ylab = "Probability",
     main = "Clinical cancer")

# Plot estimated clinical cancer CDF
for (val in v_cols) {
  if (val == "value") {
    lines(df_prevalence$age_midpoint, l_p_clinical[[val]], col = "red", type = "l")
  } else {
    lines(df_prevalence$age_midpoint, l_p_clinical[[val]], col = "red", type = "l", lty = 2)
  }
}

# Add legend
legend("bottomright", legend=c("True", paste0("Estimate with ", (1 - alpha)*100, "% CI")),
       col = c("black", "red"), lty = c(1, 1))

# Plot true preclinical cancer CDF
plot(seq(0, max(v_ages), 0.5), p_preclinical_true, type = "l",
     xlim = c(0, max(v_ages)), ylim = c(0, max(df_prevalence$p_ci_ub)),
     xlab = "Age", ylab = "Probability",
     main = "Prelinical cancer")

# Plot preclinical cancer prevalence for comparison
lines(df_prevalence$age_midpoint, df_prevalence$value, col = "blue")
points(df_prevalence$age_midpoint, df_prevalence$value, col = "blue")
arrows(df_prevalence$age_midpoint, df_prevalence$ci_lb, 
       df_prevalence$age_midpoint, df_prevalence$ci_ub, length = 0.05, angle = 90, code = 3,
       col = "blue")

# Plot estimated preclinical cancer CDF
for (val in v_cols) {
  if (val == "value") {
    lines(df_prevalence$age_midpoint, df_prevalence[[paste("p", val, sep = "_")]], col = "red", type = "l")
  } else {
    lines(df_prevalence$age_midpoint, df_prevalence[[paste("p", val, sep = "_")]], col = "red", type = "l", lty = 2)
  }
}

# Add legend
legend("bottomright", legend=c("True CDF", paste0("Estimated CDF with ", (1 - alpha)*100, "% CI"), "Prevalence"),
       col = c("black", "red", "blue"), lty = c(1, 1, 1))

# Plot log of true preclinical cancer CDF
plot(seq(0, max(v_ages), 0.5), log(p_preclinical_true), type = "l",
     xlim = c(0, max(v_ages)),
     xlab = "Age", ylab = "Log probability",
     main = "Prelinical cancer")

# Plot log of preclinical cancer prevalence for comparison
lines(df_prevalence$age_midpoint, log(df_prevalence$value), col = "blue")
points(df_prevalence$age_midpoint, log(df_prevalence$value), col = "blue")
arrows(df_prevalence$age_midpoint, log(df_prevalence$ci_lb), 
       df_prevalence$age_midpoint, log(df_prevalence$ci_ub), length = 0.05, angle = 90, code = 3,
       col = "blue")

# Plot estimated preclinical cancer CDF
for (val in v_cols) {
  if (val == "value") {
    lines(df_prevalence$age_midpoint, log(df_prevalence[[paste("p", val, sep = "_")]]), col = "red", type = "l")
  } else {
    lines(df_prevalence$age_midpoint, log(df_prevalence[[paste("p", val, sep = "_")]]), col = "red", type = "l", lty = 2)
  }
}

# Add legend
legend("bottomright", legend=c("True CDF", paste0("Estimated CDF with ", (1 - alpha)*100, "% CI"), "Prevalence"),
       col = c("black", "red", "blue"), lty = c(1, 1, 1))


##### 3.3 Fit Weibull distribution to time to preclinical cancer

# Calculate Weibull x transformation
df_prevalence <- df_prevalence %>%
  mutate(x_transformed = log(age_midpoint))

# Calculate Weibull y transformation
for (val in v_cols) {# Add spline predictions to preclinical cancer prevalence to calculate preclinical cancer cumulative probability
  df_prevalence[[paste("y", val, sep = "_")]] <- log(-log(1 - df_prevalence[[paste("p", val, sep = "_")]]))
}

# Plot transformations to validate Weibull fit - should be linear
for (val in v_cols) {
  if (val == v_cols[1]) {
    plot(df_prevalence$x_transformed, df_prevalence[[paste("y", val, sep = "_")]], 
         xlab = "log(age)", ylab = "log(-log(1 - CDF))",
         main = "Transformations for Weibull regression",
         ylim = range(df_prevalence[paste("y", v_cols, sep = "_")]))
  } else {
    points(df_prevalence$x_transformed, df_prevalence[[paste("y", val, sep = "_")]])
  }
}

# Stack y data
df_prevalence_long <- rbind(
  df_prevalence %>%
    select(x_transformed, y_transformed = y_ci_lb, se) %>%
    mutate(type = "ci_lb"),
  df_prevalence %>%
    select(x_transformed, y_transformed = y_ci_ub, se) %>%
    mutate(type = "ci_ub"))

# Get Weibull estimates using weighted linear regression
fit_lm_mean <- lm(y_value ~ x_transformed, data = df_prevalence, 
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

# Plot true vs. fitted Weibull parameters
curve(pweibull(x, d_time_H_P$params$shape, d_time_H_P$params$scale), from = 0, to = 10,
      xlab = "Time to preclinical cancer", ylab = "Density", main = "Weibull fit to preclinical cancer data")
lines(0:10, pweibull(0:10, shape_H_P, scale_H_P), col = "red")
lines(0:10, pweibull(0:10, shape_H_P_ci[1], scale_H_P_ci[1]), col = "red", lty = 2)
lines(0:10, pweibull(0:10, shape_H_P_ci[2], scale_H_P_ci[2]), col = "red", lty = 2)


##### 3.4 Fit exponential distribution to time from preclinical to clinical cancer

# Create initial exponential distribution estimate for time from preclinical to clinical cancer
# with rate = 1 / age difference between the clinical and preclinical distribution quantiles
# at the maximum probability observed in the clinical distribution
age_P_max_C <- query_distr("q", max(l_p_clinical$value), d_time_H_P_est$distr, d_time_H_P_est$params)
rate_P_C_init <- 1/(max(df_prevalence$age_midpoint) - age_P_max_C)
d_time_P_C_est <- list(distr = "exp", params = list(rate = rate_P_C_init))

# Create objective function for difference between predicted and observed clinical cancer CDF
obj_fn_clinical <- function(rate_P_C, val) {
  # Update distribution for time from preclinical to clinical
  d_time_P_C_est$params$rate <- rate_P_C
  
  # Calculate difference between observed and estimated clinical cancer CDF
  diff_cdf <- l_p_clinical[[val]] - sapply(df_prevalence$age_midpoint, function(t) {
    # Calculates sum of the two distributions
    integrate(function(u)
      query_distr("p", t - u, d_time_H_P_est$distr, d_time_H_P_est$params) * 
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


##### 3.5 Check coverage of targets

# Create dataframe of parameter priors
df_prior_vals <- data.frame(rbind(
  unname(c(shape_H_P, shape_H_P_ci, d_time_H_P$params$shape)),
  unname(c(scale_H_P, scale_H_P_ci, d_time_H_P$params$scale)),
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
    max(l_targets_temp$prevalence$value),
    max(l_targets_temp$incidence$value)))
}

# Plot coverage of targets and visually verify
par(mfrow = c(1, 2))
for (target in c("prevalence", "incidence")) {
  # Plot main estimate
  plot(l_outputs_est[[1]][[target]]$age_midpoint, l_outputs_est[[1]][[target]]$value, type = "l", col = "blue", lwd = 2,
       ylim = c(0, v_target_max[target]),
       xlab = "Age", ylab = toTitleCase(target))
  
  # Plot LHS estimates
  for (i in 1:n_samp) {
    lines(l_outputs_est[[i + 1]][[target]]$age_midpoint, l_outputs_est[[i + 1]][[target]]$value, col = "gray", lty = 2)
  }
  
  # Plot truth on top
  points(l_targets[[target]]$age_midpoint, l_targets[[target]]$value)
  arrows(l_targets[[target]]$age_midpoint, l_targets[[target]]$ci_lb, 
         l_targets[[target]]$age_midpoint, l_targets[[target]]$ci_ub, length = 0.05, angle = 90, code = 3)
}

# Verify that prior bounds cover the truth
df_priors <- df_priors %>%
  mutate(coverage = (truth >= min & truth <= max))
if (sum(df_priors$coverage == FALSE) > 0) {
  stop("Prior bounds do not cover the truth")
} else {
  print("Prior bounds cover the truth")
}
