###########################  Create Prior Distributions  ##########################
#
#  Objective: Establish data-driven prior distributions for model parameters
#  if priors do not exist
########################### <<<<<>>>>> #########################################


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

# Define nonnegative smooth function for clinical cancer incidence (note: fit_spline must be a global variable)
ir_clinical <- function(x) {
  return(pmax(0, predict(fit_spline, x)[, "fit"]))
}

#### 2. General parameters ========================================================

###### 2.1 file paths 
target_files <- list(
  prevalence = "data/prevalence_asymptomatic_cancer.csv",
  incidence = "data/incidence_symptomatic_cancer.csv",
  stage_distr = "data/stage_distr.csv")
prior_file <- "data/priors.rds"

###### 2.2 model parameters
seed_calib <- 42
max_age <- 111
n_cohort <- 100000
v_ages <- seq(0, max_age, 0.25)
alpha <- 0.01
multiplier_bounds <- 0.2

# Define parameters for calculating outcomes
# Names of lists should be the same as those of target_files
l_outcome_params <- list(
  prevalence = list(
    outcome_name = "prevalence",
    m_patients = "m_patients",
    start_var = "time_H_P", 
    end_var = "time_H_C", 
    censor_var = "time_screen_censor"),
  incidence = list(
    outcome_name = "incidence",
    m_patients = "m_patients",
    time_var = "time_H_C", 
    censor_var = "time_H_D"),
  stage_distr = list(
    outcome_name = "distr",
    m_patients = "m_patients",
    grouping_var = "stage_dx", 
    event_var = "time_H_C", 
    censor_var = "time_H_D",
    groups_expected = 1:4)
)

l_censor_vars <- list(
  time_screen_censor = c("time_H_C", "time_H_D")
)

v_cols <- c("targets", "ci_lb", "ci_ub")


#### 3. Load data  ===========================================

# Load targets
l_true_targets <- load_calibration_targets(target_files)

# Calculate approximate confidence intervals for incidence
l_true_targets$incidence <- l_true_targets$incidence %>%    
  mutate(ci_lb = pmax(0, targets - qnorm(1 - alpha/2)*se),
         ci_ub = targets + qnorm(1 - alpha/2)*se)

# Rescale incidence values by unit
for (val in c(v_cols, "se")) {
  l_true_targets$incidence[[val]] <- l_true_targets$incidence[[val]] / l_true_targets$incidence$unit
}

# Load survival data
df_surv <- load_surv_data("data/relative_survival_cancer.csv")


#### 4. Derive prior distribution for time to disease onset  ===========================================

##### 4.1 Derive average time from diagnosis to death by fitting splines

# Initialize simulation population 
m_death <- data.table()

# Sample every three values of ages for spline knots
v_knots_Dc <- df_surv[df_surv$stage == 1, "years_from_dx"]
v_knots_Dc <- c(v_knots_Dc[seq(1, length(v_knots_Dc), 4)], max_age)

# Fit splines to survival
v_cancer <- unique(df_surv$stage)
par(mfrow = c(length(v_cancer)/2, 2))
d_time_C_Dc <- list()
for (i in v_cancer) {
  # Filter survival data to stage at diagnosis
  temp_surv_data <- df_surv[df_surv$stage == i, ]
  
  # Fit spline to survival data
  fit_spline_Dc <- cobs(
    x = temp_surv_data$years_from_dx,
    y = temp_surv_data$pct_died, 
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
  
  # Create distribution data
  d_time_C_Dc[[i]] <- list(distr = "empirical", 
                           params = list(xs = v_ages, 
                                         probs = probs, 
                                         max_x = max_age))
  
  # Simulate time to death
  m_death[, paste0("time_C", i, "_Dc") := query_distr(
    "r", n_cohort, 
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
  if (l_outcome_params$incidence$censor_var == "time_H_C") {
    # Plot spline fit
    plot(fit_spline, xlab = "Age", ylab = "Incidence", main = paste("Spline fit:", label),
         xlim = c(0, max(l_true_targets$incidence$age_end) + 1))
    
    # Integrate spline to estimate cumulative hazard of clinical cancer at ages in preclinical cancer data
    chaz_clinical <- sapply(l_true_targets$prevalence$target_index,
                            function(t) integrate(ir_clinical, lower = 0, upper = t)[["value"]])
    
    # Calculate cumulative probability of clinical cancer
    l_p_clinical[[val]] <- 1 - exp(-chaz_clinical)
  } else if (l_outcome_params$incidence$censor_var == "time_H_D") {
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
}

# Plots
par(mfrow = c(1, 1))
for (val in v_cols) {
  if (val == "targets") {
    # Plot estimated preclinical cancer CDF
    plot(l_true_targets$prevalence$target_index, l_true_targets$prevalence$p_targets,
         type = "l",
         xlab = "Age", ylab = "Probability")
    
    # Plot estimated clinical cancer CDF
    lines(l_true_targets$prevalence$target_index, l_p_clinical[[val]], col = "red", type = "l")
  } else {
    lines(l_true_targets$prevalence$target_index, l_p_clinical[[val]], col = "red", type = "l", lty = 2)
    lines(l_true_targets$prevalence$target_index, l_true_targets$prevalence[[paste("p", val, sep = "_")]], col = "black", type = "l", lty = 2)
  }
}

# Plot preclinical cancer prevalence for comparison
arrows(l_true_targets$prevalence$target_index, l_true_targets$prevalence$ci_lb, 
       l_true_targets$prevalence$target_index, l_true_targets$prevalence$ci_ub, length = 0.05, angle = 90, code = 3,
       col = "blue")

# Add legend
legend("bottomright", legend=c("Clinical CDF", "Preclinical CDF", "Prevalence"),
       col = c("black", "red", "blue"), lty = c(1, 1))


##### 4.3 Fit Weibull distribution to time to preclinical cancer

# Calculate Weibull x transformation
l_true_targets$prevalence <- l_true_targets$prevalence %>%
  mutate(x_transformed = log(target_index))

# Calculate Weibull y transformation
for (val in v_cols) {# Add spline predictions to preclinical cancer prevalence to calculate preclinical cancer cumulative probability
  l_true_targets$prevalence[[paste("y", val, sep = "_")]] <- log(-log(1 - l_true_targets$prevalence[[paste("p", val, sep = "_")]]))
}

# Plot transformations to validate Weibull fit - should be linear
for (val in v_cols) {
  if (val == v_cols[1]) {
    plot(l_true_targets$prevalence$x_transformed, l_true_targets$prevalence[[paste("y", val, sep = "_")]], 
         xlab = "log(age)", ylab = "log(-log(1 - CDF))",
         main = "Transformations for Weibull regression",
         ylim = range(l_true_targets$prevalence[paste("y", v_cols, sep = "_")]))
  } else {
    points(l_true_targets$prevalence$x_transformed, l_true_targets$prevalence[[paste("y", val, sep = "_")]])
  }
}

# Stack y data
l_true_targets$prevalence_long <- rbind(
  l_true_targets$prevalence %>%
    select(x_transformed, y_transformed = y_ci_lb, se) %>%
    mutate(type = "ci_lb"),
  l_true_targets$prevalence %>%
    select(x_transformed, y_transformed = y_ci_ub, se) %>%
    mutate(type = "ci_ub"))

# Get Weibull estimates using weighted linear regression
fit_lm_mean <- lm(y_targets ~ x_transformed, data = l_true_targets$prevalence, 
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
df_priors <- readRDS(prior_file)
df_priors[df_priors$var_id == "d_time_H_P.shape", c("min", "max")] <- as.list(shape_H_P_bounds)
df_priors[df_priors$var_id == "d_time_H_P.scale", c("min", "max")] <- as.list(scale_H_P_bounds)

saveRDS(df_priors, prior_file)
