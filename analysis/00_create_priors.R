###########################  Create prior distributions  ##########################
#
#  Objective: Refine data-driven prior distributions for model parameters
########################### <<<<<>>>>> #########################################

rm(list = ls()) # Clean environment
options(scipen = 999) # View data without scientific notation

#### 1.Libraries and functions  ==================================================

###### 1.1 Load packages
library(readxl)
library(tools)
library(tidyverse)
library(assertthat)
library(data.table)
library(cobs) # For fitting constrained B-splines
library(sandwich) # For heteroskedasticity-robust linear model standard errors
library(CVXR) # For deconvolution
library(flexsurv) # If using log-logistic

###### 1.2 Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)



#### 2. General parameters ========================================================

###### 2.1 Configurations
# Load configs
file_configs <- file.path("configs", "configs_bladder.yaml")
configs <- load_configs(file_configs)

# Extract relevant parameters from configs
params_model <- configs$params_model
params_priors <- configs$params_priors
params_calib <- configs$params_calib
file_targets <- configs$params_calib$file_targets

# Load prior parameters to global environment
list2env(params_priors, envir = .GlobalEnv)


#### 3. Load data  ===========================================

# Set seed
set.seed(params_calib$seed_calib)

# Load model parameters
l_params_model <- do.call(load_model_params, c(
  params_model,
  list(seed = params_calib$seed_calib)
))

# Load targets
df_targets <- load_calibration_targets(file_targets)

# Read prior data frame
df_priors <- read.csv(params_calib$file_prior)

# Extract relevant targets and create containers for CDFs
l_targets <- list()
x_states <- list()
cdf_states <- list()
v_target_indices <- c()
for (state in names(v_state_targets)) {
  if (!is.null(v_state_targets[[state]])) {
    # Select targets
    l_targets[[state]] <- df_targets %>%
      filter(target_groups == v_state_targets[[state]])
    
    # Append target indices
    v_target_indices <- c(v_target_indices, unique(l_targets[[state]][[var_index]]))
    
    # Process incidence data
    if (params_calib$l_params_outcome[[v_state_targets[[state]]]][["outcome_type"]] == "incidence") {
      # Rescale incidence values by unit
      for (val in c(v_cols)) {
        l_targets[[state]][[val]] <- l_targets[[state]][[val]] / params_calib$l_params_outcome[[v_state_targets[[state]]]]$lit_params$rate_unit
      }
    }
    
    # Create container for CDF and x-values
    x_states[[state]] <- list()
    cdf_states[[state]] <- list()
  }
}

# De-duplicate target indices
v_target_indices <- sort(unique(v_target_indices))

# Set variables dependent on parameters
max_age <- l_params_model$max_age
v_ages <- seq(0, max_age)
n_states <- length(cdf_states)

# Get variable for disease onset
var_onset <- paste0("time_H_", l_params_model$v_states[2])

# Get variable for censoring incidence
var_censor <- params_calib$l_params_outcome$incidence$lit_params$censor_var


#### 4. Derive prior distribution for time to disease onset  ===========================================

##### 4.1 Derive time-to-event distributions from targets using splines
# Loop over targets, lower bounds, and upper bounds
x_states[[1]] <- v_target_indices
for (val in v_cols) {
  cdf_states[[1]][[val]] <- with(l_targets[[1]], {
    # Calculate CDF of clinical cancer
    incidence_rate_to_cdf(
      x_vals = get(var_index),
      y_vals = get(val),
      x_pred = x_states[[1]],
      censored_at_event = ifelse(var_censor == "time_H_C", TRUE, FALSE),
      constraints = "increase"
    )
  })
  
  # Calculate probabilities for time to preclinical cancer by adding and scaling preclinical prevalence and clinical probability
  x_states[[2]] <- v_target_indices[v_target_indices %in% l_targets[[2]][[var_index]]]
  cdf_clin <- cdf_states[[1]][[val]][v_target_indices %in% l_targets[[2]][[var_index]]] # Subset to preclinical target indices
  if (!"condition_var" %in% names(params_calib$l_params_outcome[[v_state_targets[[2]]]][["lit_params"]])) {
    # If preclinical target is prevalence out of total population, use directly
    cdf_states[[2]][[val]] <- l_targets[[2]][[val]] * (1 - cdf_clin) + cdf_clin
  } else {
    # Otherwise if preclinical target is conditional on cancer onset before death
    # (e.g., proportion of cancer cases that are incidental), CDF of preclinical
    # times proportion of known cases equals CDF of clinical
    cdf_states[[2]][[val]] <- cdf_clin / (1 - l_targets[[2]][[val]])
  }
  
  # Calculate probabilities for time to precancerous lesion by adding and scaling lesion prevalence and preclinical probability
  # Note - multiplied by clinical CDF rather than preclinical CDF because lesion screening study includes preclinical cases in the denominator
  if (params_model$lesion_state == T) {
    # If indices are same for preclinical and lesion targets, use directly
    if (all.equal(l_targets[[2]][[var_index]], l_targets[[3]][[var_index]]) == TRUE) {
      cdf_preclin <- cdf_states[[2]][[val]]
    } else {
      # Otherwise, fit spline for preclinical CDF
      cdf_preclin <- with(l_targets[[2]], {
        # Fit spline
        spline_cdf <- fit_spline(
          x_vals = get(var_index),
          y_vals = get(val),
          constraints = "increase"
        )
        
        # Predict on lesion target indices
        cdf_preclin <- predict(spline_cdf, l_targets[[3]][[var_index]])[, "fit"]
      })
    }
    
    # Calculate CDF of lesion onset
    x_states[[3]] <- v_target_indices[v_target_indices %in% l_targets[[3]][[var_index]]]
    cdf_clin <- cdf_states[[1]][[val]][v_target_indices %in% l_targets[[3]][[var_index]]] # Subset to lesion target indices
    cdf_states[[3]][[val]] <- l_targets[[3]][[val]] * (1 - cdf_clin) + cdf_preclin
  }
}

# Find maximum upper bound for plotting
max_ub <- max(sapply(cdf_states, function(x) max(x[[v_cols[3]]])))

# Plot fitted CDFs from targets
for (i in seq(cdf_states)) {
  # Subset to state CDF and x-values for plotting
  state <- cdf_states[[i]]
  
  # Plot estimated CDF
  if (i == 1) {
    plot(x_states[[i]], state[[v_cols[1]]], col = v_colors[i], ylim = c(0, max_ub),
         xlab = "Age", ylab = "Probability", main = "CDFs of time to event")
  } else {
    points(x_states[[i]], state[[v_cols[1]]], col = v_colors[i])
  }
  
  # Plot lower and upper bounds
  arrows(x_states[[i]], state[[v_cols[2]]], 
         x_states[[i]], state[[v_cols[3]]], 
         col = v_colors[i], 
         length = 0.05, angle = 90, code = 3)
}


##### 4.2 Fit distribution for time to onset

# Get distribution for onset variable
distr_onset <- unique(df_priors[df_priors$var_name == paste0("d_", var_onset), "var_distr"])

# Set max for cure model
if ("cure_max" %in% df_priors[df_priors$var_name == paste0("d_", var_onset), "param_name"]) {
  cure_max <- df_priors[df_priors$var_name == paste0("d_", var_onset) & df_priors$param_name == "cure_max", "param_val"]
} else {
  cure_max <- 1
}

# Loop over disease states
params_onset <- list()
for (val in v_cols) {
  # Fit linear model to transformed CDF
  params_onset[[val]] <- with(l_targets[[n_states]], {
    do.call(paste0("fit_", distr_onset), list(x_states[[n_states]], pmin(1, cdf_states[[n_states]][[val]] / cure_max)))
  })
}

# Get min and max parameters for priors 
params_onset_reshaped <- list(
  shape = range(sapply(params_onset, function(x) x$shape)),
  scale = range(sapply(params_onset, function(x) x$scale))
)

# Multiply by bounds multiplier
params_onset_scaled <- lapply(params_onset_reshaped, 
                              function(x) x*c(1 - multiplier_bounds, 1 + multiplier_bounds))

# Plot distribution fit for onset
for (val in v_cols) {
  lines(v_ages, 
        query_distr("p", v_ages, distr = distr_onset, 
                    params = list(shape = params_onset[[val]]$shape, 
                                  scale = params_onset[[val]]$scale, 
                                  cure_max = cure_max)), 
        lty = ifelse(val == v_cols[1], 1, 2),
        col = v_colors[n_states])
}


##### 4.3 Fit distributions for times between subsequent states
# Fit distributions to time from lesion to cancer onset
if (params_model$lesion_state == T) {
  # Set index of cancer onset state from v_state_targets
  state_P <- 2
  
  # Loop over targets, lower, and upper bounds
  params_L_P <- list()
  for (val in v_cols) {
    # Set bound to take for previous state (same for target, opposite for lower and upper bound)
    if (val == v_cols[1]) {
      val_prev = v_cols[1]
    } else if (val == v_cols[2]) {
      val_prev = v_cols[3]
    } else if (val == v_cols[3]) {
      val_prev = v_cols[2]
    }
    
    # Perform deconvolution using CDFs of lesion and cancer onset to get PDF of time between states
    pdf_L_P <- deconvolve(x_target = x_states[[state_P]],
                          y_target = cdf_states[[state_P]][[val]],
                          fn_cdf_1 = function(x) query_distr("p", x, distr = distr_onset, params = list(shape = params_onset[[1]]$shape, scale = params_onset[[1]]$scale, cure_max = cure_max)),
                          delta = delta,
                          penalty_jump = penalty1,
                          penalty_flip = penalty2)
    
    # Fit distribution
    params_L_P[[val]] = do.call(paste0("fit_", distr_onset), list(x_vals = pdf_L_P$x, y_vals = cumsum(pdf_L_P$pdf)))
  }
  
  # Get min and max parameters for priors 
  params_L_P_reshaped <- list(
    shape = range(sapply(params_L_P, function(x) x$shape)),
    scale = range(sapply(params_L_P, function(x) x$scale))
  )
  
  # Multiply by bounds multiplier
  params_L_P_scaled <- lapply(params_L_P_reshaped, 
                              function(x) x*c(1 - multiplier_bounds, 1 + multiplier_bounds))
} else {
  # Perform deconvolution using CDFs of lesion and cancer onset to get PDF of time between states
  pdf_P_C <- deconvolve(x_target = x_states[["C"]],
                        y_target = cdf_states[["C"]][[1]],
                        fn_cdf_1 = function(x) query_distr("p", x, distr = distr_onset, params = list(shape = params_onset[[1]]$shape, scale = params_onset[[1]]$scale, cure_max = cure_max)),
                        delta = delta,
                        penalty_jump = penalty1,
                        penalty_flip = penalty2)
  plot(pdf_P_C$x, pdf_P_C$pdf)
}


##### 4.4 Update prior distribution
# Update prior dataframe
for (param in names(params_onset_scaled)) {
  df_priors[df_priors$var_id == paste0("d_", var_onset, ".", param), c("min", "max")] <- as.list(params_onset_scaled[[param]])
  
  # Add mean if there is an initial guess column
  if ("param_val" %in% colnames(df_priors)) {
    df_priors[df_priors$var_id == paste0("d_", var_onset, ".", param), "param_val"] <- params_onset[[1]][[param]]
  }
}

# Overwrite prior file
write.csv(df_priors, file = params_calib$file_priors, row.names = FALSE)
