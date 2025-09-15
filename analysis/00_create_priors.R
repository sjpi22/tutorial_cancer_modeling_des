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

###### 1.2 Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)



#### 2. General parameters ========================================================

###### 2.1 Configurations
# Load configs
file_configs <- file.path("configs", "configs_colorectal.yaml")
configs <- load_configs(file_configs)

# Extract relevant parameters from configs
params_model <- configs$params_model
params_calib <- configs$params_calib
file_targets <- configs$params_calib$file_targets

###### 2.2 Other parameters
# Health states working backwards from latest to earliest and associated target group names
v_state_targets <- c(C = "incidence_clin",
                     P = "prevalence_preclin",
                     L = "prevalence_lesion")
conf_level <- 0.95 # For generating bounds
multiplier_bounds <- 0.2
v_cols <- c("targets", "ci_lb", "ci_ub")
v_colors <- c("green", "red", "orange")
var_index <- "target_index"

# Parameters for deconvolution
delta <- 1
penalty1 <- 0 # Penalty for large changes (spikes) in discrete values of PDF
penalty2 <- 1 # Penalty for changes in sign in PDF


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

# Extract relevant targets and create containers for CDFs
l_targets <- list()
x_states <- list()
cdf_states <- list()
v_target_indices <- c()
for (state in names(v_state_targets)) {
  if (!is.null(v_state_targets[state])) {
    # Select targets
    l_targets[[state]] <- df_targets %>%
      filter(target_groups == v_state_targets[state])
    
    # Append target indices
    v_target_indices <- c(v_target_indices, unique(l_targets[[state]][[var_index]]))
    
    # Process incidence data
    if (params_calib$l_params_outcome[[v_state_targets[state]]][["outcome_type"]] == "incidence") {
      # Rescale incidence values by unit
      for (val in c(v_cols)) {
        l_targets[[state]][[val]] <- l_targets[[state]][[val]] / params_calib$l_params_outcome[[v_state_targets[state]]]$lit_params$rate_unit
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
for (val in v_cols) {
  x_states[[1]] <- v_target_indices
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
  cdf_states[[2]][[val]] <- l_targets[[2]][[val]] * (1 - cdf_clin) + cdf_clin
  
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


##### 4.2 Fit Weibull distribution for time to onset

# Loop over disease states
params_onset <- list()
for (val in v_cols) {
  # Fit Weibull linear model
  params_onset[[val]] <- with(l_targets[[n_states]], {
    fit_weibull(x_states[[n_states]],
                cdf_states[[n_states]][[val]])
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

# Plot fitted CDFs from targets
for (i in seq(cdf_states)) {
  # Subset to state CDF and x-values for plotting
  state <- cdf_states[[i]]

  # Plot estimated CDF
  if (i == 1) {
    plot(x_states[[i]], state[[v_cols[1]]], ylim = c(0, 1), col = v_colors[i],
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

# Plot Weibull fit for onset
for (val in v_cols) {
  lines(v_ages, 
        pweibull(v_ages, 
                 shape = params_onset[[val]]$shape, 
                 scale = params_onset[[val]]$scale), 
        lty = ifelse(val == v_cols[1], 1, 2),
        col = v_colors[n_states])
}

##### 4.3 Update prior distribution
# Update prior dataframe
df_priors <- read.csv(params_calib$file_prior)
for (param in names(params_onset_scaled)) {
  df_priors[df_priors$var_id == paste0("d_time_H_", l_params_model$v_states[2], ".", param), c("min", "max")] <- as.list(params_onset_scaled[[param]])
}

# Overwrite prior file
write.csv(df_priors, file = params_calib$file_priors, row.names = FALSE)
