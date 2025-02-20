###########################  Run IMABC  ##########################
#
#  Objective: Program to run IMABC based on vignette at 
#  https://github.com/c-rutter/imabc
#  See here for documentation: https://github.com/c-rutter/imabc/tree/a58a3b7c8db18948ff87fb6be55c6175399f41a2
########################### <<<<<>>>>> #########################################

rm(list = ls()) # Clean environment
options(scipen = 999) # View data without scientific notation

#### 1.Libraries and functions  ==================================================

###### 1.1 Load packages
# Run once
# install.packages("imabc")

library(imabc)
library(MASS)
library(data.table)
library(foreach)
library(parallel)
library(truncnorm)
library(lhs)
library(methods)
library(tidyverse)
library(doParallel)
library(assertthat)

###### 1.2 Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)


#### 2. General parameters ========================================================

###### 2.1 Configurations
# Run file to process configurations
source("configs/process_configs.R")

# Extract relevant parameters from configs
file_params_calib <- configs$paths$file_params_calib

# Get list of relevant output file paths and load to global environment
l_filepaths <- update_config_paths("files_imabc", configs$paths)
list2env(l_filepaths, envir = .GlobalEnv)

# Load IMABC parameters from configs file
list2env(configs$params_imabc, envir = .GlobalEnv)


#### 3. Pre-processing actions  ===========================================

# Load model and calibration parameters
l_params_calib <- readRDS(file_params_calib)

# Define priors
param_df <- with(l_params_calib, {
  data.frame(
    name_var = prior_map$var_id,
    dist_var = prior_map$distr,
    min = prior_map$min,
    max = prior_map$max
  )
})

priors <- as.priors(
  param_df,
  parameter_name = "name_var", dist_base_name = "dist_var"
)

# Calculate current and stopping bounds for targets given alpha values
target_map <- l_params_calib$df_targets %>%
  mutate(current_lower_bounds = targets - se*qnorm(1 - alpha_current/2),
         current_upper_bounds = targets + se*qnorm(1 - alpha_current/2),
         stopping_lower_bounds = ifelse(is.na(ci_lb), targets - se*qnorm(1 - alpha_stop/2), ci_lb),
         stopping_upper_bounds = ifelse(is.na(ci_ub), targets + se*qnorm(1 - alpha_stop/2), ci_ub))

# Define target values
target_df <- target_map %>%
  dplyr::select(target_groups, target_names, targets,
                current_lower_bounds, current_upper_bounds,
                stopping_lower_bounds, stopping_upper_bounds)

targets <- as.targets(target_df)

# Define target function
fn <- function(v_params_update) {
  v_targets <- with(l_params_calib, {
    params_to_calib_outputs(
      l_params_model = l_params_model,
      v_params_update = v_params_update,
      param_map = prior_map,
      l_outcome_params = l_outcome_params,
      l_censor_vars = l_censor_vars
    )
  })
  return(v_targets)
}

target_fun <- define_target_function(
  targets, priors, FUN = fn, use_seed = fn_use_seed
)

# IMABC core parameters
imabc_inputs <- list(
  target_fun = target_fun,
  priors = priors,
  targets = targets,
  N_start = N_start_multiplier * length(priors),
  seed = l_params_calib$l_params_model$seed)

# Add optional args
imabc_inputs <- c(imabc_inputs,
                  optional_args)

# Set number of cores to use
registerDoParallel(cores = detectCores(logical = TRUE) - l_params_calib$n_cores_reserved_local)


#### 4. Run IMABC  ===========================================

# Calibrate model with IMABC
start_time <- Sys.time()
calibration_results <- do.call(imabc, imabc_inputs)
end_time <- Sys.time()
runtime <- end_time - start_time
print(runtime)

saveRDS(calibration_results, file = file_posterior)
write.csv(data.frame(runtime = runtime), file = file_runtime, row.names = F)
