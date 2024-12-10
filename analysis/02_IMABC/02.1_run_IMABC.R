# Run IMABC: Vignette at https://github.com/c-rutter/imabc

# Run once
# install.packages("imabc")

# Load packages
library(imabc)
library(MASS)
library(data.table)
library(foreach)
library(parallel)
library(truncnorm)
library(lhs)
library(methods)
library(stats)
library(utils)
library(tidyverse)
library(doParallel)
library(assertthat)

#* Clean environment
rm(list = ls())

# Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)

#### 2. General parameters ========================================================

###### 2.1 file paths
file_params <- "data/calibration_params.rds"
outpath <- "output/calibration/IMABC"
file_imabc_params <- file.path(outpath, "params_IMABC.rds")
file_posterior <- file.path(outpath, "calibrated_posteriors_IMABC.rds")

###### 2.2 IMABC parameters 
alpha_current <- 1e-15
alpha_stop <- 0.05
fn_use_seed <- FALSE
N_start_multiplier = 1000
optional_args = list(
  N_centers = 1,
  Center_n = 100,
  N_post = 2000,
  max_iter = 1000
)


#### 3. Pre-processing actions  ===========================================

# Load model and calibration parameters
l_params_calib <- readRDS(file_params)

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
target_map <- l_params_calib$df_true_targets %>%
  mutate(current_lower_bounds = targets - se*qnorm(1-alpha_current/2),
         current_upper_bounds = targets + se*qnorm(1-alpha_current/2),
         stopping_lower_bounds = targets - se*qnorm(1-alpha_stop/2),
         stopping_upper_bounds = targets + se*qnorm(1-alpha_stop/2))

# Define targets values
target_df <- target_map %>%
  dplyr::select(target_groups, target_names, targets,
                current_lower_bounds, current_upper_bounds,
                stopping_lower_bounds, stopping_upper_bounds)

targets <- as.targets(target_df)

# Define targets function
fn <- function(v_params_update) {
  v_targets <- with(l_params_calib, {
    params_to_calib_outputs(
      l_params_all = l_params_all,
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
  seed = l_params_calib$l_params_all$seed)

# Add optional args
imabc_inputs <- c(imabc_inputs,
                  optional_args)

# Set number of cores to use
registerDoParallel(cores=detectCores(logical = TRUE) - l_params_calib$n_cores_reserved_local)


#### 4. Run IMABC  ===========================================

# Calibrate model - see here for documentation: https://github.com/c-rutter/imabc/tree/a58a3b7c8db18948ff87fb6be55c6175399f41a2
start_time <- Sys.time()
calibration_results <- do.call(imabc, imabc_inputs)
end_time <- Sys.time()
print(end_time - start_time)

print('Saving output')
saveRDS(calibration_results, file = file_posterior)
saveRDS(list(imabc_inputs = imabc_inputs,
             runtime = end_time - start_time), file = file_imabc_params)

