###########################  Load Calibration Parameters  ##########################
#
#  Objective: Program to load general calibration parameters
########################### <<<<<>>>>> #########################################

rm(list = ls()) # Clean environment
options(scipen = 999) # View data without scientific notation

#### 1.Libraries and functions  ==================================================

###### 1.1 Load packages
library(readxl)
library(tidyverse)
library(lhs)
library(assertthat)
library(data.table)

###### 1.2 Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)


#### 2. General parameters ========================================================

###### 2.1 Configurations
# Run file to process configurations
source("configs/process_configs.R")

###### 2.1 file paths
outpath <- paths$output
file_params <- file.path(outpath, "calibration", "calibration_params.rds")

###### 2.3 Monte Carlo analysis parameters
n_init <- 10000
n_mc_reps <- 50
mc_multiplier <- 2


#### 3. Load parameters  ===========================================

# Load model parameters
l_params_init <- do.call(load_model_params, c(
  params_model,
  list(seed = params_calib$seed_calib)
))

# Load calibration parameters
l_params_calib <- do.call(load_calib_params, c(
  l_params_model = list(l_params_init),
  params_calib
))


#### 4. Monte Carlo error analysis  ===========================================

# Set Monte Carlo simulation parameters
l_params_mc <- l_params_calib$l_params_model
l_params_mc$n_cohort <- n_init
l_params_mc$seed <- NULL # Remove seed that is reset every time model is run
set.seed(l_params_calib$l_params_model$seed) # Set seed externally

df_res_mc <- matrix(nrow = n_mc_reps, ncol = nrow(l_params_calib$df_true_targets))
for (i in 1:n_mc_reps) {
  df_res_mc[i, ] <- params_to_calib_outputs(l_params_mc, 
                                            l_outcome_params = l_params_calib$l_outcome_params,
                                            l_censor_vars = l_params_calib$l_censor_vars)
}

# Get column-wise mean and SD
mean_mc <- colMeans(df_res_mc)
sd_mc <- apply(df_res_mc, 2, sd)

# Calculate required sample size
n_target <- n_init * (sd_mc / l_params_calib$df_true_targets$se)^2

# Set final N as multiple of maximum required sample size, rounded up to second highest digit
n_final <- mc_multiplier * max(n_target)
n_digits <- floor(log(n_final, base = 10))
n_final <- ceiling(n_final / 10^(n_digits - 1)) * 10^(n_digits - 1)

# Update parameters
l_params_calib$l_params_model$n_cohort <- n_final


#### 5. Save parameters  ===========================================

saveRDS(l_params_calib, file = file_params)
