###########################  BayCANN outputs   #########################################
#
#  Objective: Script to generate calibration target and decision outputs for 
#  BayCANN calibrated parameters
########################### <<<<<>>>>> ##############################################

rm(list = ls()) # Clean environment
options(scipen = 999) # View data without scientific notation

#### 1.Libraries and functions  ==================================================

###### 1.1 Load packages
library(tidyverse)
library(doParallel)
library(foreach)
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

# Extract calibration parameters from configs
file_params_calib <- configs$paths$file_params_calib

# Get list of BayCANN output file paths and load to global environment
l_filepaths <- update_config_paths("files_baycann", configs$paths)
list2env(l_filepaths, envir = .GlobalEnv)


#### 3. Pre-processing actions  ===========================================

# Load model and calibration parameters
l_params_calib <- readRDS(file_params_calib)

# Load BayCANN calibrated parameters
calibrated_params_baycann <- read.csv(file_posterior) %>%
  dplyr::select(-lp) %>% # Remove last non-parameter column
  as.matrix()

# Set number of cores to use
registerDoParallel(cores = detectCores(logical = TRUE) - l_params_calib$n_cores_reserved_local)

# Unit for LYG 
unit_lyg <- 1000
age_min <- 40

# Run base case model
m_cohort <- run_base_model(l_params_calib$l_params_model)

# Create censor variables
m_cohort$patient_level[, time_screen_censor := pmin(time_H_C, time_H_D, na.rm=T)]

# Calculate base case life years
res_base <- calc_lifeyears(m_cohort$patient_level, 
                           sum_var = "time_H_D",
                           censor_var = "time_screen_censor",
                           age_min = age_min)

# Calculate base case diagnostic tests
ct_tests_base <- m_cohort$patient_level[time_screen_censor > age_min & time_H_C < time_H_D, .(N_diag = .N)]

# Run screening scenario
run_screening_counterfactual(m_cohort, l_params_calib$l_params_model,
                             50, 80, 10, c(L = 0.6, P = 0.9), 0.85,
                             overwrite = TRUE)

# Calculate screening scenario outcomes
res_screen <- calc_lifeyears(m_cohort$patient_level, 
                             sum_var = "time_H_D",
                             censor_var = "time_screen_censor",
                             age_min = age_min)

# Calculate screening and diagnostic tests
ct_tests_screen <- colSums(m_cohort$patient_level[, .SD, .SDcols = patterns("ct_")], na.rm = T)

# Calculate LYG
lyg <- (res_screen$time_total - res_base$time_total) / res_base$N * unit_lyg
print(lyg)


#### 4. Generate BayCANN outputs  ===========================================
# Run model for each input parameter sample and get corresponding targets
stime <- system.time({
  m_outputs <- foreach(
    i=1:nrow(calibrated_params_baycann), 
    .combine=rbind, 
    .inorder=FALSE, 
    .packages=c("data.table","tidyverse")) %dopar% {
      # Get row of parameters and calculate outputs
      v_params_update <- calibrated_params_baycann[i,]
      v_calib_outputs <- with(l_params_calib, {
        params_to_calib_outputs(
          l_params_all = l_params_all,
          v_params_update = v_params_update,
          param_map = prior_map,
          l_outcome_params = l_outcome_params,
          l_censor_vars = l_censor_vars
        )
      })
      # Call item to save
      t(v_calib_outputs)
    }
})
print(stime)
closeAllConnections()

# Save model outputs
saveRDS(m_outputs, file = file_outputs)

