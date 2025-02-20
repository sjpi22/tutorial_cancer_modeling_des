###########################  Generate BayCANN Sample  ##########################
#
#  Objective: Program to simulate parameter inputs and model outputs for 
#  BayCANN model calibration
########################### <<<<<>>>>> #########################################

rm(list = ls()) # Clean environment
options(scipen = 999) # View data without scientific notation

#### 1.Libraries and functions  ==================================================

###### 1.1 Load packages
library(data.table)
library(tidyverse)
library(lhs)
library(doParallel)
library(foreach)
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
l_filepaths <- update_config_paths("files_baycann", configs$paths)
list2env(l_filepaths, envir = .GlobalEnv)

# Load BayCANN parameters from configs file
list2env(configs$params_baycann$params_sampling, envir = .GlobalEnv)


#### 3. Pre-processing  ===========================================

# Load model and calibration parameters
l_params_calib <- readRDS(file_params_calib)

# Set seed
set.seed(l_params_calib$seed_calib)

# Set number of cores to use
registerDoParallel(cores = detectCores(logical = TRUE) - l_params_calib$n_cores_reserved_local)


#### 4. Generate random set of inputs  ===========================================

m_param_samp <- with(l_params_calib, {
  # Get number of params to calibrate and number of samples
  n_param <- nrow(prior_map)
  n_samp <- n_samp_per_param * n_param
    
  # Sample unit Latin Hypercube
  m_lhs_unit <- randomLHS(n_samp, n_param)
  
  # Rescale to min/max of each parameter
  m_param_samp <- matrix(nrow = n_samp, ncol = n_param)
  for (i in 1:n_param) {
    m_param_samp[, i] <- qunif(
      m_lhs_unit[, i],
      min = prior_map$min[i],
      max = prior_map$max[i])
  }
  colnames(m_param_samp) <- prior_map$var_id
  
  return(m_param_samp)
}
)


#### 5. Generate corresponding outputs  ===========================================

# Run model for each input parameter sample and get corresponding outputs with parallel processing
stime <- system.time({
  m_outputs <- foreach(
    i=1:nrow(m_param_samp), 
    .combine=rbind, 
    .inorder=TRUE, 
    .packages=c("data.table","tidyverse")) %dopar% {
      # Get row of parameters and calculate outputs
      v_params_update <- m_param_samp[i,]
      v_calib_outputs <- with(l_params_calib, {
        params_to_calib_outputs(
          l_params_model = l_params_model,
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
  
# Set column names
colnames(m_outputs) <- l_params_calib$df_target$target_names

# Check for any NaN
validate_that(
  sum(sapply(m_param_samp, function(x) any(is.nan(x)))) == 0, 
  msg = "Parameters include NaN")
validate_that(
  sum(sapply(m_outputs, function(x) any(is.nan(x)))) == 0, 
  msg = "Outputs include NaN")

# Save parameter sample and corresponding outputs
saveRDS(list(m_param_samp = m_param_samp, 
             m_calib_outputs = m_outputs,
             runtime = stime), 
        file = file_sample)
