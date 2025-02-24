###########################  Generate IMABC outputs   ##########################
#
#  Objective: Script to generate decision outputs for IMABC calibrated
#  parameters
########################### <<<<<>>>>> #########################################

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

# Extract relevant parameters from configs
file_params_calib <- configs$paths$file_params_calib
file_plot_labels <- configs$paths$file_plot_labels

# Get list of relevant output file paths and load to global environment
l_filepaths <- update_config_paths("files_imabc", configs$paths)
list2env(l_filepaths, envir = .GlobalEnv)

# Load coverage analysis parameters from configs file
list2env(configs$params_coverage, envir = .GlobalEnv)


#### 3. Pre-processing actions  ===========================================

# Load model and calibration parameters
l_params_calib <- readRDS(file_params_calib)

# Set flag to count diagnostic tests in base case scenario
l_params_calib$l_params_model$fl_count_tests <- TRUE

# Load IMABC posteriors
l_outputs <- readRDS(file_posterior)
m_params <- l_outputs$good_parm_draws %>%
  dplyr::select(l_params_calib$prior_map$var_id)

# Set decision outcome parameters
l_outcome_params <- params_screen$l_outcome_base

# Set screening test and strategy parameters
l_screen_params <- list(test_chars = params_screen$test_chars,
                        strats = params_screen$strats)

# Set counterfactual comparison parameters
l_outcome_params_counter <- params_screen$l_outcome_counterfactual

# Set seed
set.seed(l_params_calib$l_params_model$seed, kind = "L'Ecuyer-CMRG")

# Set number of cores to use
registerDoParallel(cores = detectCores(logical = TRUE) - l_params_calib$n_cores_reserved_local)


#### 4. Generate IMABC outputs  ===========================================

# Run model for each input parameter sample and get corresponding targets
stime <- system.time({
  m_outputs <- foreach(
    i=1:nrow(m_params), 
    .combine=c, 
    .inorder=TRUE, 
    .packages=c("data.table","tidyverse")) %dopar% {
      # Get row of parameters and calculate outputs
      v_params_update <- m_params[i,]
      l_calib_outputs <- with(l_params_calib, {
        params_to_outputs(
          l_params_model = l_params_model,
          v_params_update = v_params_update,
          param_map = prior_map,
          l_outcome_params = l_outcome_params,
          l_screen_params = l_screen_params,
          l_outcome_params_counter = l_outcome_params_counter,
          l_censor_vars = l_censor_vars,
          reshape_output = FALSE
        )
      })
      # Call item to save
      list(l_calib_outputs)
    }
})
print(stime)
closeAllConnections()

# Save model outputs
saveRDS(list(m_outputs = m_outputs,
             runtime = stime), file = file_outputs)
