###########################  View screening data generation process  ##########################
#
#  Objective: Script to look at model and screening results
########################### <<<<<>>>>> #########################################

rm(list = ls()) # Clean environment
options(scipen = 999) # View data without scientific notation

#### 1.Libraries and functions  ==================================================

###### 1.1 Load packages
library(tidyverse)
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
params_screen <- configs$params_screen

# Get list of relevant output file paths and load to global environment
l_filepaths <- update_config_paths("files_imabc", configs$paths)
list2env(l_filepaths, envir = .GlobalEnv)

# Load coverage analysis parameters from configs file
list2env(configs$params_coverage, envir = .GlobalEnv)

# Other parameters
strat <- c(1, 7) # Screening strategy/strategies to run
strat <- c(7) # Screening strategy/strategies to run
paramset <- 1 # Parameter set(s) to run


#### 3. Pre-processing actions  ===========================================

# Load model and calibration parameters
l_params_calib <- readRDS(file_params_calib)

# Set flag to count diagnostic tests in base case scenario
l_params_calib$l_params_model$fl_count_tests <- TRUE

# Load IMABC posteriors
l_posteriors <- readRDS(file_posterior)
m_params <- l_posteriors$good_parm_draws %>%
  dplyr::select(l_params_calib$prior_map$var_id)

# Set decision outcome parameters
l_params_calib$l_params_outcome <- params_screen$l_outcome_base

# Set screening test and strategy parameters
l_params_screen <- list(test_chars = params_screen$test_chars,
                        strats = params_screen$strats[strat],
                        surveil = params_screen$surveil)

# Set counterfactual comparison parameters
l_params_outcome_counter <- params_screen$l_outcome_counterfactual

# Set seed
set.seed(l_params_calib$l_params_model$seed)


#### 4. Generate outputs  ===========================================

# Run model for parameter set and get corresponding targets
v_params_update <- m_params[paramset,]
l_calib_outputs <- with(l_params_calib, {
  params_to_outputs(
    l_params_model = l_params_model,
    v_params_update = v_params_update,
    param_map = prior_map,
    l_params_outcome = l_params_outcome,
    l_params_screen = l_params_screen,
    l_params_outcome_counter = l_params_outcome_counter,
    l_censor_vars = l_censor_vars,
    reshape_output = FALSE,
    individual_data = TRUE
  )
})

# Manually examine outputs
