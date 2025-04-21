###########################  Generate ground truth outputs   ###################
#
#  Objective: Script to generate decision outputs for ground truth
#  parameters; should be run after establishing the model sample size 
#  in analysis/01_load_calibration_params.R, which
########################### <<<<<>>>>> ##############################################

rm(list = ls()) # Clean environment
options(scipen = 999) # View data without scientific notation

#### 1.Libraries and functions  ==================================================

###### 1.1 Load packages
library(readxl)
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
params_model <- configs$params_model
file_params_calib <- configs$paths$file_params_calib
params_screen <- configs$params_screen
params_calib <- configs$params_calib
seed <- params_calib$seed_calib

# Get list of BayCANN output file paths and load to global environment
l_filepaths <- update_config_paths("files_baycann", configs$paths)
list2env(l_filepaths, envir = .GlobalEnv)

###### 2.2 File paths
file_true_params <- file.path("_ground_truth", "true_params.xlsx")

###### 2.3 Other parameters
n_cores_reserved_local <- 2
n_reps <- 100


#### 3. Pre-processing actions  ===========================================

# Load calibration parameters
l_params_calib <- readRDS(file_params_calib)

# Load ground truth model parameters (with file_surv set to NULL as survival data is generated in this script)
l_params_model <- do.call(load_model_params, c(
  modifyList(params_model,
             list(file.surv = NULL),
             keep.null = T),
  list(seed = NULL,
       n_cohort = l_params_calib$l_params_model$n_cohort,
       file.distr = file_true_params)
))

# Map variables to parameters for tuning - make dataframe of all parameters with "src = unknown"
param_map <- make_param_map(l_params_model)

# Set flag to count diagnostic tests in base case scenario
l_params_model$fl_count_tests <- TRUE

# Set base case outcome parameters (include calibration target parameters)
l_params_outcome_base <- params_screen$l_outcome_base

# Set screening test and strategy parameters
l_params_screen <- list(test_chars = params_screen$test_chars,
                        strats = params_screen$strats,
                        surveil = params_screen$surveil)

# Set screening outcome parameters
l_params_outcome_screen <- params_screen$l_outcome_base

# Set counterfactual comparison parameters
l_params_outcome_counter <- params_screen$l_outcome_counterfactual

# Set censor variables
l_censor_vars <- params_calib$l_censor_vars

# Set seed
set.seed(seed, kind = "L'Ecuyer-CMRG")

# Set number of cores to use
registerDoParallel(cores = detectCores(logical = TRUE) - n_cores_reserved_local)


#### 4. Generate BayCANN outputs  ===========================================

# Run model for each input parameter sample and get corresponding targets
# Test sample size with longitudinally vs. cross-sectionally calculated vs. other outcomes
stime <- system.time({
  l_outputs <- foreach(
    i=1:n_reps, 
    .combine=c, 
    .inorder=FALSE, 
    .packages=c("data.table","tidyverse")) %dopar% {
      # Get row of parameters and calculate outputs
      l_calib_outputs <- params_to_outputs(
        l_params_model = l_params_model,
        l_params_outcome = l_params_outcome_base,
        l_params_screen = l_params_screen,
        l_params_outcome_screen = l_params_outcome_screen,
        l_params_outcome_counter = l_params_outcome_counter,
        l_censor_vars = l_censor_vars,
        reshape_output = FALSE
      )
      # Call item to save
      list(l_calib_outputs)
    }
})
print(stime)
closeAllConnections()

# Save model outputs
saveRDS(l_outputs, file = file_outputs)
