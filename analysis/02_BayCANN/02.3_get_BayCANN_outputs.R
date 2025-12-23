###########################  Generate BayCANN outputs   #########################################
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
# User-set default config
config_default <- "colorectal"

# Set configs to load
if(!is.na(as.integer(Sys.getenv("SLURM_NTASKS_PER_NODE")))) { # If running on cluster
  # Import arguments from Terminal
  args <- commandArgs(trailingOnly = TRUE)
  
  # Set config version to argument value or default if not provided
  if (length(args) > 0) {
    config_version <- args[1]
  } else {
    config_version <- config_default
  }
} else { # If running locally, set configs to default
  config_version <- config_default
}

# Load configs
file_configs <- file.path("configs", paste0("configs_", config_version, ".yaml"))
configs <- load_configs(file_configs)
print(paste("config path:", file_configs))

# Extract relevant parameters from configs
file_params_calib <- configs$paths$file_params_calib
params_screen <- configs$params_screen

# Get list of BayCANN output file paths and load to global environment
l_filepaths <- update_config_paths("files_baycann", configs$paths)
list2env(l_filepaths, envir = .GlobalEnv)


#### 3. Pre-processing actions  ===========================================

# Load model and calibration parameters
l_params_calib <- readRDS(file_params_calib)

# Set flag to count diagnostic tests in base case scenario
l_params_calib$l_params_model$fl_count_tests <- TRUE

# Load BayCANN posteriors
m_params <- read.csv(file_posterior) %>%
  dplyr::select(-lp) %>% # Remove last non-parameter column
  as.matrix()

# Check for internal validation only flag
if (is.null(params_screen$internal_val_only)) {
  internal_val_only <- FALSE
} else {
  internal_val_only <- params_screen$internal_val_only
}

# If only performing internal validation, set screening parameters to NULL
if (internal_val_only) {
  l_params_outcome_base <- l_params_outcome
  l_params_screen <- NULL
  l_params_outcome_screen <- NULL
  reshape_output <- TRUE
} else {
  # Set base case outcome parameters (include calibration target parameters)
  l_params_outcome_base <- c(l_params_calib$l_params_outcome,
                             params_screen$l_outcome_base)
  
  # Set screening test and strategy parameters
  l_params_screen <- list(test_chars = params_screen$test_chars,
                          strats = params_screen$strats,
                          surveil = params_screen$surveil)
  
  # Set screening outcome parameters
  l_params_outcome_screen <- params_screen$l_outcome_base
  
  # Set counterfactual comparison parameters
  l_params_outcome_counter <- params_screen$l_outcome_counterfactual
  
  # Keep output as list
  reshape_output <- FALSE
}

# Set seed
seed <- l_params_calib$l_params_model$seed
set.seed(seed, kind = "L'Ecuyer-CMRG")
l_params_calib$l_params_model$seed <- NULL

# Set number of cores to use
if(!is.na(as.integer(Sys.getenv("SLURM_NTASKS_PER_NODE")))) {
  # If using Sherlock, use environment variable to set the number of cores to use
  registerDoParallel(cores = (Sys.getenv("SLURM_NTASKS_PER_NODE")))
  print("Running on Sherlock")
} else {
  # If running locally, use all available cores except for reserved ones
  registerDoParallel(cores = detectCores(logical = TRUE) - l_params_calib$n_cores_reserved_local)
  print("Running locally")
}

# Show the number of parallel workers to be used
print(paste("# parallel workers:", getDoParWorkers())) 


#### 4. Generate BayCANN outputs  ===========================================

# Run model for each input parameter sample and get corresponding targets
stime <- system.time({
  l_outputs <- foreach(
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
          l_params_outcome = l_params_outcome_base,
          l_params_screen = l_params_screen,
          l_params_outcome_screen = l_params_outcome_screen,
          l_params_outcome_counter = l_params_outcome_counter,
          l_censor_vars = l_censor_vars,
          reshape_output = reshape_output
        )
      })
      # Call item to save
      list(l_calib_outputs)
    }
})
print(stime)
closeAllConnections()

# Save model outputs
saveRDS(list(l_outputs = l_outputs,
             runtime = stime), file = file_outputs)
