###########################  Process configuration file  ################
#
#  Objective: Process parameters from configuration file
########################### <<<<<>>>>> ##############################################

# Load libraries
library(yaml)

# Read configuration file
configs <- read_yaml(file.path("configs", "configs.yaml")) # Load configuration parameters
list2env(configs, envir = .GlobalEnv) # Assign items within configs to global environment

# Make updates for lesion model if necessary
if ("lesion_state" %in% names(params_model)) {
  if (params_model$lesion_state == T) {
    params_calib$l_outcome_params <- c(params_calib$l_outcome_params_lesion, params_calib$l_outcome_params)
  }
}