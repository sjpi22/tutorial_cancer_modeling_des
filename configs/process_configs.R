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
    for (item_name in names(params_calib$lesion_state_true)) {
      params_calib[[item_name]] <- c(params_calib$lesion_state_true[[item_name]], params_calib[[item_name]])
    }
  }
}

# Remove calibration parameters specific to lesion model
params_calib$lesion_state_true <- NULL