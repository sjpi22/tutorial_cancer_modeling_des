###########################  Process configuration file  ################
#
#  Objective: Process parameters from configuration file
########################### <<<<<>>>>> ##############################################

# Load libraries
library(yaml)

# Read configuration file
configs <- read_yaml(file.path("configs", "configs.yaml")) # Load configuration parameters

# Make updates for lesion model if necessary
if ("lesion_state" %in% names(configs$params_model)) {
  if (configs$params_model$lesion_state == T) {
    for (item_name in names(configs$params_calib$lesion_state_true)) {
      configs$params_calib[[item_name]] <- c(configs$params_calib$lesion_state_true[[item_name]], configs$params_calib[[item_name]])
    }
  }
}

# Remove calibration parameters specific to lesion model
configs$params_calib$lesion_state_true <- NULL

# Define function to update output file paths
update_config_paths <- function(label, paths) {
  # Initialize list of file paths
  l_filepaths <- list()
  
  # Iterate over output types i
  for (i in names(paths[[label]])) {
    # Get root directory of output type
    root_dir <- paths$root_dir[[i]]
    
    # Get prefix for file variable
    if (i %in% c("none", "output")) {
      file_prefix <- "file_"
    } else if (i == "figs") {
      file_prefix <- "file_fig_"
    } else {
      stop("Invalid label")
    }
    
    # Iterate over file names j
    for (j in names(paths[[label]][[i]])) {
      # Create full file path
      if (i == "none") {
        filepath <- paths[[label]][[i]][[j]]
      } else {
        filepath <- file.path(root_dir, paths[[label]][[i]][[j]])
      }
      
      # Append to list of file paths
      l_filepaths[[paste0(file_prefix, j)]] <- filepath
    }
  }
  # Output list of file paths
  return(l_filepaths)
}