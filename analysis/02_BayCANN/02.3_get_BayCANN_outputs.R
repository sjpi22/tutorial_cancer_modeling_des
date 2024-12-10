###########################  BayCANN outputs   #########################################
#
#  Objective: Script to generate calibration target outputs for BayCANN calibrated parameters
########################### <<<<<>>>>> ##############################################


#### 1.Libraries and functions  ==================================================
#* Clean environment
rm(list = ls())

library(tidyverse)
library(doParallel)
library(foreach)
library(data.table)

###### 1.1 Load functions =================================================
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)


#### 2. General parameters ========================================================

###### 2.1 file paths
file_params <- "data/calibration_params.rds"
outpath <- "output/calibration/BayCANN"
file_posterior <- file.path(outpath, "calibrated_posteriors_BayCANN.csv")
file_outputs <- file.path(outpath, "calibration_outputs_BayCANN.rds")


#### 3. Pre-processing actions  ===========================================

# Load model and calibration parameters
l_params_calib <- readRDS(file_params)

# Load BayCANN calibrated parameters
calibrated_params_baycann <- read_csv(file_posterior) %>%
  dplyr::select(-lp) %>% # Remove last non-parameter column
  as.matrix()

# Set number of cores to use
registerDoParallel(cores=detectCores(logical = TRUE) - l_params_calib$n_cores_reserved_local)


################################################################################
###  Generate BayCANN outputs                                          ###
################################################################################
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

