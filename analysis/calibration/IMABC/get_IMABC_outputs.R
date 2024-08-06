###########################  IMABC outputs   #########################################
#
#  Objective: Script to generate calibration target outputs for IMABC calibrated parameters
########################### <<<<<>>>>> ##############################################


#### 1.Libraries and functions  ==================================================
library(tidyverse)
library(dplyr)
library(doParallel)
library(foreach)

###### 1.1 Load functions =================================================
#* Clean environment
rm(list = ls())

# Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)


#### 2. General parameters ========================================================

###### 2.1 file paths 
prior_file <- "data/priors.RData"

targets_files <- list(prevalence = "data/prevalence_asymptomatic_cancer.csv",
                      incidence = "data/incidence_symptomatic_cancer.csv",
                      stage_distr = "data/stage_distr.csv")

path_baycann <- "output/calibrated_posteriors_BayCANN.csv"
outpath <- "output/calibration_targets_BayCANN.csv"

###### 2.1 model parameters 

# Load default data
l_params_all <- load_default_params()


#### 3. Pre-processing actions  ===========================================

# Set random seed
set.seed(l_params_all$seed)

# Load parameter mapping
load(prior_file) 

# Load targets
l_true_targets <- recursive_read_csv(targets_files)

# Reshape true targets and get mapping
l_true_reshaped <- reshape_calib_targets(l_true_targets, output_se = TRUE, output_map = TRUE)

# Load BayCANN calibrated parameters
calibrated_params_baycann <- read_csv(path_baycann) %>%
  dplyr::select(-lp) %>% # Remove last non-parameter column
  as.matrix()

# Get vector of ages for prevalence and incidence
v_ages_prevalence <- get_age_range(l_true_targets$prevalence)
v_ages_incidence <- get_age_range(l_true_targets$incidence)
v_ages <- list(prevalence = v_ages_prevalence,
               incidence = v_ages_incidence)

# Set number of cores to use
if(run_parallel) {
  registerDoParallel(cores=detectCores(logical = TRUE) - 2)  
}


################################################################################
###  Generate BayCANN outputs                                          ###
################################################################################
# Run model for each input parameter sample and get corresponding targets
if(run_parallel) {
  # Parallel processing
  stime <- system.time({
    out_calib_targets <- foreach(i=1:nrow(calibrated_params_baycann), .combine=rbind, 
                                 .inorder=FALSE, 
                                 .packages=c("data.table","tidyverse")) %dopar% {
                                   
                                   # Get row of parameters and calculate targets
                                   v_params_update <- calibrated_params_baycann[i,]
                                   v_calib_targets <- params_to_calib_targets(l_params_all, 
                                                                              v_params_update, 
                                                                              prior_map,
                                                                              v_ages,
                                                                              verbose = verbose)
                                   # Call item to save
                                   t(v_calib_targets)
                                 }
  })
  
  print(stime)
  
  closeAllConnections()
} else {
  
  out_calib_targets <- data.frame()
  start_time <- Sys.time()
  for (i in 1:nrow(calibrated_params_baycann)) {
    
    v_params_update <- calibrated_params_baycann[i,]
    v_calib_targets <- params_to_calib_targets(l_params_all, v_params_update, prior_map,
                                               v_ages, verbose = verbose)
    
    out_calib_targets <- rbind(out_calib_targets, t(v_calib_targets))
  }
  end_time <- Sys.time()
  print(paste('Simulation time:', end_time - start_time))
  
}

# Save the unscaled posterior samples
write.csv(out_calib_targets,
          file = outpath,
          row.names = FALSE)
