###########################  IMABC outputs   ###################################
#
#  Objective: Script to regenerate calibration target outputs for IMABC 
# calibrated parameters with cohort size consistent with targets
########################### <<<<<>>>>> #########################################


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
target_files <- list(prevalence = "data/prevalence_asymptomatic_cancer.csv",
                     incidence = "data/incidence_symptomatic_cancer.csv",
                     stage_distr = "data/stage_distr.csv")
inpath <- 'output/calibration/IMABC/IMABC_outputs.rds'
outpath <- 'output/calibration/IMABC'
outfile <- 'IMABC_targets_resampled.rds'

###### 2.2 model parameters 
n_cohort <- 100000
seed_calib <- 42

#### 3. Pre-processing actions  ===========================================

# Load model parameters
l_params_init <- load_default_params()

# Load calibration parameters
l_params_calib <- load_calib_params(l_params_init,
                                    target_files = target_files,
                                    n_cohort_calib = n_cohort,
                                    seed_calib = seed_calib,
                                    outpath = outpath)

# Load calibration outputs and extract parameter samples
calibration_outputs <- readRDS(inpath)
m_param_samp <- calibration_outputs$good_parm_draws %>%
  dplyr::select(l_params_calib$prior_map$var_id)

# Set number of cores to use
if(is.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE"))) {
  # use the environment variable SLURM_NTASKS_PER_NODE to set
  # the number of cores to use
  registerDoParallel(cores=(Sys.getenv("SLURM_NTASKS_PER_NODE")))
} else {
  registerDoParallel(cores=detectCores(logical = TRUE) - l_params_calib$n_cores_reserved_local)  
}

################################################################################
###  Generate IMABC outputs                                          ###
################################################################################

# Run model for each input parameter sample and get corresponding targets
m_outputs <- with(l_params_calib, {
  m_outputs <- param_sample_to_outputs(
    m_param_samp, 
    fn = 'params_to_calib_targets', 
    param_arg_name = 'v_params_update',
    fn_other_args = list(
      l_params_all = l_params_all,
      param_map = prior_map,
      v_ages = v_ages),
    run_parallel = TRUE)
  return(m_outputs)
}
)

# Save the posterior outputs
saveRDS(m_outputs, file = file.path(outpath, outfile))
