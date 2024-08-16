###########################  Generate BayCANN Sample  ##########################
#
#  Objective: Program to simulate parameter inputs and model outputs for 
#  BayCANN model calibration
########################### <<<<<>>>>> #########################################


#### 1.Libraries and functions  ==================================================
# Clean environment
rm(list = ls())

library(readxl)
library(data.table)
library(tidyverse)
library(lhs)
library(doParallel)
library(foreach)
library(assertthat)

###### 1.1 Load  functions =================================================
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)

#### 2. General parameters ========================================================

###### 2.1 file paths 
target_files <- list(prevalence = "data/prevalence_asymptomatic_cancer.csv",
                     incidence = "data/incidence_symptomatic_cancer.csv",
                     stage_distr = "data/stage_distr.csv")
outpath <- 'output/calibration/BayCANN'

###### 2.2 model parameters
n_cohort_calib <- 500000
seed_calib <- 42
n_samp <- 2000

#### 3. Pre-processing actions  ===========================================

# Load model parameters
l_params_init <- load_default_params()
verbose <- FALSE

# Load calibration parameters
l_params_calib <- load_calib_params(l_params_init,
                                    target_files = target_files,
                                    n_cohort_calib = n_cohort_calib,
                                    seed_calib = seed_calib,
                                    outpath = outpath)

# Set number of cores to use
if(is.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE"))) {
  # use the environment variable SLURM_NTASKS_PER_NODE to set
  # the number of cores to use
  registerDoParallel(cores=(Sys.getenv("SLURM_NTASKS_PER_NODE")))
} else {
  registerDoParallel(cores=detectCores(logical = TRUE) - l_params_calib$n_cores_reserved_local)  
}

#### 4. Generate random set of inputs  ===========================================

m_param_samp <- with(l_params_calib, {
  # Get number of params to calibrate
  n_param <- nrow(prior_map)
  
  # Sample unit Latin Hypercube
  m_lhs_unit <- randomLHS(n_samp, n_param)
  
  # Rescale to min/max of each parameter
  m_param_samp <- matrix(nrow = n_samp, ncol = n_param)
  for (i in 1:n_param) {
    m_param_samp[, i] <- qunif(m_lhs_unit[, i],
                               min = prior_map$min[i],
                               max = prior_map$max[i])
  }
  colnames(m_param_samp) <- prior_map$var_id
  
  return(m_param_samp)
}
)

#### 5. Generate corresponding outputs  ===========================================

# Run model for each input parameter sample and get corresponding targets
m_outputs <- with(l_params_calib, {
  m_outputs <- param_sample_to_outputs(
    m_param_samp, 
    fn = 'params_to_calib_targets', 
    param_arg_name = 'v_params_update',
    fn_other_args = list(
      l_params_all = l_params_all,
      param_map = prior_map,
      v_ages = v_ages))
  return(m_outputs)
}
)

closeAllConnections()

# Check for any NaN
validate_that(
  sum(sapply(m_param_samp, function(x) any(is.nan(x)))) == 0, 
  msg = 'Parameters include NaN')
validate_that(
  sum(sapply(m_outputs, function(x) any(is.nan(x)))) == 0, 
  msg = 'Outputs include NaN')

#### 6. Save data files  ===========================================

saveRDS(list(m_param_samp = m_param_samp, 
             m_calib_outputs = m_outputs), 
        file = file.path(outpath, 'BayCANN_sample.rds'))
