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

# Create directory if it does not exist
dir.create(file.path(data_outpath), showWarnings = FALSE)

# Load model parameters
l_params_init <- load_default_params()

# Load calibration parameters
l_params_calib <- load_calib_params(l_params_init,
                                    target_files = target_files,
                                    n_cohort_calib = n_cohort_calib,
                                    seed_calib = seed_calib,
                                    outpath = outpath)

# Set seed
set.seed(l_params_all$seed)

# Set prior distributions (set uniform for everything)
prior_map <- readRDS(file.path(data_inpath, 'priors.rds'))

# Get number of params to calibrate
n_param <- nrow(prior_map)

# Set number of cores to use
if(is.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE"))) {
  # use the environment variable SLURM_NTASKS_PER_NODE to set
  # the number of cores to use
  registerDoParallel(cores=(Sys.getenv("SLURM_NTASKS_PER_NODE")))
} else {
  registerDoParallel(cores=detectCores(logical = TRUE) - l_params_calib$n_cores_reserved_local)  
}

#### 4. Generate random set of inputs  ===========================================

# Sample unit Latin Hypercube
m_lhs_unit <- randomLHS(n_samp, n_param)

# Rescale to min/max of each parameter
m_BayCANN_param_samp <- matrix(nrow = n_samp, ncol = n_param)
for (i in 1:n_param) {
  m_BayCANN_param_samp[, i] <- qunif(m_lhs_unit[, i],
                                     min = prior_map$prior_min[i],
                                     max = prior_map$prior_max[i])
}
colnames(m_BayCANN_param_samp) <- prior_map$var_id


#### 5. Generate corresponding outputs  ===========================================

# Run model for each input parameter sample and get corresponding targets
m_BayCANN_calib_outputs <- data.frame()

# Parallel processing
stime <- system.time({
  m_BayCANN_calib_outputs <- foreach(i=1:n_samp, .combine=rbind, 
                                     .inorder=FALSE, 
                                     .packages=c("data.table","tidyverse")) %dopar% {
                                       
                                       # Get row of parameters and calculate targets
                                       v_params_update <- m_BayCANN_param_samp[i,]
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


# Check for any NaN
assertthat::validate_that(sum(sapply(m_BayCANN_param_samp, function(x) any(is.nan(x)))) == 0, 
                          msg = 'Parameters include NaN')
assertthat::validate_that(sum(sapply(m_BayCANN_calib_outputs, function(x) any(is.nan(x)))) == 0, 
                          msg = 'Outputs include NaN')

#### 6. Save data files  ===========================================

save(m_BayCANN_param_samp, m_BayCANN_calib_outputs, file = file.path(data_outpath, 'BayCANN_sample.RData'))
