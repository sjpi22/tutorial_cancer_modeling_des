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
data_inpath <- data_outpath <- 'data/'

###### 2.2 modifiable parameters
# Control variables for running on cluster and/or parallelized
run_parallel <- TRUE

# For debugging and viewing outputs
debug_small <- FALSE
debug_large <- FALSE
print_increment = 0.05

# Number of samples
if (debug_small) {
  n_samp <- 10
  n_cohort = 1000
} else if (debug_large) {
  n_samp <- 100
  n_cohort <- 100000
} else {
  n_samp <- 2000
  n_cohort <- 100000
}

# Set number of cores to use
if(run_parallel) {
  registerDoParallel(cores=detectCores(logical = TRUE) - 2)  
}


#### 3. Pre-processing actions  ===========================================

# Create directory if it does not exist
dir.create(file.path(data_outpath), showWarnings = FALSE)

# Load default data
l_params_all <- load_default_params()

# Make cohort small for testing
l_params_all <- update_param_list(l_params_all,
                                  list(n_cohort = n_cohort,
                                       v_strats = l_params_all$v_strats[1]))

# Set seed
set.seed(l_params_all$seed)

# Set prior distributions (set uniform for everything)
load(paste0(data_inpath, 'priors.RData'))

# Get number of params to calibrate
n_param <- nrow(prior_map)

# Load calibration targets for vectors of ages
true_prevalence <- read.csv(file = paste0(data_inpath, 'prevalence_asymptomatic_cancer.csv'))
true_incidence_cancer <- read.csv(file = paste0(data_inpath, 'incidence_symptomatic_cancer.csv'))

# Get vector of ages for prevalence and incidence
v_ages_prevalence <- get_age_range(true_prevalence)
v_ages_incidence <- get_age_range(true_incidence_cancer)
v_ages <- list(prevalence = v_ages_prevalence,
               incidence = v_ages_incidence)

#### 4. Generate random set of inputs  ===========================================

# Sample unit Latin Hypercube
m_lhs_unit <- randomLHS(n_samp, n_param)

# Rescale to min/max of each parameter
m_param_samp <- matrix(nrow = n_samp, ncol = n_param)
for (i in 1:n_param) {
  m_param_samp[, i] <- qunif(m_lhs_unit[, i],
                             min = prior_map$prior_min[i],
                             max = prior_map$prior_max[i])
}
colnames(m_param_samp) <- prior_map$var_id


#### 5. Generate corresponding outputs  ===========================================

# Run model for each input parameter sample and get corresponding targets
out_calib_targets <- data.frame()
verbose <- FALSE

if(run_parallel) {
  # Parallel processing
  stime <- system.time({
    out_calib_targets <- foreach(i=1:n_samp, .combine=rbind, 
                                 .inorder=FALSE, 
                                 .packages=c("data.table","tidyverse")) %dopar% {
                                   
                                   # Get row of parameters and calculate targets
                                   v_params_update <- m_param_samp[i,]
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
  start_time <- Sys.time()
  for (i in 1:n_samp) {
    # Print index for progress
    if (debug_small) {
      print('================')
      print(paste('Simulation', i))
      verbose = TRUE
    } else if (debug_large) {
      if (i <= 3) {
        print('================')
        print(paste('Simulation', i))
        verbose = TRUE
        if (i == 3) {
          end_time <- Sys.time()
          print(paste('Simulation time:', end_time - start_time))
        }
      } else verbose = FALSE
    } else {
      # For progress, print every 5% of the way
      if (round(i %% (n_samp * print_increment)) == 1) {
        print(paste0(round(i/n_samp * 100, 1), '% of simulations generated'))
      }
    }
    
    # Get row of parameters and calculate targets
    v_params_update <- m_param_samp[i,]
    v_calib_targets <- params_to_calib_targets(l_params_all, v_params_update, prior_map,
                                               v_ages,
                                               verbose = verbose)
    
    # Append target vector to dataframe
    out_calib_targets <- rbind(out_calib_targets, t(v_calib_targets))
  }
  end_time <- Sys.time()
  print(paste('Simulation time:', end_time - start_time))
}

# Check for any NaN
assertthat::validate_that(sum(sapply(m_param_samp, function(x) any(is.nan(x)))) == 0, 
                          msg = 'Parameters include NaN')
assertthat::validate_that(sum(sapply(out_calib_targets, function(x) any(is.nan(x)))) == 0, 
                          msg = 'Outputs include NaN')

#### 6. Save data files  ===========================================

if(!debug_small & !debug_large) {
  print('Saving output')
  save(m_param_samp, out_calib_targets, file = paste0(data_outpath, 'calibration_sample.RData'))
}
