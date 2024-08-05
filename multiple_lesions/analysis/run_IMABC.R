# Run IMABC: Vignette at https://github.com/c-rutter/imabc

# Run once
# install.packages("imabc")

# Load packages
library(imabc)
library(MASS)
library(data.table)
library(foreach)
library(parallel)
library(truncnorm)
library(lhs)
library(methods)
library(stats)
library(utils)
library(readxl)
library(tidyverse)
library(doParallel)

#* Clean environment
rm(list = ls())

# Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)

#### 2. General parameters ========================================================
debug <- FALSE
run_slurm <- TRUE
run_parallel <- TRUE

###### 2.1 file paths 
sample_file <- "data/calibration_sample.RData"

targets_files <- list(prevalence = list(
  a = "data/prevalence_lesion_a.csv",
  b = "data/prevalence_lesion_b.csv"),
  incidence = "data/incidence_cancer.csv",
  stage_distr = "data/stage_distr.csv")

###### 2.1 model parameters 

# Load default data
l_params_all <- load_default_params()

# Update defaults
l_params_all <- update_param_list(l_params_all,
                                  list(n_cohort = 100000,
                                       v_strats = l_params_all$v_strats[1]))

#### 3. Pre-processing actions  ===========================================

# Set random seed
set.seed(l_params_all$seed)

# Set number of cores to use
if(run_slurm) {
  # use the environment variable SLURM_NTASKS_PER_NODE to set
  # the number of cores to use
  registerDoParallel(cores=(Sys.getenv("SLURM_NTASKS_PER_NODE")))
  outpath <- 'sherlock/output/imabc'
  
  # Create directory if it does not exist
  dir.create(file.path(outpath))
} else if(run_parallel) {
  registerDoParallel(cores=min(detectCores(logical = TRUE), 6) - 2)  
  outpath <- 'output/imabc'
  
  # Create directory if it does not exist
  dir.create(file.path(outpath))
}

# File path for debug
if (debug) {
  outpath <- paste0(outpath, '/debug')
  
  # Create directory if it does not exist
  dir.create(file.path(outpath))
}

# Load parameter mapping
load(sample_file) # Only need parameter mapping param_map here

# If debug, use narrow range around true values
if (debug) {
  # Resave param map
  model_param_map <- param_map
  
  # Load true parameters
  true_param_path <- 'disease_simulator/true_param_map_consistent.RData'
  load(true_param_path)
  model_param_map$prior_min <- param_map$param_val * 0.9
  model_param_map$prior_max <- param_map$param_val * 1.2
  
  param_map <- model_param_map
}

# Load targets
l_true_targets <- recursive_read_csv(targets_files)

# Reshape true targets and get mapping
l_true_reshaped <- reshape_calib_targets(l_true_targets, output_se = TRUE, output_map = TRUE)

# Get vector of ages for prevalence
v_ages_prevalence <- list()
v_ages_prevalence[['a']] <- get_age_range(l_true_targets$prevalence[['a']])
v_ages_prevalence[['b']] <- get_age_range(l_true_targets$prevalence[['b']])

# Get vector of ages for incidence
v_ages_incidence <- get_age_range(l_true_targets$incidence)


#### 4. IMABC inputs  ===========================================


# Define Model Parameters/Priors
param_df <- data.frame(
  name_var = param_map$var_id,
  dist_var = param_map$prior_distr,
  min = param_map$prior_min,
  max = param_map$prior_max
)

priors <- as.priors(
  param_df,
  parameter_name = "name_var", dist_base_name = "dist_var"
)

# Multiplier for SD to get bounds
if (debug) {
  alpha_current = rep(0.0001, length(l_true_reshaped$v_targets))
  alpha_current[l_true_reshaped$target_map$target_groups == 'Cancer incidence'] <- 1e-30
  alpha_current[l_true_reshaped$target_map$target_groups == 'Stage at diagnosis'] <- 1e-30
  
  alpha_stop = rep(0.001, length(l_true_reshaped$v_targets))
  alpha_stop[l_true_reshaped$target_map$target_groups == 'Cancer incidence'] <- 1e-10
  alpha_stop[l_true_reshaped$target_map$target_groups == 'Stage at diagnosis'] <- 1e-15
} else {
  alpha_current = rep(0.0001, length(l_true_reshaped$v_targets))
  alpha_current[l_true_reshaped$target_map$target_groups == 'Cancer incidence'] <- 1e-15
  alpha_current[l_true_reshaped$target_map$target_groups == 'Stage at diagnosis'] <- 1e-30
  
  alpha_stop = rep(0.05, length(l_true_reshaped$v_targets))
  alpha_stop[l_true_reshaped$target_map$target_groups == 'Cancer incidence'] <- 1e-7
  alpha_stop[l_true_reshaped$target_map$target_groups == 'Stage at diagnosis'] <- 1e-15
}

# Define Target Values
target_df <- data.frame(
  target_groups = make.names(l_true_reshaped$target_map$target_groups),
  target_names = names(l_true_reshaped$v_targets),
  targets = unname(l_true_reshaped$v_targets),
  current_lower_bounds = pmax(0, 
                              unname(l_true_reshaped$v_targets) - qnorm((1-alpha_current) + alpha_current/2) * unname(l_true_reshaped$v_se)),
  current_upper_bounds = unname(l_true_reshaped$v_targets) + qnorm((1-alpha_current) + alpha_current/2) * unname(l_true_reshaped$v_se),
  stopping_lower_bounds = pmax(0, 
                               unname(l_true_reshaped$v_targets) - qnorm((1-alpha_stop) + alpha_stop/2) * unname(l_true_reshaped$v_se)),
  stopping_upper_bounds = unname(l_true_reshaped$v_targets) + qnorm((1-alpha_stop) + alpha_stop/2) * unname(l_true_reshaped$v_se)
)

targets <- as.targets(target_df)

# Define Target Function
fn <- function(v_params_update) {
  v_targets <- params_to_calib_targets(l_params_all, v_params_update, param_map, 
                                       v_ages_prevalence, v_ages_incidence)
  return(v_targets)
}

target_fun <- define_target_function(
  targets, priors, FUN = fn, use_seed = FALSE
)

# Calibrate model - see here for documentation: https://github.com/c-rutter/imabc/tree/a58a3b7c8db18948ff87fb6be55c6175399f41a2
calibration_results <- imabc(
  priors = priors,
  targets = targets,
  target_fun = target_fun,
  seed = l_params_all$seed,
  N_start = 1000 * length(priors),
  N_centers = 10,
  Center_n = 1000,
  output_directory = outpath,
  N_post = 5000,
  verbose = TRUE,
  max_iter = 1000,
  validate_run = TRUE
)

if(debug) {
  print('Saving debug output')
  save(calibration_results, file = paste0(outpath, 'IMABC_calibration_debug.RData'))
} else {
  print('Saving output')
  save(calibration_results, file = paste0(outpath, 'IMABC_calibration.RData'))
}

# Getting error: Error in imabc(priors = priors, targets = targets, target_fun = target_fun,  : 
# No valid parameters to work from.