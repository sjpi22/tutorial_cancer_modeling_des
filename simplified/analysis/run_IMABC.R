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
run_parallel <- TRUE

###### 2.1 file paths 
targets_files <- list(prevalence = "data/prevalence_asymptomatic_cancer.csv",
  incidence = "data/incidence_symptomatic_cancer.csv",
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
if(run_parallel) {
  registerDoParallel(cores=detectCores(logical = TRUE) - 2)  
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
load('data/priors.RData')

# Load targets
l_true_targets <- recursive_read_csv(targets_files)

# Reshape true targets and get mapping
l_true_reshaped <- reshape_calib_targets(l_true_targets, output_se = TRUE, output_map = TRUE)

# Get vector of ages for prevalence and incidence
v_ages_prevalence <- get_age_range(l_true_targets$prevalence)
v_ages_incidence <- get_age_range(l_true_targets$incidence)
v_ages <- list(prevalence = v_ages_prevalence,
               incidence = v_ages_incidence)


#### 4. IMABC inputs  ===========================================

# Define Model Parameters/Priors
param_df <- data.frame(
  name_var = prior_map$var_id,
  dist_var = prior_map$prior_distr,
  min = prior_map$prior_min,
  max = prior_map$prior_max
)

priors <- as.priors(
  param_df,
  parameter_name = "name_var", dist_base_name = "dist_var"
)

# Multiplier for SD to get bounds
alpha_current = rep(0.0001, length(l_true_reshaped$v_targets))
alpha_current[l_true_reshaped$target_map$target_groups == 'Cancer incidence'] <- 1e-15
alpha_current[l_true_reshaped$target_map$target_groups == 'Stage at diagnosis'] <- 1e-15

alpha_stop = rep(0.05, length(l_true_reshaped$v_targets))
alpha_stop[l_true_reshaped$target_map$target_groups == 'Cancer incidence'] <- 1e-7
alpha_stop[l_true_reshaped$target_map$target_groups == 'Stage at diagnosis'] <- 1e-7


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
  v_targets <- params_to_calib_targets(l_params_all, v_params_update, prior_map, 
                                       v_ages)
  return(v_targets)
}

target_fun <- define_target_function(
  targets, priors, FUN = fn, use_seed = FALSE
)

# Calibrate model - see here for documentation: https://github.com/c-rutter/imabc/tree/a58a3b7c8db18948ff87fb6be55c6175399f41a2
start_time <- Sys.time()
calibration_results <- imabc(
  priors = priors,
  targets = targets,
  target_fun = target_fun,
  seed = l_params_all$seed,
  # N_start = 1000 * length(priors),
  N_start = 1000,
  N_centers = 1,
  Center_n = 100,
  # output_directory = outpath,
  N_post = 1000,
  verbose = TRUE,
  max_iter = 1000
  # validate_run = TRUE
)
end_time <- Sys.time()
print(end_time - start_time) # 29.41799 mins

# If output_directory is included, Getting error: Error in length(dots) == 2 && class(dots[[1]]) != "list" : 
# 'length = 2' in coercion to 'logical(1)'

if(debug) {
  print('Saving debug output')
  save(calibration_results, file = paste0(outpath, '/IMABC_calibration_debug.RData'))
} else {
  print('Saving output')
  save(calibration_results, file = paste0(outpath, '/IMABC_calibration.RData'))
}
