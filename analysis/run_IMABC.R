# Run IMABC: Vignette at https://github.com/c-rutter/imabc

# Run once
install.packages("imabc")

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

#* Clean environment
rm(list = ls())

# Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)

#### 2. General parameters ========================================================

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

# Load parameter mapping
load(sample_file) # Only need parameter mapping param_map here

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



# Define Target Values


df <- data.frame(
  target_groups = make.names(l_true_reshaped$target_map$target_groups),
  target_names = names(l_true_reshaped$v_targets),
  targets = unname(l_true_reshaped$v_targets),
  current_lower_bounds = unname(l_true_reshaped$v_targets) - 10*qnorm(l_params_all$conf_level + (1-l_params_all$conf_level)/2) * unname(l_true_reshaped$v_se),
  current_upper_bounds = unname(l_true_reshaped$v_targets) + 10*qnorm(l_params_all$conf_level + (1-l_params_all$conf_level)/2) * unname(l_true_reshaped$v_se),
  stopping_lower_bounds = unname(l_true_reshaped$v_targets) - 10*unname(l_true_reshaped$v_se),
  stopping_upper_bounds = unname(l_true_reshaped$v_targets) + 10*unname(l_true_reshaped$v_se)
)
targets <- as.targets(df)

# Define Target Function
fn <- function(v_params_update) {
  v_targets <- params_to_calib_targets(l_params_all, v_params_update, param_map, 
                                       v_ages_prevalence, v_ages_incidence)
  return(v_targets)
}

target_fun <- define_target_function(
  targets, priors, FUN = fn, use_seed = FALSE
)

# Calibrate model
calibration_results <- imabc(
  priors = priors,
  targets = targets,
  target_fun = target_fun,
  seed = l_params_all$seed,
  # N_start = 2000,
  # N_centers = 2,
  # Center_n = 500,
  N_start = 500,
  N_centers = 1,
  Center_n = 100,
  N_cov_points = 50,
  N_post = 100,
  verbose = TRUE,
  max_iter = 10
)

# Getting error: Error in imabc(priors = priors, targets = targets, target_fun = target_fun,  : 
# No valid parameters to work from.
