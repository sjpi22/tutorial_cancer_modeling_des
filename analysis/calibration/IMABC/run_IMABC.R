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

###### 2.1 file paths 
target_files <- list(prevalence = "data/prevalence_asymptomatic_cancer.csv",
  incidence = "data/incidence_symptomatic_cancer.csv",
  stage_distr = "data/stage_distr.csv")

###### 2.12 model parameters 
n_cohort_calib <- 500000
seed_calib <- 42
outpath <- 'output/calibration/IMABC'


#### 3. Pre-processing actions  ===========================================

# Load model parameters
l_params_init <- load_default_params()

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

# Reshape true targets and get mapping
l_true_reshaped <- reshape_calib_targets(l_params_calib$l_true_targets, output_se = TRUE, output_map = TRUE)

#### 4. IMABC inputs  ===========================================

# Define Model Parameters/Priors
param_df <- with(l_params_calib, {
  data.frame(
    name_var = prior_map$var_id,
    dist_var = prior_map$distr,
    min = prior_map$min,
    max = prior_map$max
  )
})

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
  v_targets <- with(l_params_calib, {
    params_to_calib_targets(l_params_all, v_params_update, prior_map, 
                            v_ages)
  })
  return(v_targets)
}

target_fun <- define_target_function(
  targets, priors, FUN = fn, use_seed = FALSE
)

# Calibrate model - see here for documentation: https://github.com/c-rutter/imabc/tree/a58a3b7c8db18948ff87fb6be55c6175399f41a2
start_time <- Sys.time()
calibration_results <- with(l_params_calib, {
  imabc(
    priors = priors,
    targets = targets,
    target_fun = target_fun,
    seed = l_params_all$seed,
    N_start = 1000 * length(priors),
    N_centers = 1,
    Center_n = 100,
    # output_directory = outpath,
    N_post = 2000,
    verbose = TRUE,
    max_iter = 1000
    # validate_run = TRUE
  )
})
end_time <- Sys.time()
print(end_time - start_time) # 29.41799 mins

# If output_directory is included, Getting error: Error in length(dots) == 2 && class(dots[[1]]) != "list" : 
# 'length = 2' in coercion to 'logical(1)'

print('Saving output')
saveRDS(calibration_results, file = file.path(outpath, 'IMABC_outputs.rds'))
saveRDS(l_params_calib, file = file.path(outpath, 'IMABC_params.rds'))
