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

###### 2.2 model parameters 
n_cohort_calib <- 500000
seed_calib <- 42

###### 2.3 IMABC parameters 
outpath <- 'output/calibration/IMABC'
alpha_current <- c(1e-4, 1e-15, 1e-15)
alpha_stop = c(0.05, 1e-7, 1e-7)
fn_use_seed <- FALSE
optional_args = list(
  N_centers = 1,
  Center_n = 100,
  N_post = 2000
)

#### 3. Pre-processing actions  ===========================================

# Load model parameters
l_params_init <- load_default_params()

# Load calibration parameters
l_params_calib <- load_calib_params(l_params_init,
                                    target_files = target_files,
                                    n_cohort_calib = n_cohort_calib,
                                    seed_calib = seed_calib,
                                    outpath = outpath)

# Add IMABC-specific parameters
l_params_calib <- add_IMABC_params(l_params_calib,
                                   alpha_current = alpha_current,
                                   alpha_stop = alpha_stop,
                                   use_seed = fn_use_seed,
                                   optional_args = optional_args)

# Set number of cores to use
if(is.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE"))) {
  # use the environment variable SLURM_NTASKS_PER_NODE to set
  # the number of cores to use
  registerDoParallel(cores=(Sys.getenv("SLURM_NTASKS_PER_NODE")))
} else {
  registerDoParallel(cores=detectCores(logical = TRUE) - l_params_calib$n_cores_reserved_local)  
}

#### 4. Run IMABC  ===========================================

# Calibrate model - see here for documentation: https://github.com/c-rutter/imabc/tree/a58a3b7c8db18948ff87fb6be55c6175399f41a2
start_time <- Sys.time()
calibration_results <- with(l_params_calib, {
  do.call(imabc, imabc_inputs)
})
end_time <- Sys.time()
print(end_time - start_time)

print('Saving output')
saveRDS(calibration_results, file = file.path(outpath, 'IMABC_outputs.rds'))
saveRDS(l_params_calib, file = file.path(outpath, 'IMABC_params.rds'))
