###########################  Load Calibration Parameters  ##########################
#
#  Objective: Program to load general calibration parameters
########################### <<<<<>>>>> #########################################


#### 1.Libraries and functions  ==================================================
# Clean environment
rm(list = ls())

library(readxl)
library(tidyverse)
library(lhs)
library(assertthat)
library(data.table)

###### 1.1 Load functions =================================================
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)


#### 2. General parameters ========================================================

###### 2.1 file paths 
target_files <- list(
  prevalence = "data/prevalence_asymptomatic_cancer.csv",
  incidence = "data/incidence_symptomatic_cancer.csv",
  stage_distr = "data/stage_distr.csv")
prior_file <- "data/priors.rds"
outpath <- "data"
file_params <- file.path(outpath, 'calibration_params.rds')

###### 2.2 model parameters
n_cohort_calib <- 200000
seed_calib <- 42

# Define parameters for calculating outcomes
# Names of lists should be the same as those of target_files
l_outcome_params <- list(
  prevalence = list(
    outcome_name = "prevalence",
    m_patients = "m_patients",
    start_var = "time_H_P", 
    end_var = "time_H_C", 
    censor_var = "time_screen_censor"),
  incidence = list(
    outcome_name = "incidence",
    m_patients = "m_patients",
    time_var = "time_H_C", 
    censor_var = "time_H_D"),
  stage_distr = list(
    outcome_name = "distr",
    m_patients = "m_patients",
    grouping_var = "stage_dx", 
    event_var = "time_H_C", 
    censor_var = "time_H_D",
    groups_expected = 1:4)
)

l_censor_vars <- list(
  time_screen_censor = c("time_H_C", "time_H_D")
)


#### 3. Load and save parameters  ===========================================

# Load model parameters
l_params_init <- load_model_params()

# Load calibration parameters
l_params_calib <- load_calib_params(
  l_params_init,
  target_files = target_files,
  l_outcome_params = l_outcome_params,
  l_censor_vars = l_censor_vars,
  prior_file = prior_file,
  n_cohort_calib = n_cohort_calib,
  seed_calib = seed_calib,
  outpath = outpath
)

saveRDS(l_params_calib, file = file_params)
