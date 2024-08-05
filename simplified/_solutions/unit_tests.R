# Unit tests

################################################################################
# Setup
################################################################################
# Clear workspace
rm(list = ls())

# Options
options(scipen=999)

# Load packages
library(readxl)
library(data.table)
library(tidyverse)
library(survival)
library(survminer)
library(assertthat)

# Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)

################################################################################
# Parameters
################################################################################

#### Modifiable parameters ####
n_cohort <- 100000 # Number to simulate in cohort
n_screen_sample <- 10000

# Randomization 
seed <- 1 # Random seed for generating data

# Outcome reporting
v_ages_prevalence <- seq(30, 80, 10) # Age for lesion prevalence
v_ages <- list(prevalence = v_ages_prevalence,
               incidence = c(v_ages_prevalence, 100)) # Age for cancer incidence 
v_time_surv <- seq(0, 10) # Times from event to calculate relative survival

################################################################################
# Test character cancer stages
################################################################################

# Load default data
l_params_all <- load_default_params(v_cancer = c('a', 'b'),
                                    file.surv = NULL)

# Add true survival distribution of exponential from diagnosis
l_params_all$time_Ca_Dc$distr <- "exp"
l_params_all$time_Ca_Dc$params <- list(rate = 0.1)
l_params_all$time_Cb_Dc$distr <- "exp"
l_params_all$time_Cb_Dc$params <- list(rate = 0.3)

# Map variables to parameters for tuning - make dataframe of all parameters with "src = unknown"
param_map <- make_param_map(l_params_all)

# Set "true" parameters
v_param_update <- c(0.42, 0.25, 0.5, 3, 300)

# Update params
l_params_all <- update_param_from_map(l_params_all, v_param_update, param_map)
l_params_all <- update_param_list(l_params_all,
                                  list(seed = seed,
                                       n_cohort = n_cohort,
                                       v_strats = l_params_all$v_strats[1]))

# Update parameter map
param_map$param_val <- v_param_update

#### Initialize population and disease natural history ####
results <- run_model(l_params_all)
results_noscreening <- results[['None']]

################################################################################
# Test more than 2 stages
################################################################################

# Update cancer stages
l_params_all$v_cancer = c('a', 'b', 'c', 'd')

# Add progression within preclinical stage variables
l_params_all$time_Pb_Pc$distr <- "exp"
l_params_all$time_Pb_Pc$params <- list(rate = 0.2)
l_params_all$time_Pc_Pd$distr <- "exp"
l_params_all$time_Pc_Pd$params <- list(rate = 0.3)

# Add detection variables
l_params_all$time_Pc_C$distr <- "exp"
l_params_all$time_Pc_C$params <- list(rate = 0.4)
l_params_all$time_Pd_C$distr <- "exp"
l_params_all$time_Pd_C$params <- list(rate = 0.6)

# Add true survival distribution of exponential from diagnosis
l_params_all$time_Cc_Dc$distr <- "exp"
l_params_all$time_Cc_Dc$params <- list(rate = 0.4)
l_params_all$time_Cd_Dc$distr <- "exp"
l_params_all$time_Cd_Dc$params <- list(rate = 0.5)

#### Initialize population and disease natural history ####
results <- run_model(l_params_all)
results_noscreening <- results[['None']]

# Check sums
are_equal(results_noscreening[stage_dx == 'a'], results_noscreening[time_Pa_Pb > time_Pa_C])
are_equal(results_noscreening[stage_dx == 'b'], results_noscreening[time_Pb_Pc > time_Pb_C])
are_equal(results_noscreening[stage_dx == 'c'], results_noscreening[time_Pc_Pd > time_Pc_C])

# are_equal(results_noscreening[stage_dx == 'a', time_H_Dc], results_noscreening[(time_H_P])
