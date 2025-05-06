###########################  Unit test: Epi functions   ################
#
#  Objective: Run unit tests for calculating epidemiological summary statistics,
#  including prevalence and lesion count
########################### <<<<<>>>>> ##############################################

rm(list = ls()) # Clean environment
options(scipen = 999) # View data without scientific notation

#### 1.Libraries and functions  ==================================================

###### 1.1 Load packages
library(testthat)
library(tidyverse)
library(readxl)
library(data.table)

###### 1.2 Load functions

# Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)

# Function to get expected value of a distribution
get_ev <- function(distr, lb, ub) {
  ev <- integrate(function(u) u*query_distr("d", u, distr$distr, distr$params),
                  lb, ub)$value
  return(ev)
}

# Function to get expected location of a distribution
get_el <- function(distr, lb, ub) {
  ev <- get_ev(distr, lb, ub)
  cdf_diff <- query_distr("p", ub, distr$distr, distr$params) - query_distr("p", lb, distr$distr, distr$params)
  return(ev / cdf_diff)
}


#### 2. General parameters ========================================================

###### 2.1 Configurations
# Run file to process configurations
source("configs/process_configs.R")
list2env(configs, envir = .GlobalEnv)

###### 2.2 File paths
path_truth <- "_ground_truth"
path_data <- paths$data
file_distr <- file.path(path_truth, "true_params.xlsx")

###### 2.3 Other parameters
seed <- 2025 # Random seed
conf_level <- 0.95
n_sim <- 40
time_1_2 <- 5 # Artificial time from disease onset to next state
v_outcomes_cs <- c("prevalence", "num_lesions") # Outcome types to calculate cross-sectionally


#### 3. Pre-processing  ===========================================

# Set seed
set.seed(seed)

# Load initial model parameters
l_params_init <- do.call(load_model_params, c(
  params_model
))

# Load calibration parameters
l_params_calib <- do.call(load_calib_params, c(
  l_params_model = list(l_params_init),
  params_calib
))

# Extract model parameters and remove decision model seed
l_params_model <- l_params_calib$l_params_model
l_params_model$seed <- NULL

# Get distribution for disease onset
var_onset <- paste0("time_H_", l_params_model$v_states[2])
d_var_onset <- l_params_model[[paste0("d_", var_onset)]]

# Get distribution for end of first disease state
var_end <- paste0("time_H_", l_params_model$v_states[3])
d_var_end <- l_params_model[[paste0("d_", var_end)]]


#### 4. Prevalence  ===========================================

full_summ_prevalence <- data.table()
full_summ_prevalence_bounded <- data.table()
full_summ_nlesions <- data.table()
for (i in 1:n_sim) {
  # Simulate cohort
  m_cohort <- run_base_model(l_params_model)
  
  # Separate patient and lesion data as necessary
  if (is.data.table(m_cohort)) {
    m_patients <- m_cohort
  } else {
    m_patients <- m_cohort$patient_level
    m_lesions <- m_cohort$lesion_level
  }
  
  # Add artificial upper bounds for tests
  m_patients[, `:=` (time_H_Inf = Inf,
                     time_H_2 = get(var_onset) + time_1_2)]
  
  # Create censor variables
  if (!is.null(params_calib$l_censor_vars)) {
    for (dt in names(params_calib$l_censor_vars)) {
      for (varname in names(params_calib$l_censor_vars[[dt]])) {
        get(dt)[, (varname) := do.call("pmin", c(.SD, na.rm = TRUE)),
                .SDcols = params_calib$l_censor_vars[[dt]][[varname]]]
      }
    }
  }
  
  ###### 4.1 Test that prevalence with unbounded disease duration is close to distribution of disease onset
  # Calculate prevalence (longitudinal)
  summ_prevalence <- calc_prevalence(
    m_patients, 
    var_onset, 
    "time_H_Inf", 
    "time_H_Inf", 
    v_ages = l_params_calib$l_params_outcome$prevalence$lit_params$v_ages)
  
  # Calculate prevalence (cross-sectional)
  summ_prevalence_cs <- calc_prevalence_cs(
    m_patients, 
    var_onset, 
    "time_H_Inf", 
    "time_H_Inf",
    v_ages = l_params_calib$l_params_outcome$prevalence$lit_params$v_ages)
  
  # Merge longitudinal and cross-sectional
  summ_prevalence <- merge(summ_prevalence,
                           summ_prevalence_cs, 
                           by = c("age_start", "age_end"),
                           suffixes = c("_long", "_cs"))
  
  # Calculate true prevalence
  summ_prevalence[, pct_true := mapply(function(lb, ub) {
    integrate(function(u) query_distr("p", u, d_var_onset$distr, d_var_onset$params), lb, ub)$value / (ub - lb)
  }, age_start, age_end)]
  
  # Append to full data table
  full_summ_prevalence <- rbind(full_summ_prevalence, summ_prevalence)
  
  ###### 4.2 Test prevalence with bounded disease duration
  # Calculate prevalence (longitudinal)
  summ_prevalence_bounded <- calc_prevalence(
    m_patients, 
    var_onset, 
    "time_H_2", 
    "time_H_Inf", 
    v_ages = l_params_calib$l_params_outcome$prevalence$lit_params$v_ages)
  
  # Calculate prevalence (cross-sectional)
  summ_prevalence_cs_bounded <- calc_prevalence_cs(
    m_patients, 
    var_onset, 
    "time_H_2",
    "time_H_Inf",
    v_ages = l_params_calib$l_params_outcome$prevalence$lit_params$v_ages)
  
  # Merge longitudinal and cross-sectional
  summ_prevalence_bounded <- merge(summ_prevalence_bounded,
                                   summ_prevalence_cs_bounded, 
                                   by = c("age_start", "age_end"),
                                   suffixes = c("_long", "_cs"))
  
  # Calculate theoretical prevalence within age range
  summ_prevalence_bounded[, pct_true := mapply(function(lb, ub) {
    integrate(function(u) (pmin(u + time_1_2, ub) - u)*query_distr("d", u, d_var_onset$distr, d_var_onset$params), lb, ub)$value
  }, age_start, age_end) +
    mapply(function(lb, ub) {
      integrate(function(u) (u + time_1_2 - ub)*query_distr("d", u, d_var_onset$distr, d_var_onset$params), lb, ub)$value
    }, age_start - time_1_2, age_start)]
  summ_prevalence_bounded[, pct_true := pct_true/(age_end - age_start)]
  
  # Append to full data table
  full_summ_prevalence_bounded <- rbind(full_summ_prevalence_bounded, summ_prevalence_bounded)
  
  ###### 4.3 Compare longitudinal and cross-sectional lesion count
  if (params_model$lesion_state == T) {
    # Calculate lesion count (longitudinal)
    outcome <- "n_lesions"
    summ_nlesions <- do.call(
      calc_nlesions,  
      c(lapply(l_params_calib$l_params_outcome[[outcome]][["get_params"]], get, envir = sys.frame(sys.parent(0))), 
        l_params_calib$l_params_outcome[[outcome]][["lit_params"]]))
    
    # Calculate lesion count (cross-sectional)
    summ_nlesions_cs <- do.call(
      calc_nlesions_cs,  
      c(lapply(l_params_calib$l_params_outcome[[outcome]][["get_params"]], get, envir = sys.frame(sys.parent(0))), 
        l_params_calib$l_params_outcome[[outcome]][["lit_params"]]))
    
    # Merge longitudinal and cross-sectional
    summ_nlesions <- merge(summ_nlesions,
                           summ_nlesions_cs, 
                           by = "n_lesions",
                           suffixes = c("_long", "_cs"))
    
    # Append to full data table
    full_summ_nlesions <- rbind(full_summ_nlesions, summ_nlesions)
    
  }
}

# (4.1) Check whether unbounded cross-sectional CI contains longitudinal and true prevalence
full_summ_prevalence[, `:=` (contained_long = value_long >= ci_lb & value_long <= ci_ub,
                             contained_true = pct_true >= ci_lb & pct_true <= ci_ub)]

test_that("Unbounded prevalence is close to distribution of disease onset", {
  expect_true(mean(full_summ_prevalence$contained_long) >= conf_level)
  expect_true(mean(full_summ_prevalence$contained_true) >= conf_level)
})

# (4.2) Check whether cross-sectional CI contains longitudinal and expected prevalence
full_summ_prevalence_bounded[, `:=` (contained_long = value_long >= ci_lb & value_long <= ci_ub,
                                     contained_true = pct_true >= ci_lb & pct_true <= ci_ub)]

test_that("Prevalence is close to incidence x duration", {
  expect_true(mean(full_summ_prevalence_bounded$contained_long) >= conf_level)
  expect_true(mean(full_summ_prevalence_bounded$contained_true)  >= conf_level)
})

# (4.3) Check whether cross-sectional lesion count is close to longitudinal
full_summ_nlesions[, `:=` (contained_long = value_long >= ci_lb & value_long <= ci_ub)]

test_that("Longitudinal and cross-sectional lesion count are close", {
  expect_true(mean(full_summ_nlesions$contained_long) >= conf_level)
})
