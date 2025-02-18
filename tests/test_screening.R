###########################  Unit test: Screening   ################
#
#  Objective: Run unit tests for screening evaluation
########################### <<<<<>>>>> ##############################################

rm(list = ls()) # Clean environment
options(scipen = 999) # View data without scientific notation

#### 1.Libraries and functions  ==================================================

###### 1.1 Load packages
library(tidyverse)
library(testthat)

###### 1.2 Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)


#### 2. General parameters ========================================================

set.seed(123)
eps_pct <- 0.01 # tolerated error
alpha <- 0.05 # p-value threshold

# Screening parameters
age_screen_start <- 45
age_screen_end <- 75
test_interval <- 10
p_spec <- 0.8

#### 3. Pre-processing  ===========================================

# Test expected number of screens and false positives in healthy state
v_time_H_L <- seq(age_screen_start - test_interval, age_screen_end + 2 * test_interval, test_interval / 2) + eps_pct

# Generate data
m_patient_base <- data.table(pt_id = length(v_time_H_L),
                             time_H_L = v_time_H_L,
                             time_H_P = v_time_H_L + test_interval)
m_patient_base[, time_H_C := time_H_P + runif(.N, min = 1, max = test_interval)]
m_patient_base[, `:=` (time_H_Dc = time_H_C + runif(.N, min = 1, max = test_interval),
                       time_H_Do = max(v_time_H_P))]
m_patient_base[, `:=` (time_H_D = pmin(time_H_Do, time_H_Dc),
                       time_screen_censor = pmin(time_H_Do, time_H_C),
                       n_screen_exp_L = pmax(0, floor((pmin(age_screen_end, time_H_L) - age_screen_start) / test_interval) + 1),
                       n_screen_exp_P = pmax(0, floor((pmin(age_screen_end, time_H_P) - age_screen_start) / test_interval) + 1))]

# Test number of screens
test_that("Match expected number of screens in healthy state", {
  # Create copy of data
  m_patient_screen <- copy(m_patient_base)
  
  # Equality with precancerous lesion state
  simulate_screening_H(m_patient_screen,
                       var_onset = "time_H_L",
                       age_screen_start = age_screen_start,
                       age_screen_end = age_screen_end,
                       int_screen = test_interval,
                       p_spec = p_spec)
  expect_equal(m_patient_screen[, n_screen_exp_L], m_patient_screen[, ct_screen])
  
  # Equality without precancerous lesion state
  simulate_screening_H(m_patient_screen,
                       var_onset = "time_H_P",
                       age_screen_start = age_screen_start,
                       age_screen_end = age_screen_end,
                       int_screen = test_interval,
                       p_spec = p_spec)
  expect_equal(m_patient_screen[, n_screen_exp_P], m_patient_screen[, ct_screen])
})


# Load default parameters
l_params_init <- load_model_params(
  lesion_state = TRUE
)

# Make parameter map and change parameters to be more realistic
df_params <- make_param_map(l_params_init)
df_params$param_val <- c(2, 75, 1/20, 1/10, 1/7, 1/35, 1/6, 1/16, 1/5, 1/10, 1/2)

# Update parameter list
l_params_model <- update_param_from_map(l_params_init, df_params$param_val, df_params)

# Simulate natural history model
res <- run_base_model(l_params_model)


#### 4. Tests  ===========================================





# Test expected number of screens, detections, and false positives in lesion state

# Test expected number of screens and detections in preclinical cancer state

# # Run screening 
# m_times_screen <- copy(m_times)
# run_screening_counterfactual(    age_screen_start = 45,
#                                  age_screen_end = 75,
#                                  int_screen = 10,
#                                  p_sens = list(P = 0.9),
#                                  p_spec = 0.95,
#                                  l_params_all = l_params_calib$l_params_all,
#                                  m_times)


### Test for lesion
# Load model parameters
l_params_init <- load_model_params(
  lesion_state = TRUE
)

l_params_init$d_time_H_L$params = list(shape = 2, scale = 75)

# Generate sample data
res <- run_base_model(l_params_init)

run_screening_counterfactual(    age_screen_start = 45,
                                 age_screen_end = 75,
                                 int_screen = 10,
                                 p_sens = list(L = 0.7,
                                               P = 0.9),
                                 p_spec = 0.95,
                                 l_params_all = l_params_init,
                                 m_times = res$patient_level,
                                 m_lesion = res$lesion_level)