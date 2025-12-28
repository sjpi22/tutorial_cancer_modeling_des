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
library(readxl)
library(data.table)

###### 1.2 Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)


#### 2. General parameters ========================================================

###### 2.1 Configurations
# Load configs
file_configs <- file.path("configs", "configs_simulated.yaml")
configs <- load_configs(file_configs)
list2env(configs, envir = .GlobalEnv)

# Extract relevant parameters from configs
file_params_calib <- configs$paths$file_params_calib
params_screen <- configs$params_screen

# Get list of relevant output file paths and load to global environment
l_filepaths <- update_config_paths("files_imabc", configs$paths)
list2env(l_filepaths, envir = .GlobalEnv)

# Load coverage analysis parameters from configs file
list2env(configs$params_coverage, envir = .GlobalEnv)

# Other parameters
file_times <- "tests/test_screening_data.xlsx"
file_true_params <- file.path("_ground_truth", "params_model.rds")
nreps_var <- 500 # Replicates for testing variation
p_alpha <- 0.05 # p-value threshold


#### 3. Pre-processing  ===========================================

# Load times and set key
m_patient_base <- read_excel(file_times) %>% setDT()
m_lesion_base <- read_excel(file_times, sheet = 2) %>% setDT()
setkey(m_patient_base, pt_id)
setkey(m_lesion_base, pt_id)

# Load ground truth parameters
l_params_true <- readRDS(file_true_params)
l_params_true$fl_count_tests <- TRUE

# Load model and calibration parameters
l_params_calib <- readRDS(file_params_calib)

# Set flag to count diagnostic tests in base case scenario
l_params_calib$l_params_model$fl_count_tests <- TRUE

# Load IMABC posteriors
l_posteriors <- readRDS(file_posterior)
m_params <- l_posteriors$good_parm_draws %>%
  dplyr::select(l_params_calib$prior_map$var_id)

# Set decision outcome parameters
l_params_calib$l_params_outcome <- params_screen$l_outcome_base

# Set screening test and strategy parameters
l_params_screen <- list(
  test_chars = list(
    confirm = list(
      p_sens = list(
        L = 1,
        P = 1
      ),
      p_spec = 1)),
  strats = list(
    strat_gold = list(
      mod = "confirm",
      age_screen_start = 50,
      age_screen_stop = 70,
      int_screen = 10
    ),
    strat_screen = list(
      mod = "confirm",
      age_screen_start = 50,
      age_screen_stop = 70,
      int_screen = 1
    )
  )
)

# Set counterfactual comparison parameters
l_params_outcome_counter <- params_screen$l_outcome_counterfactual

# Set seed
set.seed(l_params_calib$l_params_model$seed)


#### 4. Generate outputs  ===========================================

# Set censor date
m_patient_base[, `:=` (time_screen_censor = pmin(time_H_D, time_H_C))]

# Test screening in healthy state
test_that("Match expected number of screens in healthy state", {
  # Equality with (L) or without (P) precancerous lesion state
  for (dx in c("L", "P")) {
    # 1 - gold standard test, 2 - screening test
    for (strat in 1:2) {
      for (spec in c(0, 1)) {
        # Modify specificity
        l_params_screen$test_chars$confirm$p_spec <- spec
        
        # Create copy of data
        m_patient_screen <- copy(m_patient_base)
        
        # Simulate screening in healthy state
        simulate_screening_H(m_patient_screen,
                             var_onset = paste0("time_H_", dx),
                             l_params_strategy = l_params_screen$strats[[strat]],
                             l_test_chars = l_params_screen$test_chars)
        
        # Screening test count
        expect_equal(m_patient_screen[, get(paste0("ct_", l_params_screen$strats[[strat]]$int_screen, "_H", dx))], m_patient_screen[, ct_tests_screen])
        
        # False positive count
        if (spec == 1) {
          expect_equal(m_patient_screen$ct_tests_positive, rep(0, nrow(m_patient_screen))) # Expect 0 positive tests if perfect specificity
        } else {
          expect_equal(m_patient_screen$ct_tests_positive, m_patient_screen[, ct_tests_screen]) # Expect every screening is positive if zero specificity
        }
      }
    }
  }
})

# Test total screens (no-lesion cancer)
test_that("Match expected number of total screens for non-lesion cancer", {
  dx <- "P" # Only for non-lesion cancer
  # 1 - gold standard test, 2 - screening test
  for (strat in 1:2) {
    for (spec in c(0, 1)) {
      # Modify specificity
      l_params_screen$test_chars$confirm$p_spec <- spec
      
      # Create copy of data
      m_patient_screen <- copy(m_patient_base)
      
      # Run screening without surveillance
      run_screening_counterfactual(m_patient_screen,
                                   l_params_model = modifyList(l_params_calib$l_params_model, list(v_states = c("H", "P", "C", "D"))),
                                   l_params_strategy = l_params_screen$strats[[strat]],
                                   l_params_tests = l_params_screen$test_chars)
      
      # Screening and base tests
      expect_equal(m_patient_screen[, get(paste("ct", l_params_screen$strats[[strat]]$int_screen, "s", dx, sep = "_"))], m_patient_screen$ct_tests_screen)
      expect_equal(m_patient_screen[, get(paste("ct", l_params_screen$strats[[strat]]$int_screen, "b", dx, sep = "_"))], m_patient_screen$ct_tests_base)
      
      # False positive count
      if (spec == 1) {
        expect_equal(m_patient_screen[, get(paste("ct", l_params_screen$strats[[strat]]$int_screen, "pos", dx, sep = "_"))], m_patient_screen$ct_tests_positive) # Expect only positive test if developed lesion before death and perfect specificity
      } else {
        expect_equal(m_patient_screen$ct_tests_screen, m_patient_screen$ct_tests_positive) # Expect every screening is positive if zero specificity
      }
    }
  }
})

# Test screening with lesion state
test_that("Match expected number of screens in lesion state", {
  dx <- "L"
  # 1 - gold standard test, 2 - screening test
  for (strat in 1:2) {
    for (spec in c(0, 1)) {
      # Modify specificity
      l_params_screen$test_chars$confirm$p_spec <- spec
      
      # Create copy of data
      m_patient_screen <- copy(m_patient_base)
      m_lesion_screen <- copy(m_lesion_base)
      
      # Combine with lesion data
      m_cohort <- list(patient_level = m_patient_screen,
                       lesion_level = m_lesion_screen)
      
      # Run screening with no surveillance
      run_screening_counterfactual(m_cohort,
                                   l_params_model = l_params_calib$l_params_model,
                                   l_params_strategy = l_params_screen$strats[[strat]],
                                   l_params_tests = l_params_screen$test_chars)
      
      # Screening and base test count
      expect_equal(m_patient_screen[, get(paste("ct", l_params_screen$strats[[strat]]$int_screen, "s", dx, sep = "_"))], m_patient_screen[, ct_tests_screen])
      expect_equal(m_patient_screen[, get(paste("ct", l_params_screen$strats[[strat]]$int_screen, "b", dx, sep = "_"))], m_patient_screen[, ct_tests_base])
      
      # Detection time
      expect_equal(m_cohort$lesion_level[[paste("test", l_params_screen$strats[[strat]]$int_screen, "time_detected", sep = "_")]], m_cohort$lesion_level$time_detected)
      
      # False positive count
      if (spec == 1) {
        expect_equal(m_patient_screen[, get(paste("ct", l_params_screen$strats[[strat]]$int_screen, "pos", dx, sep = "_"))], m_patient_screen$ct_tests_positive) # Expect only positive test if developed lesion before death and perfect specificity
      } else {
        expect_equal(m_patient_screen$ct_tests_screen, m_patient_screen$ct_tests_positive) # Expect every screening is positive if zero specificity
      }
    }
  }
})

# Test imperfect screening
test_that("Match expected number of screens with variation", {
  # Set strategy
  strat <- 1
  
  # Set imperfect sensitivity
  p_sens_L <- 0.5
  p_sens_P <- 0.7
  l_params_screen$test_chars$confirm$p_sens <- list(
    L = p_sens_L,
    P = p_sens_P
  )
  
  # Set patients to extract
  v_pts <- c(3,4)
  
  # Select only 3rd and 4th row of data
  m_patient_subset <- list(
    patient_level = m_patient_base[pt_id %in% v_pts, ],
    lesion_level = m_lesion_base[pt_id %in% v_pts, ]
  )
  
  # Create repeated version of data
  m_cohort <- list()
  for (i in names(m_patient_subset)) {
    # Stack copies of patient and lesion data
    m_cohort[[i]] <- do.call(rbind, replicate(nreps_var, m_patient_subset[[i]], simplify=FALSE))
    
    # Save original patient ID and update patient IDs to be unique
    setnames(m_cohort[[i]], "pt_id", "pt_id_orig")
  }
  
  # Set new patient ID
  m_cohort$patient_level[, pt_id := 1:.N]
  m_cohort$lesion_level[, pt_id := rep(1:nrow(m_cohort$patient_level), rep(m_patient_subset$lesion_level[, .N, by = pt_id]$N, nrow(m_cohort$patient_level)/2))]
  
  # Set keys
  setkey(m_cohort[["patient_level"]], pt_id)
  setkey(m_cohort[["lesion_level"]], pt_id, lesion_id)
  
  # Run screening with no surveillance
  run_screening_counterfactual(m_cohort,
                               l_params_model = l_params_calib$l_params_model,
                               l_params_strategy = l_params_screen$strats[[strat]],
                               l_params_tests = l_params_screen$test_chars
  )
  
  # Check summary statistics for patients
  ct_cols <- c("time_H_P", "time_H_C", "time_H_D", "ct_tests_screen", "ct_tests_positive", "ct_tests_base")
  m_patient_summ_P <- m_cohort$patient_level[, .(ct = .N), by = c("pt_id_orig", "time_H_P")]
  m_patient_summ_C <- m_cohort$patient_level[, .(ct = .N), by = c("pt_id_orig", "time_H_C", "ct_tests_base")]
  m_patient_summ_C[, max_time := max(time_H_C), by = pt_id_orig]
  m_patient_summ_D <- m_cohort$patient_level[, .(avg = mean(time_H_D)), by = c("pt_id_orig")]
  m_patient_summ_ct <- m_cohort$patient_level[, .(ct = .N), by = c("pt_id_orig", "ct_tests_screen")]
  
  # Checks for patient 3
  expect_equal(m_patient_summ_P[pt_id_orig == v_pts[1], ct], nreps_var) # No variation in time_H_P
  expect_equal(m_patient_summ_ct[pt_id_orig == v_pts[1], ct], nreps_var) # No variation in test count
  expect_gt(binom.test(m_patient_summ_C[pt_id_orig == v_pts[1] & time_H_C == 50, ct], nreps_var, p_sens_P)$p.value, p_alpha) # Some variation in time_H_C
  
  # Checks for patient 4
  expect_equal(sort(m_patient_summ_P[pt_id_orig == v_pts[2], time_H_P]), sort(unique(c(NA, m_patient_subset$lesion_level[pt_id %in% v_pts[2], time_H_Pj]))))
  
  # Checks for both - base test should only be applied for cancer diagnosis at highest age
  expect_equal(nrow(m_patient_summ_C[ct_tests_base %in% 1 & time_H_C != max_time]), 0)
  expect_equal(sum(m_patient_summ_D[, avg] <= m_patient_subset$patient_level$time_H_D), 0) # Some variation in time_H_D
})

# Test that base data is unchanged for multiple screening strategies
test_that("Data reversion after screening modifications is successful", {
  # Change strategies to be 3 of the same strategy
  l_params_screen$strats <- list(strat1 = l_params_screen$strats$strat_gold,
                                 strat2 = l_params_screen$strats$strat_gold,
                                 strat3 = l_params_screen$strats$strat_gold)
  
  # Set seed
  l_params_screen$seed <- 2025
  
  # Run model with true parameters
  # If params_to_outputs() internal unit test passed, there should be no errors
  l_calib_outputs <- with(l_params_calib, {
    params_to_outputs(
      l_params_model = l_params_true,
      param_map = prior_map,
      l_params_outcome = l_params_outcome,
      l_params_screen = l_params_screen,
      l_params_outcome_counter = l_params_outcome_counter,
      l_censor_vars = l_censor_vars,
      reshape_output = FALSE,
      unit_test = TRUE
    )
  })
  
  # Outputs for every test are the same
  for (i in 2:length(l_calib_outputs$outputs_screen)) {
    expect_equal(l_calib_outputs$outputs_screen[[1]], l_calib_outputs$outputs_screen[[i]])
  }
})
