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
# Run file to process configurations
source("configs/process_configs.R")

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
nreps_var <- 500 # Replicates for testing variation


#### 3. Pre-processing  ===========================================

# Load times and set key
m_patient_base <- read_excel(file_times) %>% setDT()
m_lesion_base <- read_excel(file_times, sheet = 2) %>% setDT()
setkey(m_patient_base, pt_id)
setkey(m_lesion_base, pt_id)

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
      p_spec = 1,
      type = "direct"),
    screen = list(
      p_sens = list(
        L = 1,
        P = 1
      ),
      p_spec = 1,
      type = "indirect")),
  strats = list(
    strat_gold = list(
      mod = "confirm",
      age_screen_start = 50,
      age_screen_stop = 71,
      int_screen = 10
    ),
    strat_screen = list(
      mod = "screen",
      age_screen_start = 50,
      age_screen_stop = 71,
      int_screen = 5,
      mod_conf = "confirm",
      int_conf = 10
    )
  ),
  surveil = params_screen$surveil
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
      # Create copy of data
      m_patient_screen <- copy(m_patient_base)
      
      # Simulate screening in healthy state
      simulate_screening_H(m_patient_screen,
                           var_onset = paste0("time_H_", dx),
                           l_params_strategy = l_params_screen$strats[[strat]],
                           l_test_chars = l_params_screen$test_chars)
      
      # Screening test count
      if (strat == 1) {
        expect_equal("ct_tests_screen" %in% names(m_patient_screen), FALSE)
        expect_equal(m_patient_screen[, get(paste0("ct_", strat, "_c_H", dx))], m_patient_screen[, ct_tests_confirm])
      } else {
        expect_equal(m_patient_screen[, get(paste0("ct_", strat, "_s_H", dx))], m_patient_screen[, ct_tests_screen])
        expect_equal(m_patient_screen[, ct_tests_confirm], rep(0, nrow(m_patient_base)))
      }
    }
  }
})

# Test total screens (no-lesion cancer)
test_that("Match expected number of total screens for non-lesion cancer", {
  dx <- "P" # Only for non-lesion cancer
  # 1 - gold standard test, 2 - screening test
  for (strat in 1:2) {
    # Create copy of data
    m_patient_screen <- copy(m_patient_base)
    
    # Run screening without surveillance
    run_screening_counterfactual(m_patient_screen,
                                 l_params_model = modifyList(l_params_calib$l_params_model, list(v_states = c("H", "P", "C", "D"))),
                                 l_params_strategy = l_params_screen$strats[[strat]],
                                 l_params_test = l_params_screen$test_chars,
                                 l_params_surveil = NULL)
    
    # Screening test count
    if (strat == 1) {
      expect_equal("ct_tests_screen" %in% names(m_patient_screen), FALSE)
    } else {
      expect_equal(m_patient_screen[, get(paste0("ct_", strat, "_s_P"))], m_patient_screen[, ct_tests_screen])
    }
    
    # Confirmatory and base tests
    expect_equal(m_patient_screen[, get(paste0("ct_", strat, "_c_P"))], m_patient_screen[, ct_tests_confirm])
    expect_equal(m_patient_screen[, get(paste0("ct_", strat, "_b_P"))], m_patient_screen[, ct_tests_base])
  }
})

# Test screening during lesion state
test_that("Match expected number of screens in lesion state", {
  dx <- "L"
  # With and without surveillance
  for (survstrat in c("nosurv", "surv")) {
    # 1 - gold standard test, 2 - screening test
    for (strat in 1:2) {
      # Create copy of data
      m_patient_screen <- copy(m_patient_base)
      m_lesion_screen <- copy(m_lesion_base)
      
      # Combine with lesion data
      m_cohort <- list(patient_level = m_patient_screen,
                       lesion_level = m_lesion_screen)
      
      # Set surveillance parameters
      if (survstrat == "nosurv") {
        params_surv <- NULL
      } else {
        params_surv <- l_params_screen$surveil
      }
      
      # Run screening with no surveillance
      run_screening_counterfactual(m_cohort,
                                   l_params_model = l_params_calib$l_params_model,
                                   l_params_strategy = l_params_screen$strats[[strat]],
                                   l_params_test = l_params_screen$test_chars,
                                   l_params_surveil = params_surv)
      
      # Screening test count
      if (strat == 1) {
        expect_equal("ct_tests_screen" %in% names(m_patient_screen), FALSE)
      } else {
        expect_equal(m_patient_screen[, get(paste("ct", strat, survstrat, "s", dx, sep = "_"))], m_patient_screen[, ct_tests_screen])
      }
      
      # Confirmatory and base tests
      expect_equal(m_patient_screen[, get(paste("ct", strat, survstrat, "c", dx, sep = "_"))], m_patient_screen[, ct_tests_confirm])
      expect_equal(m_patient_screen[, get(paste("ct", strat, survstrat, "b", dx, sep = "_"))], m_patient_screen[, ct_tests_base])
      
      # Detection time
      expect_equal(m_cohort$lesion_level[[paste("test", strat, survstrat, "time_detected", sep = "_")]], m_cohort$lesion_level$time_detected)
    }
  }
})

# Test variation of screening
test_that("Match expected number of screens with variation", {
  # Stack copies of patient and lesion data
  m_patient_rep <- do.call(rbind, replicate(nreps_var, m_patient_base[,.SD, .SDcols = !patterns("ct_|purpose")], simplify=FALSE))
  m_lesion_rep <- do.call(rbind, replicate(nreps_var, m_lesion_base[,.SD, .SDcols = !patterns("ct_|purpose")], simplify=FALSE))
  
  # Save original patient ID and update patient IDs to be unique
  m_patient_rep$pt_id_orig <- m_patient_rep$pt_id
  m_lesion_rep$pt_id_orig <- m_lesion_rep$pt_id
  m_patient_rep$pt_id <- m_patient_rep$pt_id + nrow(m_patient_base)*rep(seq(0, nreps_var-1), each = nrow(m_patient_base))
  m_lesion_rep$pt_id <- m_lesion_rep$pt_id + nrow(m_patient_base)*rep(seq(0, nreps_var-1), each = nrow(m_lesion_base))

  # Set patient ID as key
  setkey(m_patient_rep, pt_id)
  setkey(m_lesion_rep, pt_id)
  
  # # manual
  # dx <- "L"
  # survstrat <- "surv"
  # strat <- 1
  
  # For cancer with and without lesion state
  for (dx in c("L", "P")) {
    # Set model parameters
    if (dx == "L") {
      l_params_model_temp <- l_params_calib$l_params_model
    } else {
      l_params_model_temp <- modifyList(l_params_calib$l_params_model, list(v_states = c("H", "P", "C", "D")))
    }
    
    # With and without surveillance
    for (survstrat in c("nosurv", "surv")) {
      # Skip surveillance for preclinical
      if (dx == "P" & survstrat == "surv") next
      
      # 1 - gold standard test, 2 - screening test
      for (strat in 1:2) {
        # Create copy of data
        m_patient_screen <- copy(m_patient_rep)
        m_lesion_screen <- copy(m_lesion_rep)
        
        # Combine with lesion data
        m_cohort <- list(patient_level = m_patient_screen,
                         lesion_level = m_lesion_screen)
        
        # Set surveillance parameters
        if (survstrat == "nosurv") {
          params_surv <- NULL
        } else {
          params_surv <- l_params_screen$surveil
        }
        
        # Run screening with no surveillance
        run_screening_counterfactual(m_cohort,
                                     l_params_model = l_params_model_temp,
                                     l_params_strategy = l_params_screen$strats[[strat]],
                                     l_params_test = params_screen$test_chars, # Original distributions for screening parameters
                                     l_params_surveil = params_surv)
        
        # Consolidate results at original patient ID level
        ct_cols <- grep("ct_tests_", names(m_patient_screen), value = TRUE)
        m_patient_summ <- m_patient_screen[, {
          means <- lapply(.SD, function(x) mean(x, na.rm = TRUE))
          sds   <- lapply(.SD, function(x) sd(x, na.rm = TRUE))
          setNames(c(means, sds), c(paste0("mean_", ct_cols), paste0("sd_", ct_cols)))
        }, by = pt_id_orig, .SDcols = ct_cols]
        
        # Set label for variable in Excel
        if (dx == "L") {
          label_var <- paste(strat, survstrat, sep = "_")
        } else {
          label_var <- strat
        }
        
        # Screening test count
        # if (strat == 1) {
        #   expect_equal("ct_tests_screen" %in% names(m_patient_screen), FALSE)
        # } else {
        #   expect_equal(m_patient_screen[, get(paste("ct", label_var, "s", dx, sep = "_"))], m_patient_screen[, ct_tests_screen])
        # }
        # 
        # # Confirmatory and base tests
        # expect_equal(m_patient_screen[, get(paste("ct", label_var, "c", dx, sep = "_"))], m_patient_screen[, ct_tests_confirm])
        # expect_equal(m_patient_screen[, get(paste("ct", label_var, "b", dx, sep = "_"))], m_patient_screen[, ct_tests_base])
        # 
        # # Detection time
        # expect_equal(m_cohort$lesion_level[[paste("test", label_var, "time_detected", sep = "_")]], m_cohort$lesion_level$time_detected)
      }
    }
  }
})
