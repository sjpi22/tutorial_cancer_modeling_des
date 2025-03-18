#' Generate time to event data for base case (no screening) decision model
#'
#' @param l_params_model List with all parameters of decision model
#' 
#' @return Data table of time to event data for each individual
#' 
#' @export
run_base_model <- function(l_params_model) {
  # Set seed
  if(!is.null(l_params_model$seed)) {
    set.seed(l_params_model$seed)
  }
  
  # Initialize matrix of patient data
  m_patients <- initialize_cohort(l_params_model)
  
  # Simulate baseline characteristics
  simulate_baseline_data(m_patients, l_params_model)
  
  # Simulate disease (lesion and/or cancer) onset
  simulate_disease_onset(m_patients, l_params_model)
  
  # Simulate additional lesion onset and lesion progression
  if("L" %in% l_params_model$v_states) {
    m_lesions <- simulate_additional_lesions(m_patients, l_params_model)
    simulate_lesion_progression(m_patients, m_lesions, l_params_model)
  }
  
  # Simulate cancer progression
  simulate_cancer_progression(m_patients, l_params_model)
  
  # Simulate cancer mortality
  simulate_cancer_mortality(m_patients, l_params_model)
  
  # Compile overall mortality outcomes from background and cancer data
  calc_mortality_outcomes(m_patients)
  
  # Count confirmatory tests for clinical cancer cases if flag is not NULL
  if (!is.null(l_params_model$fl_count_tests)) {
    m_patients[time_H_C < time_H_D, ct_tests_conf := 1]
  }
  
  # Wrap no screening results
  if("L" %in% l_params_model$v_states) {
    res <- list(patient_level = m_patients,
                lesion_level = m_lesions)
  } else {
    res <- m_patients
  }
  
  return(res)
}


# Initialize time to event matrix for patient cohort with patient IDs and baseline strategy
initialize_cohort <- function(l_params_model) {
  m_patients_init <- with(as.list(l_params_model), {
    m_patients <- data.table(
      pt_id = 1:n_cohort
    )
    setkey(m_patients, pt_id)
    return(m_patients)
  }
  )
  return(m_patients_init)
}


# Simulate baseline characteristics
simulate_baseline_data <- function(m_patients, l_params_model) {
  with(l_params_model, {
    # Sample time to death from other causes for multiple sexes
    if (length(sex) > 1) {
      # Sample male / female
      m_patients[, male := query_distr(
        "r", .N, 
        d_male$distr, 
        d_male$params
      )]
      
      # Sample time to death from other causes if male
      m_patients[male == 1, time_H_Do := query_distr(
        "r", .N, 
        d_time_H_Do[["male"]]$distr, 
        d_time_H_Do[["male"]]$params
      )]
      
      # Sample time to death from other causes if female
      m_patients[male == 0, time_H_Do := query_distr(
        "r", .N, 
        d_time_H_Do[["female"]]$distr, 
        d_time_H_Do[["female"]]$params
      )]
    } else { # Sample time to death from other causes for single sex case
      # Assign male / female
      m_patients[, male := (sex == "male")]
      
      # Sample time to death from other causes if female
      m_patients[, time_H_Do := query_distr(
        "r", .N, 
        d_time_H_Do[[sex]]$distr, 
        d_time_H_Do[[sex]]$params
      )]
    }
  })
  return(NULL)
}


# Simulate disease (lesion or cancer) onset
simulate_disease_onset <- function(m_patients, l_params_model) {
  with(l_params_model, {
    # If including lesion state, simulate time to lesion and cancer onset
    if("L" %in% v_states) {
      # Simulate time to first lesion onset
      m_patients[, time_H_L := query_distr(
        "r", .N, 
        d_time_H_L$distr, 
        d_time_H_L$params
      )]
      
      # Simulate number of additional lesions (only for individuals with lesion onset before death)
      # Scale rate down by mean screening duration, then stretch by duration until death from other causes
      m_patients[time_H_Do > time_H_L, n_L_add := query_distr(
        "r", .N, 
        d_n_L$distr, 
        lapply(d_n_L$params, '*', (time_H_Do - time_H_L))
      )]
      
      # Cap number of additional lesions
      m_patients[n_L_add >= n_lesions_max, n_L_add := n_lesions_max - 1]
    } else {
      # Simulate time to de novo cancer onset
      m_patients[, time_H_P := query_distr("r", .N, d_time_H_P$distr, d_time_H_P$params)]
    }
  })
  return(NULL)
}


# Simulate development of additional lesions and output lesion-level data table
simulate_additional_lesions <- function(m_patients, l_params_model) {
  m_lesions <- with(as.list(l_params_model), {
    # Subset to individuals with lesion onset before death
    m_lesions <- m_patients[time_H_L < time_H_Do]
    
    # Account for case with no lesions
    if (nrow(m_lesions) > 0) {
      # Set initial lesion ID and individual onset time variable
      m_lesions[, `:=` (lesion_id = 1, time_L_Lj = 0)]
      
      # Create new row for each additional lesion (account for case where no one has >1 lesion)
      if (nrow(m_lesions[n_L_add > 0]) > 0) {
        m_other_lesion <- m_lesions[n_L_add > 0, .(
          time_H_Do,
          time_H_L = time_H_L,
          lesion_id = seq(2, n_L_add + 1)
        ), 
        pt_id]
        
        # Simulate onset time of each additional lesion
        m_other_lesion[, time_L_Lj := query_distr(
          "r", .N, 
          d_time_L_Lj$distr, 
          lapply(d_time_L_Lj$params, '*', (time_H_Do - time_H_L))
        )]
        
        # Append original lesion data
        m_lesions <- rbindlist(list(
          m_lesions[, .SD, .SDcols = intersect(names(m_lesions), names(m_other_lesion))],
          m_other_lesion
        ))
      }
      
      # Calculate time from birth to lesion-specific onset
      m_lesions[, time_H_Lj := time_H_L + time_L_Lj]
      
      # Set pt_id as key
      setkey(m_lesions, pt_id)
    }
    return(m_lesions)
  })
  return(m_lesions)
}


# Simulate lesion progression to cancer
simulate_lesion_progression <- function(m_patients, m_lesions, l_params_model) {
  with(as.list(l_params_model), {
    # Account for case of no lesions
    if (!is.null(m_lesions)) {
      # Sample time to cancer
      m_lesions[, time_Lj_Pj := query_distr(
        "r", .N, 
        d_time_L_P$distr, 
        d_time_L_P$params)]
      
      # Calculate time from first lesion to lesion-specific cancer onset
      m_lesions[, time_L_Pj := time_L_Lj + time_Lj_Pj]
      
      # Calculate time from healthy to lesion-specific cancer onset
      m_lesions[, time_H_Pj := time_H_L + time_L_Pj]
      
      # Get earliest time to cancer onset
      m_patients_cancer <- m_lesions[, .(time_H_P = min(time_H_Pj)), by = pt_id]
      
      # Merge to patient and lesion data
      m_patients[m_patients_cancer, time_H_P := i.time_H_P]
      m_lesions[m_patients_cancer, time_H_P := i.time_H_P]
    } else {
      m_patients[, time_H_P := Inf]
    }
  })
  return(NULL)
}


# Simulate cancer stage progression by stage
simulate_cancer_progression <- function(m_patients, l_params_model) {
  with(l_params_model, {
    # Populate time to progression to next preclinical stage
    for (i in 1:(length(v_cancer)-1)) {
      var_progress <- paste0("time_P", v_cancer[i], "_P", v_cancer[i+1])
      m_patients[time_H_P < time_H_Do, (var_progress) := query_distr(
        "r", .N, 
        get(paste0("d_", var_progress))$distr, 
        get(paste0("d_", var_progress))$params
      )]
    }
    
    # Populate time to detection (clinical cancer)
    for (i in 1:length(v_cancer)) {
      var_detect <- paste0("time_P", v_cancer[i], "_C", v_cancer[i])
      m_patients[time_H_P < time_H_Do, (var_detect) := query_distr(
        "r", .N, 
        get(paste0("d_", var_detect))$distr, 
        get(paste0("d_", var_detect))$params
      )]
    }
    
    # Initialize stage of diagnosis variable
    m_patients[, stage_dx := eval(parse(text = paste0("as.", class(v_cancer), "(NA)")))]
    
    # Initialize time from cancer onset to detection (should be NA for patients without cancer onset before death from other causes)
    m_patients[time_H_P < time_H_Do, time_P_C := 0]
    
    # Loop through stages until detection occurs
    for (i in 1:(length(v_cancer)-1)) {
      # Get progression and detection variables
      var_progress <- paste0("time_P", v_cancer[i], "_P", v_cancer[i+1])
      var_detect <- paste0("time_P", v_cancer[i], "_C", v_cancer[i])
      
      # If stage at diagnosis has not been set yet and detection occurs before 
      # progression, add time to detection to running total of time from cancer 
      # onset to detection and set stage at diagnosis
      m_patients[time_H_P < time_H_Do & get(var_detect) < get(var_progress) & is.na(stage_dx), `:=` (
        stage_dx = v_cancer[i],
        time_P_C = time_P_C + get(var_detect))]
      
      # Otherwise if not detected yet, add time to detection to running total 
      # of time from cancer onset to detection
      m_patients[time_H_P < time_H_Do & is.na(stage_dx), 
                 time_P_C := time_P_C + get(var_progress)]
    }
    
    # If not detected yet by second to last stage, set last stage as stage at detection
    var_detect <- paste0("time_P", tail(v_cancer, 1), "_C", tail(v_cancer, 1))
    m_patients[time_H_P < time_H_Do & is.na(stage_dx), `:=` (
      stage_dx = tail(v_cancer, 1),
      time_P_C = time_P_C + get(var_detect))]
  })
  
  # Set time to cancer diagnosis
  m_patients[!is.na(time_P_C), time_H_C := time_H_P + time_P_C]
  return(NULL)
}


# Simulate cancer mortality by stage at diagnosis among people with cancer onset before death from other causes
simulate_cancer_mortality <- function(m_patients, l_params_model) {
  # Get indices of patients with cancer onset before death from other causes
  # Not cancer diagnosis before death from other causes to be more conservative, in case screening may convert from preclinical to clinical before death
  idx <- m_patients[time_H_P < time_H_Do, which = TRUE]
  
  # Get progression and detection variables
  if (length(idx) > 0) {
    with(l_params_model, {
      m_patients[idx, time_C_Dc := query_distr(
        "r", .N, 
        get(paste0("d_time_C", stage_dx, "_Dc"))$distr, 
        get(paste0("d_time_C", stage_dx, "_Dc"))$params
      ), by = stage_dx]
    })
    
    # Calculate death from cancer
    m_patients[idx, time_H_Dc := time_H_C + time_C_Dc]
  } else {
    m_patients[, `:=`(time_C_Dc = NA,
                      time_H_Dc = NA)]
  }
  return(NULL)
}


# Generate mortality outcomes
calc_mortality_outcomes <- function(m_patients) {
  # Join cancer data to patient-level data
  # Calculate all-cause death and cause of death
  m_patients[, time_H_D := pmin(time_H_Do, time_H_Dc, na.rm = TRUE)]
  m_patients[, fl_Dc := (time_H_Do > pmin(time_H_Dc, Inf, na.rm = TRUE))]
  
  # Calculate survival from cancer diagnosis
  m_patients[time_H_C <= time_H_D, time_C_D := time_H_D - time_H_C]
  return(NULL)
}


# Rerun version for screening (note: overwrites input data by reference)
run_screening_counterfactual <- function(
    m_cohort,
    l_params_model,
    l_params_strategy,
    l_params_tests,
    verbose = FALSE
) {
  # Separate patient and lesion data as necessary
  if (is.data.table(m_cohort)) {
    m_patients <- m_cohort
  } else {
    m_patients <- m_cohort$patient_level
    m_lesions <- m_cohort$lesion_level
  }
  
  # Extract screening and confirmatory test modalities
  v_mods <- l_params_strategy[["mod"]]
  if (!is.null(l_params_strategy[["conf_test"]])) {
    v_mods <- c(v_mods, l_params_strategy[["conf_test"]])
  }
  
  # Sample sensitivity and specificity of screening and (if applicable) confirmatory tests
  l_params_tests_sample <- list()
  for (mod in v_mods) {
    l_params_tests_sample[[mod]] <- with(l_params_tests[[mod]], {
      list(
        p_sens = lapply(d_p_sens, 
                        function (u) query_distr("r", 1, u$distr, u$params)),
        p_spec = query_distr("r", 1, d_p_spec$distr, d_p_spec$params)
      )
    })
  }
  
  #### Case 1: Simulate screening before disease ####
  # Set variable for time to disease onset
  var_onset <- paste("time", l_params_model$v_states[1], l_params_model$v_states[2], sep = "_")
  
  # Simulate screening during healthy state
  simulate_screening_H(m_patients = m_patients,
                       var_onset = var_onset,
                       l_params_strategy = l_params_strategy,
                       l_test_chars = l_params_tests_sample,
                       verbose = verbose)
  
  #### Case 2: Lesions developed by screen age, but no preclinical cancer ####
  if ('L' %in% l_params_model$v_states) {
    # Simulate lesion progression to cancer
    simulate_screening_L(m_patients = m_patients, 
                         m_lesions = m_lesions,
                         l_params_strategy = l_params_strategy,
                         l_test_chars = l_params_tests_sample,
                         verbose = verbose)
  }
  
  #### Case 3: Preclinical cancer screening and surveillance ####
  simulate_screening_P(m_patients = m_patients,
                       l_params_model = l_params_model,
                       l_params_strategy = l_params_strategy,
                       l_test_chars = l_params_tests_sample,
                       verbose = verbose)
  
  # Recalculate mortality outcomes
  calc_mortality_outcomes(m_patients)
  
  return(NULL)
}


# Simulate screening occurring during the healthy state
simulate_screening_H <- function(m_patients,
                                 var_onset,
                                 l_params_strategy,
                                 l_test_chars,
                                 verbose = FALSE
) {
  # Extract screening and confirmatory test modalities
  mod <- l_params_strategy[["mod"]]
  mod_conf <- l_params_strategy[["conf_test"]]
  
  # Extract screening test specificity and interval
  p_spec <- l_test_chars[[mod]][["p_spec"]]
  int_screen <- l_params_strategy[["int_screen"]]
  
  # Extract confirmatory test parameters if applicable
  if (!is.null(mod_conf)) {
    # Extract confirmatory test interval
    int_conf <- l_params_strategy[["int_conf"]]
    
    # If confirmatory test interval is not specified, set equal to screening interval
    if (is.null(int_conf)) {
      int_conf <- int_screen
    }
  } else { # Otherwise set NA
    int_conf <- NA
  }
  
  # Calculate age at which individual no longer screens as minimum of patient censor age and strategy stop age
  # and reset confirmatory tests to 0
  m_patients[, `:=` (time_screen_stop = pmin(time_screen_censor, l_params_strategy$age_screen_stop),
                     ct_tests_conf = 0)]
  
  # If screening interval is constant, perform simplified calculation of number of screening tests during healthy state
  # Otherwise, loop over testing rounds and assign applicable screening interval
  if (is.null(mod_conf) | int_conf %in% int_screen) {
    # Get routine screening ages
    v_screen_ages <- with(l_params_strategy, seq(age_screen_start, age_screen_stop, int_screen))
    
    # Get number of screening tests during healthy state before earliest of censor age or disease onset
    m_patients[, ct_tests_screen := findInterval(pmin(time_screen_stop, get(var_onset)), v_screen_ages, left.open = T)]
    
    # Simulate number of false positives screens
    m_patients[, `:=` (
      ct_tests_screen_FP = rbinom(
        .N,
        size = ct_tests_screen,
        prob = 1 - p_spec
      ))]
    
    # If screening test is followed by confirmatory test, set confirmatory tests equal to false positive screening tests
    if (!is.null(mod_conf)) {
      m_patients[, ct_tests_conf := ct_tests_screen_FP]
    }
    
    # Get first screening age with disease present
    m_patients[, screen_age := v_screen_ages[ct_tests_screen + 1], by = ct_tests_screen]
    
    # Reset first screening age to NA if after end date
    m_patients[screen_age >= time_screen_stop, screen_age := NA]
  } else { 
    # Initialize test counts to 0
    m_patients[, `:=` (ct_tests_screen = 0,
                       ct_tests_screen_FP = 0)]
    
    # Initialize running screening age at screening start age for individuals eligible to receive screening at that age
    m_patients[l_params_strategy$age_screen_start < time_screen_censor, screen_age := l_params_strategy$age_screen_start]
    
    # Calculate number of healthy patients remaining to screen
    if (verbose) print("Number of healthy individuals remaining to be screened:")
    n_healthy_screen <- m_patients[screen_age < pmin(get(var_onset), time_screen_stop, na.rm = T), .N]
    
    # Run while loop
    while (n_healthy_screen > 0) {
      if (verbose) print(n_healthy_screen)
      # For individuals in healthy state, increment number of screening tests and sample whether test produces a false positive result
      m_patients[screen_age < pmin(get(var_onset), time_screen_stop, na.rm = T), `:=` (
        ct_tests_screen = ct_tests_screen + 1,
        fl_positive = rbinom(
          .N,
          size = 1,
          prob = 1 - p_spec
        ))]
      
      # For false positives, increment number of false positive screening tests
      m_patients[fl_positive == 1, `:=` (ct_tests_screen_FP = ct_tests_screen_FP + 1)]
      
      # Set next screening interval and increment confirmatory tests depending on type of test used
      if (is.null(mod_conf)) { 
        m_patients[!is.na(fl_positive), `:=` (int_test_next = int_screen)]
      } else {
        # Increment confirmatory tests and set next test to confirmatory test interval for false positives
        m_patients[fl_positive == 1, `:=` (ct_tests_conf = ct_tests_conf + 1,
                                           int_test_next = int_conf)]
        
        # Set next test to screening interval for true negatives
        m_patients[fl_positive == 0, `:=` (int_test_next = int_screen)]
      }
      
      # Update screen age, reset to NA if after stop date, and reset temporary variables
      m_patients[!is.na(fl_positive), screen_age := screen_age + int_test_next]
      m_patients[screen_age >= time_screen_stop, screen_age := NA]
      m_patients[, `:=` (fl_positive = NA,
                         int_test_next = NA)]
      
      # Recalculate number of healthy patients remaining to screen
      n_healthy_screen <- m_patients[screen_age < pmin(get(var_onset), time_screen_stop, na.rm = T), .N]
    }
  }
  
  return(NULL)
}


# Simulate screening occurring during the precancerous lesion state
# Assumes that screening during healthy state has been completed for population,
# which initializes the screen_age variable
simulate_screening_L <- function(m_patients,
                                 m_lesions,
                                 var_onset,
                                 l_params_strategy,
                                 l_test_chars,
                                 verbose = FALSE
) {  
  # Extract screening and confirmatory test modalities
  mod <- l_params_strategy[["mod"]]
  mod_conf <- l_params_strategy[["conf_test"]]
  
  # Extract screening test sensitivity, specificity, and interval
  p_sens <- l_test_chars[[mod]][["p_sens"]][["L"]]
  p_spec <- l_test_chars[[mod]][["p_spec"]]
  int_screen <- l_params_strategy[["int_screen"]]
  
  # Extract confirmatory test parameters if applicable
  if (!is.null(mod_conf)) {
    p_sens_conf <- l_test_chars[[mod_conf]][["p_sens"]][["L"]]
    p_spec_conf <- l_test_chars[[mod_conf]][["p_spec"]]
    
    if (!is.null(mod_conf)) {
      # Extract confirmatory test interval
      int_conf <- l_params_strategy[["int_conf"]]
      
      # If confirmatory test interval is not specified, set equal to screening interval
      if (is.null(int_conf)) {
        int_conf <- int_screen
      }
    } else { # Otherwise set NA
      int_conf <- NA
    }
  }
  
  # Merge first screen age within lesion state and minimum of cancer 
  # onset and screening stop date as lesion screening censor date (to be updated after every screening)
  setkey(m_lesions, pt_id)
  m_lesions[m_patients[screen_age >= time_H_L & 
                         screen_age < pmin(time_H_P, time_screen_stop, na.rm = T)], `:=` (
                           time_lesion_censor = pmin(i.time_H_P, i.time_screen_stop, na.rm = T),
                           screen_age = i.screen_age)]
  
  # Among patients who will receive screening during the lesion state,
  # initialize flag for whether lesion was removed
  m_lesions[!is.na(screen_age), `:=` (fl_removed = 0)]
  
  # Get number of lesions available for screening
  n_lesion_screen <- m_lesions[!is.na(screen_age), .N]
  
  # Loop through screening tests until there are no more lesions to screen
  if (verbose) print("Number of lesions remaining to screen:")
  while (n_lesion_screen > 0) {
    if (verbose) print(n_lesion_screen)
    
    ###### 2.0 Flag eligible screeners and lesions
    # Initialize flags for eligible screeners and lesions
    m_lesions[, `:=` (fl_screen = NA,
                      fl_present = NA)]
    
    # Flag individuals eligible for screening at current screening time
    m_lesions[screen_age < time_lesion_censor, fl_screen := 1]
    
    # Flag lesions that are present in eligible screeners at current screen time
    m_lesions[fl_screen == 1 & time_H_Lj <= screen_age & fl_removed == 0,
              fl_present := 1]
    
    ###### 2.1 Flag whether eligible lesions would be detected
    # Sample whether eligible lesions produce positive test result
    m_lesions[fl_present == 1, `:=` (
      fl_positive = rbinom(
        .N,
        size = 1,
        prob = p_sens
      ))]
    
    # Apply downstream effect of positive result
    if (is.null(mod_conf)) {
      # If the screening test is also the confirmatory test, set positive lesions as detected/removed
      m_lesions[fl_present == 1 & fl_positive == 1, `:=` (fl_removed = 1)]
    } else {
      # If there is a follow-on confirmatory test, apply confirmatory tests 
      # based on whether test is applied to only positive lesions (targeted is TRUE) or 
      # all present lesions if there is at least one positive lesion (targeted is FALSE)
      if (l_test_chars[[mod]][["targeted"]]) { # Apply confirmatory tests to positive lesions only
        # Sample whether positive lesions are removed
        m_lesions[fl_present == 1 & fl_positive == 1, `:=` (
          fl_removed = rbinom(
            .N,
            size = 1,
            prob = p_sens_conf
          ))]
      } else { # Apply confirmatory tests to all lesions if at least one produces positive result
        # Check whether any lesion led to positive screening test
        m_detect <- m_lesion[fl_present == 1,
                             .(fl_positive_any = max(fl_positive)),
                             by = pt_id]
        
        # Merge positive screening test flag
        m_lesion[m_detect, fl_positive_any := i.fl_positive_any]
        
        # Sample removal with confirmatory test
        m_lesion[fl_present == 1 & fl_positive_any == 1,`:=` (
          fl_removed = rbinom(
            .N,
            size = 1,
            prob = p_sens_conf
          ))]
      }
    }
    
    # For removed lesions, set time to cancer onset to Inf, and set detection time
    m_lesions[fl_present == 1 & fl_removed == 1, `:=` (
      time_H_Pj = Inf,
      time_detected = screen_age)]
    
    # Check number of eligible lesions and removed lesions among individuals screened at this round
    # and update time to onset of preclinical cancer
    m_lesions_summary <- m_lesions[fl_screen == 1, 
                                   .(ct_eligible = sum(fl_present, na.rm = T), # Number of lesions present
                                     fl_positive = max(fl_positive, na.rm = T), # Whether any lesions produced a positive screen result
                                     time_H_P = min(time_H_Pj, na.rm = T)), 
                                   by = pt_id]
    
    ###### 2.2 Flag false positives among screened individuals with no active lesions
    m_lesions_summary[ct_eligible == 0, fl_FP := rbinom(
      .N,
      size = 1,
      prob = 1 - p_spec)]
    
    ###### 2.3 Update screening data
    # Process positive and confirmatory test flags for merging
    m_lesions_summary[ct_eligible == 0, fl_positive := fl_FP]
    if (!is.null(mod_conf)) {
      m_lesions_summary[, fl_test_conf := fl_positive]
    } else {
      m_lesions_summary[, fl_test_conf := 0]
    }
    
    # Set next screen time
    if (is.null(mod_conf)) {
      m_lesions_summary[, int_screen_next := int_screen]
    } else {
      m_lesions_summary[, int_screen_next := ifelse(fl_positive %in% 1, int_conf, int_screen)]
    }
    
    # Update time to preclinical cancer and following times, increment number of screens, number of confirmatory tests, and screen age
    m_patients[m_lesions_summary, `:=` (time_H_P = i.time_H_P,
                                        ct_tests_screen = ct_tests_screen + 1,
                                        ct_tests_screen_FP = ct_tests_screen_FP + pmax(0, i.fl_FP, na.rm = T),
                                        ct_tests_conf = ct_tests_conf + pmax(0, i.fl_test_conf, na.rm = T),
                                        screen_age = screen_age + int_screen_next)]
    
    # Reset screening age to NA if after censor date or end age
    m_patients[screen_age >= time_screen_stop, screen_age := NA]
    
    # Merge screen age within lesion state and minimum of cancer 
    # onset and death as lesion screening censor date (to be updated after every screening)
    m_lesions[, `:=` (screen_age = NA)]
    m_lesions[m_patients[screen_age >= time_H_L & 
                           screen_age < pmin(time_H_P, time_screen_stop, na.rm = T)], `:=` (
                             time_H_P = i.time_H_P,
                             time_lesion_censor = pmin(i.time_H_P, i.time_screen_stop, na.rm = T),
                             screen_age = i.screen_age)]
    
    # Update number of lesions available for screening
    n_lesion_screen <- m_lesions[!is.na(screen_age), .N]
  }
  
  # Update mortality outcomes
  m_patients[, `:=` (time_H_C = time_H_P + time_P_C,
                     time_H_Dc = time_H_P + time_P_C + time_C_Dc)]
  calc_mortality_outcomes(m_patients)
  return(NULL)
}


# Simulate screening occurring during the preclinical cancer state
# Assumes that screening during healthy state has been completed for population,
# which initializes the screen_age variable
simulate_screening_P <- function(m_patients,
                                 l_params_model,
                                 l_params_strategy,
                                 l_test_chars,
                                 verbose = FALSE
) {
  # Extract screening and confirmatory test modalities
  mod <- l_params_strategy[["mod"]]
  mod_conf <- l_params_strategy[["conf_test"]]
  
  # Extract screening test sensitivity and interval
  p_sens <- l_test_chars[[mod]][["p_sens"]][["P"]]
  int_screen <- l_params_strategy[["int_screen"]]
  
  # Extract confirmatory test parameters if applicable
  if (!is.null(mod_conf)) {
    # Extract test sensitivity
    p_sens_conf <- l_test_chars[[mod_conf]][["p_sens"]][["P"]]
    
    # Extract confirmatory test interval
    int_conf <- l_params_strategy[["int_conf"]]
    
    # If confirmatory test interval is not specified, set equal to screening interval
    if (is.null(int_conf)) {
      int_conf <- int_screen
    }
  } else { # Otherwise set NA
    inf_conf <- NA
  }
  
  # Run while loop
  if (verbose) print("Number of individuals with preclinical cancer remaining to screen:")
  while (m_patients[!is.na(screen_age), .N] > 0) {
    if (verbose) print(m_patients[!is.na(screen_age), .N])
    # Increment number of screening tests and sample whether cancer leads to positive screen
    m_patients[!is.na(screen_age), `:=` (
      ct_tests_screen = ct_tests_screen + 1,
      fl_positive = rbinom(
        .N,
        size = 1,
        prob = p_sens
      ))]
    
    # Apply detection after positive tests
    if (is.null(mod_conf)) {
      # If positive screening tests are also the confirmatory test, flag positive cases as detected
      m_patients[fl_positive == 1, `:=` (fl_detected = 1)]
    } else { 
      # Otherwise, increment confirmatory tests, sample whether cancer is detected, and set next test to confirmatory test interval for positive screens
      m_patients[fl_positive == 1, `:=` (
        ct_tests_conf = ct_tests_conf + 1,
        fl_detected = rbinom(
          .N,
          size = 1,
          prob = p_sens_conf
        ),
        int_test_next = int_conf)]
    }
    
    # Set next test to screening interval for negatives
    m_patients[fl_positive == 0, `:=` (int_test_next = int_screen)]
    
    # Reset clinical cancer age and set next screen age to NA for newly detected cases
    m_patients[!is.na(screen_age) & fl_detected == 1, `:=` (
      time_H_C = screen_age,
      time_P_C = screen_age - time_H_P,
      screen_age = NA)]
    
    # Update screen age, reset to NA if after stop date, and reset temporary variables
    m_patients[!is.na(screen_age), screen_age := screen_age + int_test_next]
    m_patients[screen_age >= time_screen_stop, screen_age := NA]
    m_patients[, `:=` (fl_positive = NA,
                       int_test_next = NA)]
  }
  
  # For patients with detected cancer, recalculate stage at diagnosis
  # Note: Assumes at least 2 stages of cancer
  m_patients[fl_detected == 1, `:=` (stage_dx = 1,
                                     time_P_C_running = time_P1_P2)]
  
  for (stg in 2:(length(l_params_model$v_cancer) - 1)) {
    m_patients[fl_detected == 1 & time_P_C > time_P_C_running, `:=` (
      stage_dx = stage_dx + 1,
      time_P_C_running = time_P_C_running + get(paste0("time_P", stg, "_P", stg + 1)))]
  }
  m_patients[fl_detected == 1 & time_P_C > time_P_C_running, stage_dx := stage_dx + 1]
  
  # Recalculate time to death from cancer among people with screen-detected cancer
  with(l_params_model, {
    m_patients[fl_detected == 1, time_C_Dc := query_distr(
      "r", .N,
      get(paste0("d_time_C", stage_dx, "_Dc"))$distr,
      get(paste0("d_time_C", stage_dx, "_Dc"))$params
    ), by = stage_dx]
  })
  
  # Flag cancers diagnosed with symptoms and increment confirmatory test count
  m_patients[!fl_detected %in% 1 & time_H_C < time_H_D, `:=` (
    fl_symptom_detected = 1,
    ct_tests_conf = ct_tests_conf + 1)]
  
  return(NULL)
}


# Calculate life years gained (LYG) from screening
calc_lyg <- function(res_base, res_screen, unit = 1) {
  lyg <- (res_screen["time_total"] - res_base["time_total"]) / res_base["N"] * unit
  return(lyg)
}


# Calculate confirmatory tests in base case scenario
calc_ntests <- function(m_patients,
                        censor_var = "time_screen_censor",
                        age_min = 0,
                        test_var_pattern = "ct_tests_"
) {
  # Calculate confirmatory tests
  if (age_min > 0 | !is.null(age_min)) {
    ct_tests <- m_patients[get(censor_var) > age_min, colSums(.SD, na.rm = T), .SDcols = patterns("ct_")]
  } else {
    ct_tests <- m_patients[, colSums(.SD, na.rm = T), .SDcols = patterns(test_var_pattern)]
  }
  return(ct_tests)
}


# Generate outputs from individual-level cohort data
# l_params_outcome: list of outcome parameters
# l_censor_vars: list of lists of variables to use to create new variables for censoring outcomes
calc_cohort_outputs <- function(m_cohort,
                                l_params_outcome,
                                l_censor_vars = NULL
) {
  # Separate patient and lesion data as necessary
  if (is.data.table(m_cohort)) {
    m_patients <- m_cohort
  } else {
    m_patients <- m_cohort$patient_level
    m_lesions <- m_cohort$lesion_level
  }
  
  # Create censor variables
  if (!is.null(l_censor_vars)) {
    for (dt in names(l_censor_vars)) {
      for (varname in names(l_censor_vars[[dt]])) {
        get(dt)[, (varname) := do.call("pmin", c(.SD, na.rm = TRUE)),
                .SDcols = l_censor_vars[[dt]][[varname]]]
      }
    }
  }
  
  # Calculate outcomes
  l_results <- list()
  for (outcome in names(l_params_outcome)) {
    l_results[[outcome]] <- do.call(
      paste0("calc_", l_params_outcome[[outcome]][["outcome_type"]]), 
      c(lapply(l_params_outcome[[outcome]][["get_params"]], get, envir = sys.frame(sys.parent(0))), 
        l_params_outcome[[outcome]][["lit_params"]]))
  }
  
  # Return outputs
  return(l_results)
}


# Reshape outputs to single vector
reshape_outputs <- function(l_outputs, var_to_keep = "value") {
  v_outputs <- c()
  for (df in l_outputs) {
    if (var_to_keep %in% names(df)) {
      if (is.data.frame(df)) { # Extract data frame column
        v_outputs <- c(v_outputs, df[[var_to_keep]])
      } else { # Extract vector element
        v_outputs <- c(v_outputs, df[var_to_keep])
      }
    } else {
      if (is.null(dim(df))) { # Add all values of vector
        v_outputs <- c(v_outputs, df)
      } else { # Add placeholders for dataframe without column
        n_vals <- nrow(df)
        v_outputs <- c(v_outputs, rep(NA, n_vals))
      }
    }
  }
  
  return(v_outputs)
}


# Wrapper for running model and outputting vector of outputs
params_to_outputs <- function(l_params_model, 
                              v_params_update = NULL, 
                              param_map = NULL,
                              l_params_outcome, # List of base case outcomes and parameters to calculate them
                              l_params_screen = NULL, # If populated, run screening scenario with list of parameters
                              l_params_outcome_screen = NULL, # If populated, calculate different outcomes for screening scenario, otherwise use base outcomes
                              l_params_outcome_counter = NULL, # Counterfactual (base vs. screening) outcomes
                              l_censor_vars = NULL,
                              reshape_output = TRUE, 
                              individual_data = FALSE, # Output individual-level data
                              conf_level = 0.95
) {
  # Update parameters
  if (!is.null(v_params_update)) {
    if (is.null(param_map)) {
      stop("Input parameter map")
    }
    l_params_update <- update_param_from_map(l_params_model, v_params_update, param_map)
  } else {
    l_params_update <- l_params_model
  }
  
  # Run base decision model
  m_cohort <- run_base_model(l_params_update)
  
  # Get base outputs
  l_outputs <- calc_cohort_outputs(m_cohort, 
                                   l_params_outcome = l_params_outcome,
                                   l_censor_vars = l_censor_vars)
  
  # Add individual-level data to results list if necessary
  if (individual_data) {
    res <- list(m_cohort = m_cohort)
  } else {
    res <- list()
  }
  
  # Reshape outputs if necessary and add to results
  if (reshape_output) {
    l_outputs <- reshape_outputs(l_outputs)
  }
  
  # Add outputs to results with label depending on whether there will be screening results to distinguish 
  if (is.null(l_params_screen)) {
    res <- c(res, outputs = list(l_outputs))
  } else {
    res <- c(res, outputs_base = list(l_outputs))
  }
  
  # If applicable, run screening counterfactual
  if (!is.null(l_params_screen)) {
    # Set screening outcome parameters
    if (is.null(l_params_outcome_screen)) {
      l_params_outcome_screen <- l_params_outcome
    }
    
    # Initialize results list for individual-level screening data if necessary
    if (individual_data) {
      res$m_cohort_screen <- list()
    }
    
    # Loop through screening strategies
    l_outputs_screen <- list()
    for (strat in names(l_params_screen$strats)) {
      # Create copy of original data if needed (otherwise, allow data to be overwritten)
      if (individual_data == T | length(l_params_screen$strats) > 1) {
        m_cohort_screen <- copy(m_cohort)
      } else {
        m_cohort_screen <- m_cohort
      }
      
      # Generate data under screening counterfactual
      do.call(run_screening_counterfactual, 
              list(m_cohort = m_cohort_screen, 
                   l_params_model = l_params_model,
                   l_params_strategy = l_params_screen$strats[[strat]],
                   l_params_tests = l_params_screen$test_chars
              ))
      
      # Add individual-level data to results list if necessary
      if (individual_data) {
        res$m_cohort_screen[[strat]] <- m_cohort_screen
      }
      
      # Calculate screening outcomes
      l_outputs_screen[[strat]] <- calc_cohort_outputs(m_cohort_screen, 
                                                       l_params_outcome = l_params_outcome_screen,
                                                       l_censor_vars = l_censor_vars)
      
      # Calculate counterfactual comparison outcomes
      if (!is.null(l_params_outcome_counter)) {
        l_outputs_counter <- list()
        for (outcome in names(l_params_outcome_counter)) {
          # Get input outcome
          input_outcome <- l_params_outcome_counter[[outcome]][["input_outcome"]]
          
          # Calculate comparison outcomes
          l_outputs_counter[[outcome]] <- do.call(
            paste0("calc_", l_params_outcome_counter[[outcome]][["outcome_type"]]), 
            c(res_base = list(res$outputs_base[[input_outcome]]), 
              res_screen = list(l_outputs_screen[[strat]][[input_outcome]]),
              l_params_outcome_counter[[outcome]][["lit_params"]]))
        }
        # Append to strategy results list
        l_outputs_screen[[strat]] <- c(l_outputs_screen[[strat]], l_outputs_counter)
      }
      
      # Reshape outputs if necessary
      if (reshape_output) {
        for (strat in l_outputs_screen) {
          l_outputs_screen[[strat]] <- reshape_outputs(l_outputs_screen[[strat]])
        }
      }
    }
    
    # Add to results
    res <- c(res, outputs_screen = list(l_outputs_screen))
  }
  
  # Return single item or list of results if >1 items
  if (length(res) == 1) {
    return(res[[1]])
  } else {
    return(res)
  }
}

