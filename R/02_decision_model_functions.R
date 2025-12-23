################################################################################
# Natural history model functions
################################################################################

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
    m_patients[time_H_C < time_H_D, ct_tests_base := 1]
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
calc_mortality_outcomes <- function(m_patients, 
                                    idx_subset = NULL # Row indices to subset to (optional)
) {
  if (is.null(idx_subset)) {
    # Join cancer data to patient-level data
    # Calculate all-cause death and cause of death
    m_patients[, time_H_D := pmin(time_H_Do, time_H_Dc, na.rm = TRUE)]
    m_patients[, fl_Dc := (time_H_Do > pmin(time_H_Dc, Inf, na.rm = TRUE))]
    
    # Calculate survival from cancer diagnosis
    m_patients[time_H_C <= time_H_D, time_C_D := time_H_D - time_H_C]
  } else {
    # Join cancer data to patient-level data
    # Calculate all-cause death and cause of death
    m_patients[idx_subset, time_H_D := pmin(time_H_Do, time_H_Dc, na.rm = TRUE)]
    m_patients[idx_subset, fl_Dc := (time_H_Do > pmin(time_H_Dc, Inf, na.rm = TRUE))]
    
    # Calculate survival from cancer diagnosis
    idx_subset_cancer <- idx_subset[m_patients[idx_subset, which(time_H_C <= time_H_D)]]
    m_patients[idx_subset_cancer, time_C_D := time_H_D - time_H_C]
  }
  return(NULL)
}


################################################################################
# Screening model functions
################################################################################

# Rerun version for screening (note: overwrites input data by reference)
run_screening_counterfactual <- function(
    m_cohort,
    l_params_model,
    l_params_strategy,
    l_params_tests,
    reversible = TRUE,
    verbose = FALSE
) {
  # Separate patient and lesion data as necessary and save original columns
  if (is.data.table(m_cohort)) {
    m_patients <- m_cohort
    if (reversible) {
      l_data_overwritten <- list(cols_orig = copy(names(m_patients)))
    }
  } else {
    m_patients <- m_cohort$patient_level
    m_lesions <- m_cohort$lesion_level
    if (reversible) {
      l_data_overwritten <- list(cols_orig = list(
        m_patients = copy(names(m_patients)),
        m_lesions = copy(names(m_lesions))
      ))
    }
  }
  
  # Sample sensitivity and specificity of screening test
  l_params_tests_sample <- sample_test_chars(l_params_tests)
  
  #### Case 1: Simulate screening before disease ####
  # Set variable for time to disease onset
  var_onset <- paste("time", l_params_model$v_states[1], l_params_model$v_states[2], sep = "_")
  
  # Simulate screening during healthy state
  simulate_screening_H(m_patients = m_patients,
                       var_onset = var_onset,
                       l_params_strategy = l_params_strategy,
                       l_test_chars = l_params_tests_sample,
                       verbose = verbose)
  
  # Save data that may be overwritten due to screening
  if (reversible) {
    v_cols_save <- c(
      "pt_id", "time_H_P", "time_H_C", "time_H_Dc", "time_H_D", "time_P_C", "time_C_Dc", "time_C_D", "time_screen_censor", "stage_dx", "fl_Dc", "ct_tests_base"
    )
    
    if (!'L' %in% l_params_model$v_states) { # Time to cancer onset would not change if cancer is de novo
      v_cols_save <- v_cols_save[!v_cols_save %in% c("time_H_P")]
    }
    
    # Subset to patients that undergo screening beyond the healthy state
    m_patient_overwritten <- m_patients[!is.na(screen_age), .SD, .SDcols = intersect(v_cols_save, names(m_patients))]
    
    if ('L' %in% l_params_model$v_states) { 
      m_lesion_overwritten <- m_lesions[pt_id %in% m_patients[!is.na(screen_age), pt_id], .SD, .SDcols = c("pt_id", "time_H_P")]
    }
  }
  
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
  
  # Join overwritten patient and lesion data as necessary
  if (reversible) {
    if (is.data.table(m_cohort)) {
      l_data_overwritten$m_cohort <- m_patient_overwritten
    } else {
      l_data_overwritten$m_cohort <- list(
        patient_level = m_patient_overwritten,
        lesion_level = m_lesion_overwritten
      )
    }
    return(l_data_overwritten)
  } else {
    return(NULL)
  }
}


# Sample test characteristics (sensitivity and specificity)
sample_test_chars <- function(l_params_tests) {
  l_params_tests_sample <- list()
  for (mod in names(l_params_tests)) {
    l_params_tests_sample[[mod]] <- with(l_params_tests[[mod]], {
      res <- list()
      # Assign sensitivity
      if(exists("p_sens")) {
        res$p_sens <- p_sens
      } else {
        res$p_sens <- lapply(d_p_sens, 
                             function (u) query_distr("r", 1, u$distr, u$params))
      }
      
      # Assign specificity
      if(exists("p_spec")) {
        res$p_spec <- p_spec
      } else {
        res$p_spec <- query_distr("r", 1, d_p_spec$distr, d_p_spec$params)
      }
      return(res)
    })
  }
  return(l_params_tests_sample)
}


# Simulate screening occurring during the healthy state
simulate_screening_H <- function(m_patients,
                                 var_onset,
                                 l_params_strategy,
                                 l_test_chars,
                                 verbose = FALSE
) {
  # Extract test specificity and interval
  p_spec <- l_test_chars[[l_params_strategy[["mod"]]]][["p_spec"]]
  int_screen <- l_params_strategy[["int_screen"]]
  
  # Set initial censor time for cancer screening
  if (!"time_screen_censor" %in% names(m_patients)) {
    m_patients[, time_screen_censor := pmin(time_H_C, time_H_D, na.rm = T)]
  }
  
  # Initialize type of screening (routine vs. surveillance), calculate age at 
  # which individual no longer screens as minimum of patient censor age and 
  # strategy stop age
  m_patients[, `:=` (
    time_screen_stop = pmin(time_screen_censor, l_params_strategy$age_screen_stop + 1)
  )]
  
  # Get routine screening ages
  v_screen_ages <- with(l_params_strategy, seq(age_screen_start, age_screen_stop, int_screen))
  
  # Get number of screening tests during healthy state before earliest of censor age or disease onset
  m_patients[, ct_tests_screen := findInterval(pmin(time_screen_stop, get(var_onset), na.rm = T), v_screen_ages, left.open = T)]
  
  # Simulate number of (false) positive screens
  m_patients[, `:=` (
    ct_tests_positive = rbinom(
      .N,
      size = ct_tests_screen,
      prob = 1 - p_spec
    ))]
  
  # Get first screening age with disease present
  m_patients[, screen_age := v_screen_ages[ct_tests_screen + 1]]
  
  # Reset first screening age to NA if after end date
  m_patients[screen_age >= time_screen_stop, screen_age := NA]
  
  return(NULL)
}


# Loop through screening in a healthy state to calculate false positives
#
# Variables that should already be defined in m_patients: 
# screen_age (initial running age for screening),
# time_screen_stop_H (censor time for screening healthy state), 
# p_spec (specificity of screening test),
# int_screen (screening interval),
# int_conf (confirmatory test interval)
#
# Variables that will be created in m_patients:
# temp_ct_tests_screen (count of screening tests)
# temp_ct_tests_screen_FP (count of false positive screening tests, implying confirmatory test)
run_screening_H_loop <- function(m_patients,
                                 p_spec,
                                 int_screen,
                                 idx_screen = NULL, # Row indices of patients to screen (otherwise defaults to all rows)
                                 verbose = FALSE) {
  # Set all row indices if NULL
  if (is.null(idx_screen)) {
    idx_screen <- seq(nrow(m_patients))
  }
  
  # Count number of healthy patients to screen
  n_screen_healthy <- length(idx_screen)
  
  # Initialize test counts as 0
  if (n_screen_healthy > 0) {
    m_patients[, `:=` (temp_ct_tests_screen = 0,
                       temp_ct_tests_screen_FP = 0)]
  }
  
  # Run while loop
  while (n_screen_healthy > 0) {
    if (verbose) print(n_screen_healthy)
    if (p_spec < 1) { # Imperfect specificity
      # Generate number of screening tests before first false positive
      m_patients[idx_screen, `:=` (
        ct_tests_screen_before_FP = rgeom(
          .N,
          prob = 1 - p_spec
        ))]
      
      # Calculate number of additional screening tests in healthy state
      m_patients[idx_screen, `:=` (ct_tests_screen_add = pmin(ct_tests_screen_before_FP + 1, 
                                                              (time_screen_stop_H - screen_age) %/% int_screen + ((time_screen_stop_H - screen_age) %% int_screen > 0), na.rm = T))] 
      
      # Increment number of screening tests and false positive tests for individuals with screening age before end of healthy state
      m_patients[idx_screen, `:=` (
        temp_ct_tests_screen = temp_ct_tests_screen + ct_tests_screen_add,
        temp_ct_tests_screen_FP = temp_ct_tests_screen_FP + fifelse(ct_tests_screen_add <= ct_tests_screen_before_FP, 0, 1),
        screen_age = screen_age + ct_tests_screen_add*int_screen
      )]
    } else { # Perfect specificity
      # Calculate number of additional screening tests in healthy state
      m_patients[idx_screen, `:=` (ct_tests_screen_add = (time_screen_stop_H - screen_age) %/% int_screen + ((time_screen_stop_H - screen_age) %% int_screen > 0))] 
      
      # Increment number of screening tests and false positive tests for individuals with screening age before end of healthy state
      m_patients[idx_screen, `:=` (
        temp_ct_tests_screen = temp_ct_tests_screen + ct_tests_screen_add,
        screen_age = screen_age + ct_tests_screen_add*int_screen
      )]
    }
    
    # Recalculate indices and number of healthy patients remaining to screen
    idx_screen <- idx_screen[m_patients[idx_screen, which(screen_age < time_screen_stop_H)]]
    n_screen_healthy <- length(idx_screen)
  }
  
  # Set variables that are no longer needed to NULL
  if (p_spec < 1) { # Imperfect specificity
    m_patients[, `:=` (ct_tests_screen_before_FP = NULL,
                       ct_tests_screen_add = NULL)]
  } else {
    m_patients[, `:=` (ct_tests_screen_add = NULL)]
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
  # Extract test sensitivity, specificity, and interval
  p_sens <- l_test_chars[[l_params_strategy[["mod"]]]][["p_sens"]][["L"]]
  p_spec <- l_test_chars[[l_params_strategy[["mod"]]]][["p_spec"]]
  int_screen <- l_params_strategy[["int_screen"]]
  
  # Set patient and lesion IDs as keys
  setkey(m_lesions, pt_id, lesion_id)
  
  # Set initial stop time for lesion screening (minimum of cancer onset and screening stop date, to be updated after every screening where lesions are removed)
  m_patients[!is.na(screen_age), time_screen_stop_L := pmin(time_H_P, time_screen_stop, na.rm = T)]
  
  # Extract row indices of individuals that will be screened in the lesion state
  idx_screen_pt <- m_patients[screen_age >= time_H_L & screen_age < time_screen_stop_L, which = TRUE]
  
  # Merge screening variables to lesion data and initialize flag for whether lesion was removed
  m_lesions[m_patients[idx_screen_pt], `:=` (
    screen_age = i.screen_age, # Current age at screening
    time_screen_stop_L = i.time_screen_stop_L, # Time that screening in the lesion state ends (updated if lesions are removed)
    time_screen_stop_max = pmin(time_H_Do, l_params_strategy$age_screen_stop + 1), # Maximum screening stop time due to death from other causes or screening stop age
    fl_removed = 0)]
  
  # Set patient-level screening age to NA for individuals that will be screened in the lesion state, to be updated
  m_patients[idx_screen_pt, screen_age := NA]
  
  # Initialize running variable of row indices of lesions that may be detected with screening 
  # (note: will not include lesions that develop after maximum screening stop time or that are removed)
  idx_screen_lesion <- m_lesions[!is.na(screen_age) & time_H_Lj < time_screen_stop_max, which = TRUE]
  
  # Get number of lesions available for screening
  n_screen_lesion <- length(idx_screen_lesion)
  
  # Loop through screening tests until there are no more lesions to screen
  if (verbose) print("Number of lesions remaining to screen:")
  while (n_screen_lesion > 0) {
    if (verbose) print(n_screen_lesion)
    
    ###### 2.0 Flag eligible lesions
    # Get indices of lesions that are present at current screen time
    idx_present <- idx_screen_lesion[m_lesions[idx_screen_lesion, which(time_H_Lj <= screen_age)]]
    
    if (length(idx_present) > 0) {
      ###### 2.1 Flag whether eligible lesions would be detected
      #### Sample whether eligible lesions get detected and thus removed
      m_lesions[idx_present, `:=` (
        fl_removed = rbinom(
          .N,
          size = 1,
          prob = p_sens
        ))]
      
      ###### 2.2 Process outcomes from removed lesions
      # Get row indices of newly removed lesions
      idx_removed <- idx_present[m_lesions[idx_present, which(fl_removed == 1)]]
      
      # For removed lesions, set detection time
      m_lesions[idx_removed, `:=` (time_detected = screen_age)]
      
      # Update time to preclinical cancer for patients with removed lesions
      m_removed_preclin <- m_lesions[pt_id %in% m_lesions[idx_removed, pt_id] & fl_removed == 0, .(time_H_P = min(time_H_Pj)), by = pt_id]
      
      # Merge preclinical cancer data to patient IDs of all present lesions (so IDs with 0 remaining lesions have NA preclinical onset time)
      m_update <- merge(m_lesions[idx_removed, .(pt_id = unique(pt_id))], m_removed_preclin, by = "pt_id", all = TRUE)
      setkey(m_update, pt_id)
      
      # Merge updated preclinical cancer times to lesion data and update lesion screening end time
      m_lesions[m_update, `:=` (time_H_P = i.time_H_P,
                                time_screen_stop_L = pmin(i.time_H_P, time_screen_stop_max, na.rm = T))]
      
      ###### 2.3 Update screening age and test count for patients with lesions present
      # Screening age
      m_lesions[pt_id %in% m_lesions[idx_present, pt_id], `:=` (screen_age = screen_age + int_screen)]
      
      # Test count
      m_patients[pt_id %in% m_lesions[idx_present, pt_id], ct_tests_screen := ct_tests_screen + 1]
      
      # Positive test count
      m_patients[pt_id %in% m_lesions[idx_removed, pt_id], ct_tests_positive := ct_tests_positive + 1]
    }
    
    
    ###### 2.4 Flag false positives among screened individuals with no active lesions
    # Get patient IDs with no lesions present at this round of screening
    # (this will happen if an earlier lesion was found and removed, but patient will develop another lesion in the future)
    v_pt_none <- setdiff(m_lesions[idx_screen_lesion, unique(pt_id)], m_lesions[idx_present, unique(pt_id)])
    
    if (length(v_pt_none > 0)) {
      # Extract screening info for patient
      m_none <- m_lesions[pt_id %in% v_pt_none][rowid(pt_id) == 1]
      setkey(m_none, pt_id)
      
      # Find earliest of next lesion and screening stop time
      m_none_censor <- m_lesions[pt_id %in% v_pt_none & fl_removed == 0, .(min_time_H_Lj = min(time_H_Lj)), by = pt_id]
      setkey(m_none_censor, pt_id)
      
      # Merge data tables
      m_none <- merge(m_none, m_none_censor, all = TRUE)
      
      # Merge screening info to patient data
      m_patients[m_none, `:=` (
        screen_age = i.screen_age,
        time_screen_stop_H = pmin(i.time_screen_stop_L, i.min_time_H_Lj, na.rm = T) # Censor time for healthy within lesion state
      )]
      
      # Get indices of patients for screening loop
      idx_none <- m_patients[pt_id %in% v_pt_none, which = TRUE]
      
      # Loop over healthy time in lesion state
      run_screening_H_loop(m_patients,
                           p_spec,
                           int_screen,
                           idx_none,
                           verbose)
      
      # Update screening and false positive test counts
      m_patients[idx_none, `:=` (ct_tests_screen = ct_tests_screen + temp_ct_tests_screen,
                                 ct_tests_positive = ct_tests_positive + temp_ct_tests_screen_FP)]
      
      # Set variables that are no longer needed to NULL
      m_patients[, `:=` (time_screen_stop_H = NULL,
                         temp_ct_tests_screen = NULL,
                         temp_ct_tests_screen_FP = NULL)]
      
      # Merge updated screening age to lesion-level data
      m_lesions[m_patients[idx_none], `:=` (screen_age = i.screen_age)]
    }
    
    ###### 2.5 Recalculate number of lesions remaining to screen
    # Update list of screen-eligible lesions to exclude removed lesions and patients that are no longer in lesion state
    idx_screen_lesion <- setdiff(idx_screen_lesion, idx_removed)
    idx_screen_lesion <- idx_screen_lesion[m_lesions[idx_screen_lesion, which(screen_age < time_screen_stop_L)]]
    
    # Recalculate number of lesions remaining to screen
    n_screen_lesion <- length(idx_screen_lesion)
  }
  
  ###### 2.6 Flag false positives among individuals whose screen-eligible lesions have all been removed
  # Find individuals who are still screened in lesion state but have not reached preclinical stage
  m_none_final <- m_lesions[screen_age < time_screen_stop_L][rowid(pt_id) == 1]
  
  if (nrow(m_none_final) > 0) {
    # Set screening variables
    setkey(m_none_final, pt_id)
    setnames(m_none_final, old = c("time_screen_stop_L"), new = c("time_screen_stop_H")) # Censor time for healthy within lesion state
    
    # Loop over healthy time in lesion state
    run_screening_H_loop(m_none_final,
                         p_spec,
                         int_screen,
                         idx_screen = NULL,
                         verbose)
    
    # Update screening and false positive test counts
    m_patients[m_none_final, `:=` (ct_tests_screen = ct_tests_screen + i.temp_ct_tests_screen,
                                   ct_tests_positive = ct_tests_positive + i.temp_ct_tests_screen_FP)]
    
    # Merge updated screening age to lesion-level data
    m_lesions[m_none_final, `:=` (screen_age = i.screen_age)]
  }
  
  ###### 2.7 Update patient-level data for individuals with removed lesions
  # Merge updated preclinical cancer times to patient data
  m_preclin_final <- m_lesions[fl_removed == 1][rowid(pt_id) == 1, .(pt_id, time_H_P)]
  setkey(m_preclin_final, pt_id)
  m_patients[m_preclin_final, `:=` (time_H_P = i.time_H_P)]
  
  # Update mortality outcomes
  idx_preclin <- m_patients[pt_id %in% m_preclin_final$pt_id, which = TRUE]
  m_patients[idx_preclin, `:=` (time_H_C = time_H_P + time_P_C)]
  m_patients[idx_preclin, `:=` (time_H_Dc = time_H_C + time_C_Dc)]
  calc_mortality_outcomes(m_patients, idx_subset = idx_preclin)
  
  # Update censoring variables
  m_patients[idx_preclin, `:=` (time_screen_censor = pmin(time_H_C, time_H_D, na.rm = T))]
  m_patients[idx_preclin, `:=` (time_screen_stop = pmin(time_screen_censor, l_params_strategy$age_screen_stop + 1))]
  
  ###### 2.8 Update screening details for individuals still undergoing screening in preclinical cancer state
  # Merge screening age
  m_screen_final <- m_lesions[!is.na(screen_age)][rowid(pt_id) == 1]
  m_patients[m_screen_final, `:=` (screen_age = i.screen_age)]
  
  # Reset screening age to NA if after censor date or end age
  m_patients[screen_age >= time_screen_stop, screen_age := NA]
  
  ###### 2.9 Remove unneeded variables
  m_lesions[, `:=` (screen_age = NULL,
                    time_screen_stop_L = NULL,
                    time_screen_stop_max = NULL)]
  
  return(NULL)
}


# Simulate screening occurring during the preclinical cancer state
# Assumes that screening during healthy state has been completed for population,
# which initializes the screen_age variable
simulate_screening_P <- function(m_patients,
                                 l_params_model,
                                 l_params_strategy,
                                 l_test_chars,
                                 nondecreasing = TRUE, # Should screening strictly increase life expectancy?
                                 verbose = FALSE
) {  
  # Extract test sensitivity and interval
  p_sens <- l_test_chars[[l_params_strategy[["mod"]]]][["p_sens"]][["P"]]
  int_screen <- l_params_strategy[["int_screen"]]
  
  # Initialize running count of patients remaining to be screened during preclinical cancer state
  idx_screen <- m_patients[!is.na(screen_age), which = TRUE]
  
  # Run while loop until no preclinical cancers remain to screen
  if (verbose) print("Number of individuals with preclinical cancer remaining to screen:")
  while (length(idx_screen) > 0) {
    if (verbose) print(length(idx_screen))
    # Increment test count and sample whether cancer would be detected
    m_patients[idx_screen, `:=` (
      ct_tests_screen = ct_tests_screen + 1,
      fl_detected = rbinom(
        .N,
        size = 1,
        prob = p_sens
      ))]
    
    # Reset clinical cancer age, set next screen age to NA, and increment positive count for newly detected cases
    idx_detected <- idx_screen[m_patients[idx_screen, which(fl_detected == 1)]]
    m_patients[idx_detected, `:=` (
      ct_tests_positive = ct_tests_positive + 1,
      time_H_C = screen_age,
      time_P_C = screen_age - time_H_P,
      screen_age = NA)]
    
    # Remove detected cases from screening-eligible patients
    idx_screen <- setdiff(idx_screen, idx_detected)
    
    # Update screen age, reset to NA if after stop date, update list of screen-eligible patients
    m_patients[idx_screen, screen_age := screen_age + int_screen]
    idx_censor <- idx_screen[m_patients[idx_screen, which(screen_age >= time_screen_stop)]]
    m_patients[idx_censor, screen_age := NA]
    idx_screen <- setdiff(idx_screen, idx_censor)
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
  
  # Calculate death from cancer
  # If applicable, restrict time to death from cancer to be no less than time to death without screening
  if (nondecreasing) {
    m_patients[fl_detected == 1, time_H_Dc_screen := time_H_C + time_C_Dc]
    m_patients[fl_detected == 1, time_H_Dc := pmax(time_H_Dc, time_H_Dc_screen)]
  } else {
    m_patients[fl_detected == 1, time_H_Dc := time_H_C + time_C_Dc]
  }
  
  # Flag cancers diagnosed with symptoms and increment base confirmatory test count
  if ("ct_tests_base" %in% names(m_patients)) m_patients[, ct_tests_base := NULL]
  m_patients[!fl_detected %in% 1 & time_H_C < time_H_Do, `:=` (
    fl_symptom_detected = 1,
    ct_tests_base = 1)]
  
  # Remove unneeded variables
  m_patients[, `:=` (time_P_C_running = NULL,
                     screen_age = NULL,
                     time_screen_stop = NULL)]
  
  return(NULL)
}


################################################################################
# Summary functions
################################################################################

# Calculate life years gained (LYG) from screening
calc_lyg <- function(res_base, res_screen, unit = 1) {
  lyg <- (res_screen["time_total"] - res_base["time_total"]) / res_base["N"] * unit
  return(lyg)
}


# Calculate confirmatory tests in base case scenario
calc_ntests <- function(m_patients,
                        censor_var = "time_screen_censor",
                        min_age = 0,
                        test_var_pattern = "ct_tests_"
) {
  # Calculate confirmatory tests
  if (min_age > 0 | !is.null(min_age)) {
    ct_tests <- m_patients[get(censor_var) > min_age, colSums(.SD, na.rm = T), .SDcols = patterns(test_var_pattern)]
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


################################################################################
# Wrapper functions for generating and summarizing data
################################################################################

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
                              conf_level = 0.95,
                              unit_test = FALSE
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
  
  # If unit test, make copy of initial data
  if (unit_test) {
    m_cohort_orig <- copy(m_cohort)
  }
  
  # Get base outputs
  l_outputs <- calc_cohort_outputs(m_cohort, 
                                   l_params_outcome = l_params_outcome,
                                   l_censor_vars = l_censor_vars)
  
  # If unit test, check that data has not changed
  if (unit_test) {
    if (is.data.table(m_cohort_orig)) {
      assert_that(all.equal(m_cohort_orig, m_cohort[, -c("time_screen_censor")]))
    } else {
      for (dt in names(m_cohort_orig)) {
        assert_that(all.equal(m_cohort_orig[[dt]], m_cohort[[dt]][, -c("time_screen_censor")]))
      }
    }
  }
  
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
      # Generate data under screening counterfactual (save results from base cohort that get overwritten)
      l_data_overwritten <- do.call(run_screening_counterfactual, 
                                    list(m_cohort = m_cohort, 
                                         l_params_model = l_params_model,
                                         l_params_strategy = l_params_screen$strats[[strat]],
                                         l_params_tests = l_params_screen$test_chars
                                    ))
      
      # Calculate screening outcomes from updated data
      l_outputs_screen[[strat]] <- calc_cohort_outputs(m_cohort, 
                                                       l_params_outcome = l_params_outcome_screen,
                                                       l_censor_vars = l_censor_vars)
      
      # Revert overwritten data
      m_cohort_screen <- swap_data(data_start = m_cohort, 
                                   data_swap = l_data_overwritten$m_cohort,
                                   cols_keep = l_data_overwritten$cols_orig)
      
      # Add individual-level data changed in screening scenario to results list if necessary
      if (individual_data) {
        res$m_cohort_screen[[strat]] <- m_cohort_screen
      }
      
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
    
    # If unit test, check that data has not changed
    if (unit_test) {
      if (is.data.table(m_cohort_orig)) {
        assert_that(all.equal(m_cohort_orig, m_cohort[, -c("time_screen_censor")]))
      } else {
        for (dt in names(m_cohort_orig)) {
          assert_that(all.equal(m_cohort_orig[[dt]], m_cohort[[dt]][, -c("time_screen_censor")]))
        }
      }
    }
  }
  
  # Return single item or list of results if >1 items
  if (length(res) == 1) {
    return(res[[1]])
  } else {
    return(res)
  }
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


# Swap data from specified rows and columns in list of data tables
# data_start: data table or list of data tables to be updated
# data_swap: data table or list of data tables with rows and columns to swap in
# cols_keep: vectors of strings with column names to keep in final result
# Modifies data_start in place and returns data that changed, including columns removed because they were not in cols_keep 
swap_data <- function(data_start, data_swap, cols_keep = NULL) {
  if (is.data.table(data_start)) { # For single data table
    data_replaced <- swap_data_table(data_start, data_swap, cols_keep)
  } else {
    data_replaced <- list()
    for (i in names(data_start)) {
      data_replaced[[i]] <- swap_data_table(data_start[[i]], data_swap[[i]], cols_keep[[i]])
    }
  }
  return(data_replaced)
}

# Internal function to swap data for single data table
swap_data_table <- function(data_start, data_swap, cols_keep = NULL) {
  browser()
  # Error handling
  if (!is.data.table(data_swap) | !is.data.table(data_start)) {
    stop("All items must be data tables")
  }
  
  # Set key if necessary
  if (!haskey(data_swap)) {
    setkey(data_swap, key(data_start))
  }
  
  # Extract columns to swap
  v_cols <- names(data_swap)
  v_cols_diff <- v_cols[!v_cols %in% key(data_start)]
  
  # Extract rows that will be switched
  data_replaced <- data_start[data_swap, nomatch = 0, .SD, .SDcols = v_cols]
  
  # Swap data 
  data_start[data_swap, (v_cols_diff) := mget(paste0("i.", v_cols_diff))]
  
  # Extract specified columns to data_replaced if applicable
  if (!is.null(cols_keep)) {
    # Select ID variable and columns to remove
    cols_remove <- setdiff(names(data_start), cols_keep)
    
    # Merge data_replace with removed columns
    if (length(cols_remove) > 0) {
      cols_merge <- c(key(data_start), cols_remove)
      data_removed <- data_start[, ..cols_merge]
      data_replaced <- merge(data_replaced, data_removed, all = TRUE)
      
      # Remove unwanted columns from data_start
      data_start[, (cols_remove) := NULL]
    }
  }
  
  return(data_replaced)
}