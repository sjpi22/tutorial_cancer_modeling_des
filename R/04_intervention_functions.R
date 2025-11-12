
# Rerun version for screening (note: overwrites input data by reference)
run_screening_counterfactual <- function(
    m_cohort,
    l_params_model,
    l_params_strategy,
    l_params_tests,
    l_params_surveil = NULL,
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
  
  # Extract screening and confirmatory test modalities
  v_mods <- l_params_strategy[["mod"]]
  if (!is.null(l_params_strategy[["mod_conf"]])) {
    v_mods <- c(v_mods, l_params_strategy[["mod_conf"]])
  }
  
  # Add surveillance modalities
  if (!is.null(l_params_surveil)) {
    v_mods <- c(v_mods, l_params_surveil[["mod"]])
    if (!is.null(l_params_surveil[["mod_conf"]])) {
      v_mods <- c(v_mods, l_params_surveil[["mod_conf"]])
    }
    v_mods <- unique(v_mods)
  }
  
  # Sample sensitivity and specificity of screening and (if applicable) confirmatory tests
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
  }
  
  #### Case 2: Lesions developed by screen age, but no preclinical cancer ####
  if ('L' %in% l_params_model$v_states) {
    # Simulate lesion progression to cancer
    simulate_screening_L(m_patients = m_patients, 
                         m_lesions = m_lesions,
                         l_params_strategy = l_params_strategy,
                         l_test_chars = l_params_tests_sample,
                         l_params_surveil = l_params_surveil,
                         verbose = verbose)
  }
  
  #### Case 3: Preclinical cancer screening and surveillance ####
  simulate_screening_P(m_patients = m_patients,
                       l_params_model = l_params_model,
                       l_params_strategy = l_params_strategy,
                       l_test_chars = l_params_tests_sample,
                       l_params_surveil = l_params_surveil,
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
      
      # Assign test type
      res$type <- type
      return(res)
    })
  }
  return(l_params_tests_sample)
}


# Set surveillance interval for negative case (no lesions detected)
set_int_surveil_neg <- function(l_params_surveil, l_params_strategy, l_test_types) {
  int_surveil_neg <- c()
  if (0 %in% l_params_surveil[["n_detected"]]) {
    # If case for 0 lesions is populated in n_detected and int_surveil, use as next surveillance screen interval
    int_surveil_neg["screen"] <- l_params_surveil[["int_surveil"]][[which(l_params_surveil[["n_detected"]] == 0)]]
  } else {
    # If case for 0 lesions is not populated in n_detected and int_surveil, use defaults from routine screening
    if (l_params_surveil[["mod"]] %in% l_test_types[["direct"]]) {
      # For direct test, use the direct test interval from routine screening
      if (l_params_strategy[["mod"]] %in% l_test_types[["direct"]]) {
        int_surveil_neg["screen"] <- l_params_strategy[["int_screen"]]
      } else {
        int_surveil_neg["screen"] <- l_params_strategy[["int_conf"]]
      }
    } else {
      # For indirect test, use the routine screening interval
      int_surveil_neg["screen"] <- l_params_strategy[["int_screen"]]
    }
  }
  
  # Populate surveillance confirmatory test interval if primary surveillance test is not confirmatory
  if (!l_params_surveil[["mod"]] %in% l_test_types[["direct"]]) {
    # Use confirmatory interval reported surveillance parameters if possible
    if ("int_conf" %in% names(l_params_surveil)) {
      int_surveil_neg["conf"] <- l_params_surveil[["int_conf"]]
    } else {
      # If confirmatory interval is not specified, use confirmatory test interval in screening
      if (v_mod["screen"] %in% l_test_types[["direct"]]) {
        # If primary modality is also confirmatory, use primary interval
        int_surveil_neg["conf"] <- l_params_strategy[["int_screen"]]
      } else {
        # If primary modality is not confirmatory, use confirmatory interval
        int_surveil_neg["conf"] <- l_params_strategy[["int_conf"]]
      }
    }
  }
  return(int_surveil_neg)
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
  mod_conf <- l_params_strategy[["mod_conf"]]
  
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
  
  # Set initial censor time for cancer screening
  if (!"time_screen_censor" %in% names(m_patients)) {
    m_patients[, time_screen_censor := pmin(time_H_C, time_H_D, na.rm = T)]
  }
  
  # Initialize type of screening (routine vs. surveillance), calculate age at 
  # which individual no longer screens as minimum of patient censor age and 
  # strategy stop age
  m_patients[, `:=` (
    screen_type = "screen",
    time_screen_stop = pmin(time_screen_censor, l_params_strategy$age_screen_stop + 1)
  )]
  
  # If screening interval is constant, perform simplified calculation of number of screening tests during healthy state
  # Otherwise, loop over testing rounds and assign applicable screening interval
  if (is.null(mod_conf) | int_conf %in% int_screen | p_spec == 1) {
    # Get routine screening ages
    v_screen_ages <- with(l_params_strategy, seq(age_screen_start, age_screen_stop, int_screen))
    
    # Get number of screening tests during healthy state before earliest of censor age or disease onset
    m_patients[, temp_ct_tests_screen := findInterval(pmin(time_screen_stop, get(var_onset), na.rm = T), v_screen_ages, left.open = T)]
    
    # Simulate number of false positives screens
    m_patients[, `:=` (
      ct_tests_screen_FP = rbinom(
        .N,
        size = temp_ct_tests_screen,
        prob = 1 - p_spec
      ))]
    
    # If screening test is followed by confirmatory test, set confirmatory tests equal to false positive screening tests
    if (!is.null(mod_conf)) {
      m_patients[, temp_ct_tests_conf := ct_tests_screen_FP]
    }
    
    # Get first screening age with disease present
    m_patients[, screen_age := v_screen_ages[temp_ct_tests_screen + 1]]
    
    # Reset first screening age to NA if after end date
    m_patients[screen_age >= time_screen_stop, screen_age := NA]
  } else { 
    # Initialize test strategy variables needed for screening loop
    m_patients[, `:=` (p_spec = p_spec,
                       int_screen = int_screen,
                       int_conf = int_conf)]
    
    # Initialize running screening age at screening start age for individuals eligible to receive screening at that age
    # and time that screening in the healthy state would end
    m_patients[l_params_strategy$age_screen_start < time_screen_stop, `:=` (screen_age = l_params_strategy$age_screen_start,
                                                                            time_screen_stop_H = pmin(get(var_onset), time_screen_stop, na.rm = T))]
    
    # Calculate number of healthy patients remaining to screen and save their row indices
    if (verbose) print("Number of healthy individuals remaining to be screened:")
    idx_screen <- m_patients[screen_age < time_screen_stop_H, which = TRUE]
    
    # Loop through screening tests until there are no more healthy individuals to screen
    run_screening_H_loop(m_patients,
                         idx_screen,
                         verbose)
    
    # Rename FP screening test variable
    setnames(m_patients, "temp_ct_tests_screen_FP", "ct_tests_screen_FP")
    
    # Update screen age to NA if after screening stop date
    m_patients[screen_age >= time_screen_stop, screen_age := NA]
    
    # Set number of FP as number of confirmatory tests, set variables that are no longer needed to NULL
    m_patients[, `:=` (temp_ct_tests_conf = ct_tests_screen_FP,
                       time_screen_stop_H = NULL,
                       p_spec = NULL,
                       int_screen = NULL,
                       int_conf = NULL)]
  }
  
  # Rename test variables
  setnames(m_patients, "temp_ct_tests_screen", paste0("ct_tests_", mod))
  if ("temp_ct_tests_conf" %in% names(m_patients)) {
    setnames(m_patients, "temp_ct_tests_conf", paste0("ct_tests_", mod_conf))
  }
  
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
                                 idx_screen, # Row indices of patients to screen
                                 verbose = FALSE) {
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
      screen_age = screen_age + fifelse(ct_tests_screen_add <= ct_tests_screen_before_FP, 
                                        ct_tests_screen_add*int_screen, 
                                        (ct_tests_screen_add - 1)*int_screen + int_conf)
    )]
    
    # Recalculate indices and number of healthy patients remaining to screen
    idx_screen <- idx_screen[m_patients[idx_screen, which(screen_age < time_screen_stop_H)]]
    n_screen_healthy <- length(idx_screen)
  }
  
  # Set variables that are no longer needed to NULL
  m_patients[, `:=` (ct_tests_screen_before_FP = NULL,
                     ct_tests_screen_add = NULL)]
                     
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
                                 l_params_surveil = NULL,
                                 verbose = FALSE
) {  
  # Extract screening and surveillance primary test modalities
  v_mod <- c(screen = l_params_strategy[["mod"]])
  
  # Extract screening and surveillance confirmatory test modalities
  v_mod_conf <- c(screen = ifelse(is.null(l_params_strategy[["mod_conf"]]), NA, l_params_strategy[["mod_conf"]]))
  
  # Extract modalities associated with each test type (direct, targeted, indirect)
  df_test_types <- data.frame(mod = names(l_test_chars),
                              type = sapply(l_test_chars, function(u) {
                                u[["type"]]
                              }))
  l_test_types <- split(df_test_types$mod, df_test_types$type)
  
  # Process surveillance parameters
  if (!is.null(l_params_surveil)) {
    # Extract surveillance primary and confirmatory test modalities
    v_mod["surveil"] <- l_params_surveil[["mod"]]
    v_mod_conf["surveil"] <- ifelse(is.null(l_params_surveil[["mod_conf"]]), NA, l_params_surveil[["mod_conf"]])
    
    # Create expression for surveillance interval based on number of lesions detected
    expr_surveil <- "fcase("
    for (n in 1:length(l_params_surveil[["n_detected"]])) {
      expr_surveil <- paste0(expr_surveil, "ct_removed >= ", l_params_surveil[["n_detected"]][n], ", ", l_params_surveil[["int_surveil"]][n])
      if (n < length(l_params_surveil[["n_detected"]])) {
        expr_surveil <- paste0(expr_surveil, ", ")
      }
    }
    expr_surveil <- paste0(expr_surveil, ")")
    
    # Set surveillance interval for no lesions based on routine screening intervals
    int_surveil_neg <- set_int_surveil_neg(l_params_surveil, l_params_strategy, l_test_types)
  }
  
  # Ensure that all modalities are represented in patient data counts
  for (mod in unique(c(v_mod, v_mod_conf[!is.na(v_mod_conf)]))) {
    if (!paste0("ct_tests_", mod) %in% names(m_patients)) {
      m_patients[, paste0("ct_tests_", mod) := 0]
    }
  }
  
  # Set patient and lesion IDs as keys
  setkey(m_lesions, pt_id, lesion_id)
  
  # Set initial stop time for lesion screening (minimum of cancer onset and screening stop date, to be updated after every screening where lesions are removed)
  m_patients[!is.na(screen_age), time_screen_stop_L := pmin(time_H_P, time_screen_stop, na.rm = T)]
  
  # Extract row indices of individuals that will be screened in the lesion state
  idx_screen_pt <- m_patients[screen_age >= time_H_L & screen_age < time_screen_stop_L, which = TRUE]
  
  # Merge screening variables to lesion data and initialize flag for whether lesion was removed
  m_lesions[m_patients[idx_screen_pt], `:=` (
    screen_age = i.screen_age, # Current age at screening
    screen_type = i.screen_type, # Type of regimen (screening vs. surveillance)
    modality = v_mod[screen_type], # Screening test modality
    int_test_next = l_params_strategy[["int_screen"]], # Default interval to next test (updated if confirmatory test or surveillance)
    time_screen_stop_L = i.time_screen_stop_L, # Time that screening in the lesion state ends (updated if lesions are removed)
    time_screen_stop_max = pmin(time_H_Do, l_params_strategy$age_screen_stop + 1), # Maximum screening stop time due to death from other causes or screening stop age
    fl_removed = 0)]
  
  # Add confirmatory test modality if applicable
  if (any(!is.na(v_mod_conf[unique(m_lesions$screen_type)]))) {
    m_lesions[, modality_conf := v_mod_conf[screen_type]]
  }
  
  # Initialize running variable of row indices of lesions that may be detected with screening 
  # (note: will not include lesions that develop after maximum screening stop time or that are removed)
  idx_screen_lesion <- m_lesions[!is.na(screen_age) & time_H_Lj < time_screen_stop_max, which = TRUE]
  
  # Get row indices of lesions that emerge after maximum screening stop time
  idx_postscreen <- m_lesions[!is.na(screen_age) & time_H_Lj >= time_screen_stop_max, which = TRUE]
  
  # Get number of lesions available for screening
  n_screen_lesion <- length(idx_screen_lesion)
  
  # Loop through screening tests until there are no more lesions to screen
  if (verbose) print("Number of lesions remaining to screen:")
  while (n_screen_lesion > 0) {
    if (verbose) print(n_screen_lesion)
    
    ###### 2.0 Flag eligible lesions
    # Get indices of lesions that are present at current screen time
    idx_present <- idx_screen_lesion[m_lesions[idx_screen_lesion, which(time_H_Lj <= screen_age)]]
    
    
    ###### 2.1 Flag whether eligible lesions would be detected and increment test counts
    if (length(idx_present) > 0) {
      #### Increment screening modality counts
      for (mod in unique(m_lesions[idx_present, modality])) {
        colname <- paste0("ct_tests_", mod)
        pt_subset <- m_lesions[idx_present][modality == mod, unique(pt_id)]
        m_patients[pt_id %in% pt_subset, (colname) := get(colname) + 1]
      }
      
      # Initialize vector to save patient IDs with positive noninvasive test
      pt_noninvasive_positive <- c()
      
      #### Sample whether eligible lesions produce positive test result
      m_lesions[idx_present, `:=` (
        fl_positive = rbinom(
          .N,
          size = 1,
          prob = l_test_chars[[modality]][["p_sens"]][["L"]]
        )), by = modality]
      
      # Get indices of lesions that test positive
      idx_positive <- idx_present[m_lesions[idx_present, which(fl_positive == 1)]]
      
      #### Apply downstream effect of positive result depending on test type (direct, targeted, indirect) ####
      ### If test is direct, positive lesions are detected and removed
      if ("direct" %in% names(l_test_types)) {
        # Get indices of positive lesions that underwent direct visualization test with biopsy
        idx_subset <- idx_positive[m_lesions[idx_positive, which(modality %in% l_test_types[["direct"]])]]
        
        # Update detected lesions as removed
        if (length(idx_subset > 0)) {
          m_lesions[idx_subset, `:=` (fl_removed = 1)]
        }
      }
      
      ### If test is targeted, apply confirmatory tests only to positive lesions to sample whether they are removed
      if ("targeted" %in% names(l_test_types)) {
        # Get indices of positive lesions that underwent targeted test
        idx_subset <- idx_positive[m_lesions[idx_positive, which(modality %in% l_test_types[["targeted"]])]]
        
        if (length(idx_subset > 0)) {
          # Save patient IDs with positive targeted test
          pt_noninvasive_positive <- c(pt_noninvasive_positive, m_lesions[idx_subset][modality == mod, unique(pt_id)])
          
          # Flag that confirmatory test occurred and whether positive lesions are detected and removed
          m_lesions[idx_subset, `:=` (
            fl_removed = rbinom(
              .N,
              size = 1,
              prob = l_test_chars[[modality_conf]][["p_sens"]][["L"]]
            )), by = modality_conf]
        }
      }
      
      ### If test is indirect (non-targeted), apply confirmatory tests to all lesions if at least one produces positive result
      if ("indirect" %in% names(l_test_types)) {
        # Get patient IDs of positive lesions that underwent indirect test
        pt_positive <- m_lesions[idx_positive][modality %in% l_test_types[["indirect"]], unique(pt_id)]
        
        if (length(pt_positive > 0)) {
          # Add patient IDs to vector of noninvasive positives
          pt_noninvasive_positive <- c(pt_noninvasive_positive, pt_positive)
          
          # Get indices of lesions that are present among patients with follow-on test after indirect test
          idx_subset <- idx_present[m_lesions[idx_present, which(pt_id %in% pt_positive)]]
          
          # Sample removal with confirmatory test among above lesions
          m_lesions[idx_subset,`:=` (
            fl_removed = rbinom(
              .N,
              size = 1,
              prob = l_test_chars[[modality_conf]][["p_sens"]][["L"]]
            )), by = modality_conf]
        }
      }
      
      #### Update time interval and confirmatory test count for positive non-invasive tests (implying there was confirmatory test)
      if (length(pt_noninvasive_positive) > 0) {
        # Set time after confirmatory test depending on regimen type (screen or surveillance)
        m_lesions[pt_id %in% pt_noninvasive_positive, `:=` (int_test_next = l_params_strategy[[screen_type]][["int_conf"]])]
        
        # Increment confirmatory modality counts
        for (mod in unique(m_lesions[pt_id %in% pt_noninvasive_positive, modality_conf])) {
          colname <- paste0("ct_tests_", mod)
          pt_subset <- m_lesions[pt_id %in% pt_noninvasive_positive][modality_conf == mod, unique(pt_id)]
          m_patients[pt_id %in% pt_subset, (colname) := get(colname) + 1]
        }
      }
      
      
      # Update time to preclinical cancer and following times, increment number of screens, number of confirmatory tests, and screen age
      m_patients[m_lesions_summary, `:=` (time_H_P = i.time_H_P,
                                          screen_type = i.screen_type,
                                          screen_age = screen_age + i.int_test_next,
                                          modality = i.modality,
                                          modality_conf = i.modality_conf)]
      
      # Update time for people with no removed lesions
      
      ###### 2.2 Process outcomes from removed lesions
      # Get row indices of newly removed lesions
      idx_removed <- idx_positive[m_lesions[idx_positive, which(fl_removed == 1)]]
      
      # For removed lesions, set detection time
      m_lesions[idx_removed, `:=` (time_detected = screen_age)]
      
      # Update time to preclinical cancer for patients with removed lesions
      m_removed_preclin <- m_lesions[pt_id %in% m_lesions[idx_removed, pt_id] & fl_removed == 0, .(time_H_P := min(time_H_Pj)), by = pt_id]
      
      # Merge preclinical cancer data to all IDs (so IDs with 0 remaining lesions have NA preclinical onset time)
      m_update <- merge(m_lesions[idx_removed, .(pt_id)], m_removed_preclin, by = "pt_id", all = TRUE)
      setkey(m_update, pt_id)
      
      # Merge updated preclinical cancer times to lesion data and update lesion screening end time
      m_lesions[m_update, `:=` (time_H_P = i.time_H_P,
                                time_screen_stop_L = pmin(i.time_H_P, time_screen_stop_max, na.rm = T))]
      
      # If there is surveillance, apply surveillance regimen
      if (!is.null(l_params_surveil)) {
        # Summarize number of removed lesions per patient
        m_removed_ct <- m_lesions[idx_removed, .(ct_removed = .N), by = pt_id]
        
        # Set next surveillance interval based on number of lesions detected
        m_removed_ct[, `:=` (
          int_test_next = eval(parse(text = expr_surveil))
        )]
        
        # Merge to lesion data
        m_lesions[m_removed_ct, `:=` (screen_type = "surveil",
                                      modality = v_mod["surveil"], # Screening test modality
                                      int_test_next = i.int_test_next)]
        
        # Add confirmatory test modality if applicable
        if (any(!is.na(v_mod_conf[unique(m_lesions[m_removed_ct, screen_type])]))) {
          m_lesions[m_removed_ct, modality_conf := v_mod_conf[screen_type]]
        }
      }
      
      ###### 2.3 Update screening age and default test interval for patients with lesions present
      # Screening age
      m_lesion[pt_id %in% m_lesion[idx_present], `:=` (screen_age = screen_age + int_test_next)]
      
      # TO DO: for positive noninvasive and surveillance cases, change interval back to default of none detected
      
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
      
      # Merge screening info to patient data and set test modality
      m_patients[m_none, `:=` (
        screen_age = i.screen_age,
        screen_type = i.screen_type,
        modality = v_mod[i.screen_type],
        modality_conf = v_mod_conf[i.screen_type],
        time_screen_stop_L = i.time_screen_stop_L, # Censor time for lesion state
        time_screen_stop_H = pmin(i.time_screen_stop_L, i.min_time_H_Lj, na.rm = T) # Censor time for healthy within lesion state
      )]
      
      # Get indices of patients for screening loop
      idx_none <- m_patients[pt_id %in% v_pt_none, which = TRUE]
      
      # Initialize test strategy variables needed for screening loop
      m_patients[idx_none, `:=` (p_spec = 1 - l_test_chars[[modality]][["p_spec"]])]
      
      # Set screening and confirmatory interval based on regimen type: screening
      idx_none_screen <- idx_none[m_patients[idx_none, which(screen_type == "screen")]]
      if (length(idx_none_subset) > 0) {
        m_patients[idx_none_subset, `:=` ( 
          int_screen = l_params_strategy[["int_screen"]],
          int_conf = l_params_strategy[["int_conf"]])]
      }
      
      # Set screening and confirmatory interval based on regimen type: surveillance
      idx_none_surveil <- idx_none[m_patients[idx_none, which(screen_type == "surveil")]]
      if (length(idx_none_surveil) > 0) {
        m_patients[idx_none_surveil, `:=` ( 
          int_screen = int_surveil_neg["screen"],
          int_conf = int_surveil_neg["conf"])]
      }
      
      # Loop over healthy time in lesion state
      run_screening_H_loop(m_patients,
                           idx_none,
                           verbose)
      
      # Update screening test counts
      for (mod in unique(m_patients[idx_none, modality])) {
        colname <- paste0("ct_tests_", mod)
        idx_subset <- idx_none[m_patients[idx_none, which(modality == mod)]]
        m_patients[idx_subset, (colname) := get(colname) + temp_ct_tests_screen]
      }
      
      # Update confirmatory test counts
      for (mod in unique(m_patients[idx_none, modality_conf])) {
        if (!is.na(mod)) {
          colname <- paste0("ct_tests_", mod)
          idx_subset <- idx_none[m_patients[idx_none, which(modality == mod)]]
          m_patients[idx_subset, (colname) := get(colname) + temp_ct_tests_screen_FP]
        }
      }
      
      # Update false positive test counts
      m_patients[idx_none, `:=` (ct_tests_screen_FP = ct_tests_screen_FP + temp_ct_tests_screen_FP)]
      
      # Set variables that are no longer needed to NULL
      m_patients[, `:=` (time_screen_stop_H = NULL,
                         temp_ct_tests_screen = NULL,
                         temp_ct_tests_screen_FP = NULL,
                         p_spec = NULL,
                         int_screen = NULL,
                         int_conf = NULL)]
      
      # For people who still receive screening in lesion state, merge updated screening age to lesion-level data
      idx_continue_lesion <- idx_none[m_patients[idx_none, which(screen_age < time_screen_stop_L)]]
      m_lesions[m_patients[idx_continue_lesion], `:=` (screen_age = i.screen_age)]
    }
    
    
    ###### 2.5 Recalculate number of lesions remaining to screen
    # Update list of screen-eligible lesions to exclude removed lesions and patients that are no longer in lesion state
    idx_screen_lesion <- setdiff(idx_screen_lesion, idx_removed)
    idx_screen_lesion <- idx_screen_lesion[m_lesions[idx_screen_lesion, which(screen_age < time_screen_stop_L)]]
    
    # Recalculate number of lesions remaining to screen
    n_screen_lesion <- length(idx_screen_lesion)
  }
  
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
  
  # Merge screening details for individuals still undergoing screening in preclinical cancer state
  m_screen_final <- m_lesions[!is.na(screen_age)][rowid(pt_id) == 1]
  m_patients[m_screen_final, `:=` (screen_age = i.screen_age,
                                   screen_type = i.screen_type)]
  
  # Reset screening age to NA if after censor date or end age
  m_patients[screen_age >= time_screen_stop, screen_age := NA]
  
  return(NULL)
}


# Simulate screening occurring during the preclinical cancer state
# Assumes that screening during healthy state has been completed for population,
# which initializes the screen_age variable
simulate_screening_P <- function(m_patients,
                                 l_params_model,
                                 l_params_strategy,
                                 l_test_chars,
                                 l_params_surveil = NULL,
                                 nondecreasing = TRUE, # Should screening strictly increase life expectancy?
                                 verbose = FALSE
) {  
  # Extract screening and surveillance primary test modalities
  v_mod <- c(screen = l_params_strategy[["mod"]])
  
  # Extract screening and surveillance confirmatory test modalities
  v_mod_conf <- c(screen = ifelse(is.null(l_params_strategy[["mod_conf"]]), NA, l_params_strategy[["mod_conf"]]))
  
  # Extract modalities associated with each test type (direct, targeted, indirect)
  df_test_types <- data.frame(mod = names(l_test_chars),
                              type = sapply(l_test_chars, function(u) {
                                u[["type"]]
                              }))
  l_test_types <- split(df_test_types$mod, df_test_types$type)
  
  # Process surveillance parameters
  if (!is.null(l_params_surveil)) {
    # Extract surveillance primary and confirmatory test modalities
    v_mod["surveil"] <- l_params_surveil[["mod"]]
    v_mod_conf["surveil"] <- ifelse(is.null(l_params_surveil[["mod_conf"]]), NA, l_params_surveil[["mod_conf"]])
    
    # Set surveillance interval for no lesions based on routine screening intervals
    int_surveil_neg <- set_int_surveil_neg(l_params_surveil, l_params_strategy, l_test_types)
  }
  
  # Initialize running count of patients remaining to be screened during preclinical cancer state
  idx_screen <- m_patients[!is.na(screen_age), which = TRUE]
  
  # Set screening and confirmatory modality based on routine vs. surveillance screening
  m_patients[idx_screen, `:=` (modality = v_mod[screen_type],
                               modality_conf = v_mod_conf[screen_type])]
  
  # Run while loop until no preclinical cancers remain to screen
  if (verbose) print("Number of individuals with preclinical cancer remaining to screen:")
  while (length(idx_screen) > 0) {
    if (verbose) print(length(idx_screen))
    # Sample whether cancer would produce positive test result
    m_patients[idx_screen, `:=` (
      fl_positive = rbinom(
        .N,
        size = 1,
        prob = l_test_chars[[modality]][["p_sens"]][["P"]]
      )), by = modality]
    
    # Get indices of positive results
    idx_positive <- idx_screen[m_patients[idx_screen, which(fl_positive == 1)]]
    
    #### Apply downstream effect of positive result depending on test type (direct, targeted, indirect) ####
    ### If test is direct, positive cancer is detected
    if ("direct" %in% names(l_test_types)) {
      idx_subset <- idx_positive[m_patients[idx_positive, which(modality %in% l_test_types[["direct"]])]]
      m_patients[idx_subset, `:=` (fl_conf = 0, fl_detected = 1)]
    }
    
    ### If test is targeted or indirect, apply confirmatory test to sample whether it is detected
    if ("targeted" %in% names(l_test_types) | "indirect" %in% names(l_test_types)) {
      idx_subset <- idx_positive[m_patients[idx_positive, which(modality %in% unlist(l_test_types[c("targeted", "indirect")]))]]
      if (length(idx_subset) > 0) {
        m_patients[idx_subset, `:=` (
          fl_conf = 1,
          fl_detected = rbinom(
            .N,
            size = 1,
            prob = l_test_chars[[modality_conf]][["p_sens"]][["P"]]
          )), by = modality_conf]
      }
    }
    
    # Add primary modality counts and delete flags
    for (mod in unique(v_mod)) {
      idx_subset <- idx_screen[m_patients[idx_screen, which(modality == mod)]]
      colname <- paste0("ct_tests_", mod)
      m_patients[idx_subset, (colname) := get(colname) + 1]
    }
    
    # Add confirmatory modality counts and delete flags
    if ("modality_conf" %in% names(m_patients)) {
      for (mod in unique(v_mod_conf[!is.na(v_mod_conf)])) {
        idx_subset <- idx_positive[m_patients[idx_positive, which(modality_conf == mod)]]
        colname <- paste0("ct_tests_", mod)
        m_patients[idx_subset, (colname) := get(colname) + 1]
      }
      m_patients[, `:=` (modality_conf = NULL)]
    }
    
    # Set next test time for negatives cases that were routinely screened
    m_patients[, int_test_next := NA_integer_]
    idx_undetected_screen <- idx_screen[m_patients[idx_screen, which(!fl_detected %in% 1 & screen_type == "screen")]]
    m_patients[idx_undetected_screen, `:=` (
      int_test_next = ifelse(!fl_conf %in% 1, 
                             l_params_strategy[["int_screen"]],
                             l_params_strategy[["int_conf"]])
    )]
    
    # Set next test time for negatives cases that were surveillance-screened
    idx_undetected_surveil <- idx_screen[m_patients[idx_screen, which(!fl_detected %in% 1 & screen_type == "surveil")]]
    m_patients[idx_undetected_surveil, `:=` (
      int_test_next = ifelse(!fl_conf %in% 1, 
                             int_surveil_neg["screen"],
                             int_surveil_neg["conf"])
    )]
    
    # Reset clinical cancer age and set next screen age to NA for newly detected cases
    idx_detected <- idx_screen[m_patients[idx_screen, which(fl_detected == 1)]]
    m_patients[idx_detected, `:=` (
      time_H_C = screen_age,
      time_P_C = screen_age - time_H_P,
      screen_age = NA)]
    
    # Update screen age, reset to NA if after stop date, update count of screen-eligible patients, and reset temporary variables
    m_patients[idx_screen, screen_age := screen_age + int_test_next]
    m_patients[screen_age >= time_screen_stop, screen_age := NA]
    idx_screen <- idx_screen[m_patients[idx_screen, which(!is.na(screen_age))]]
    m_patients[, `:=` (fl_positive = NULL,
                       int_test_next = NULL,
                       fl_conf = NULL)]
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
                     screen_type = NULL,
                     time_screen_stop = NULL)]
  
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
                                         l_params_tests = l_params_screen$test_chars,
                                         l_params_surveil = l_params_screen$surveil
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