
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
  l_params_tests_sample <- list()
  for (mod in v_mods) {
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
  
  #### Case 1: Simulate screening before disease ####
  # Set variable for time to disease onset
  var_onset <- paste("time", l_params_model$v_states[1], l_params_model$v_states[2], sep = "_")
  
  # Simulate screening during healthy state
  simulate_screening_H(m_patients = m_patients,
                       var_onset = var_onset,
                       l_params_strategy = l_params_strategy,
                       l_test_chars = l_params_tests_sample,
                       verbose = verbose)
  
  # Save data that will be replaced
  if (reversible) {
    v_cols_save <- c(
      "pt_id", "time_H_P", "time_H_C", "time_H_Dc", "time_H_D", "time_P_C", "time_C_Dc", "time_C_D", "stage_dx", "fl_Dc", "ct_tests_base"
    )
    
    if (!'L' %in% l_params_model$v_states) { # Time to cancer onset would not change if cancer is denovo
      v_cols_save <- v_cols_save[!v_cols_save %in% c("time_H_P")]
    }
    m_patient_overwritten <- m_patients[!is.na(screen_age), .SD, .SDcols = intersect(v_cols_save, names(m_patients))]
  }
  
  #### Case 2: Lesions developed by screen age, but no preclinical cancer ####
  if ('L' %in% l_params_model$v_states) {
    # Simulate lesion progression to cancer
    m_lesion_overwritten <- simulate_screening_L(m_patients = m_patients, 
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
    time_screen_stop = pmin(time_screen_censor, l_params_strategy$age_screen_stop)
  )]
  
  # If screening interval is constant, perform simplified calculation of number of screening tests during healthy state
  # Otherwise, loop over testing rounds and assign applicable screening interval
  if (is.null(mod_conf) | int_conf %in% int_screen) {
    # Get routine screening ages
    v_screen_ages <- with(l_params_strategy, seq(age_screen_start, age_screen_stop, int_screen))
    
    # Get number of screening tests during healthy state before earliest of censor age or disease onset
    m_patients[, ct_tests_screen := findInterval(pmin(time_screen_stop, get(var_onset), na.rm = T), v_screen_ages, left.open = T)]
    
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
                       ct_tests_screen_FP = 0,
                       ct_tests_conf = 0)]
    
    # Initialize running screening age at screening start age for individuals eligible to receive screening at that age
    # and time that screening in the healthy state would end
    m_patients[l_params_strategy$age_screen_start < time_screen_stop, `:=` (screen_age = l_params_strategy$age_screen_start,
                                                                            time_screen_end_H = pmin(get(var_onset), time_screen_stop, na.rm = T))]
    
    # Calculate number of healthy patients remaining to screen
    if (verbose) print("Number of healthy individuals remaining to be screened:")
    n_healthy_screen <- m_patients[screen_age < time_screen_end_H, .N]
    
    # Run while loop
    while (n_healthy_screen > 0) {
      if (verbose) print(n_healthy_screen)
      # Generate number of screening tests before first false positive
      m_patients[time_screen_end_H > screen_age, `:=` (
        ct_tests_screen_before_FP = rgeom(
          .N,
          prob = 1 - p_spec
        ))]
      
      # Calculate number of additional screening tests in healthy state
      m_patients[time_screen_end_H > screen_age, `:=` (ct_tests_screen_add = pmin(ct_tests_screen_before_FP + 1, 
                                                                                  (time_screen_end_H - screen_age) %/% int_screen + ((time_screen_end_H - screen_age) %% int_screen > 0), na.rm = T))] 
      # Increment number of screening tests and false positive tests for individuals with screening age before end of healthy state
      m_patients[time_screen_end_H > screen_age, `:=` (
        ct_tests_screen = ct_tests_screen + ct_tests_screen_add,
        ct_tests_screen_FP = ct_tests_screen_FP + fifelse(ct_tests_screen_add <= pmin(Inf, ct_tests_screen_before_FP, na.rm = T), 0, 1),
        screen_age = screen_age + fifelse(ct_tests_screen_add <= pmin(Inf, ct_tests_screen_before_FP, na.rm = T), 
                                          ct_tests_screen_add*int_screen, 
                                          (ct_tests_screen_add - 1)*int_screen + int_conf)
      )]
      
      # Update screen age, reset to NA if after stop date, and reset temporary variables
      m_patients[screen_age >= time_screen_stop, screen_age := NA]
      m_patients[, `:=` (ct_tests_screen_before_FP = NULL,
                         ct_tests_screen_add = NULL)]
      
      # Recalculate number of healthy patients remaining to screen
      n_healthy_screen <- m_patients[screen_age < time_screen_end_H, .N]
    }
    
    # Set number of FP as number of confirmatory tests
    m_patients[, ct_tests_conf := ct_tests_screen_FP]
  }
  
  # Rename test variables
  if ("ct_tests_conf" %in% names(m_patients)) {
    setnames(m_patients, "ct_tests_conf", paste0("ct_tests_", mod_conf))
  }
  setnames(m_patients, "ct_tests_screen", paste0("ct_tests_", mod))
  
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
  }
  
  # Ensure that all modalities are represented in patient data counts
  for (mod in unique(c(v_mod, v_mod_conf[!is.na(v_mod_conf)]))) {
    if (!paste0("ct_tests_", mod) %in% names(m_patients)) {
      m_patients[, paste0("ct_tests_", mod) := 0]
    }
  }
  
  # Set patient and lesion IDs as keys
  setkey(m_lesions, pt_id, lesion_id)
  
  # Merge first screen age within lesion state and minimum of cancer 
  # onset and screening stop date as lesion screening censor date (to be updated after every screening)
  m_lesions[m_patients[screen_age >= time_H_L & 
                         screen_age < pmin(time_H_P, time_screen_stop, na.rm = T)], `:=` (
                           time_lesion_censor = pmin(i.time_H_P, i.time_screen_stop, na.rm = T),
                           screen_type = i.screen_type,
                           screen_age = i.screen_age)]
  
  # Initialize flag for whether lesion was removed
  m_lesions[!is.na(screen_age), `:=` (fl_removed = 0)]
  
  # Get number of people with lesions available for screening
  n_lesion_screen <- m_patients[!is.na(screen_age), .N]
  
  # Loop through screening tests until there are no more lesions to screen
  if (verbose) print("Number of lesions remaining to screen:")
  while (n_lesion_screen > 0) {
    if (verbose) print(n_lesion_screen)
    
    ###### 2.0 Flag eligible screeners and lesions, and assign modalities
    # Flag individuals eligible for screening at current screening time
    m_lesions[screen_age < time_lesion_censor, fl_screen := 1]
    
    # Flag lesions that are present in eligible screeners at current screen time
    # and initialize positive flag and confirmatory test flag to 0
    m_lesions[fl_screen == 1,
              `:=` (fl_present = (time_H_Lj <= screen_age & fl_removed == 0),
                    fl_positive = 0,
                    fl_conf = 0)]
    
    # Set modality based on routine vs. surveillance screening
    m_lesions[!is.na(screen_age), modality := v_mod[screen_type], by = screen_type]
    
    ###### 2.1 Flag whether eligible lesions would be detected
    #### Sample whether eligible lesions produce positive test result
    if (any(m_lesions$fl_present > 0, na.rm = T)) {
      m_lesions[fl_present == 1, `:=` (
        fl_positive = rbinom(
          .N,
          size = 1,
          prob = l_test_chars[[modality]][["p_sens"]][["L"]]
        )), by = modality]
    }
    
    #### Apply downstream effect of positive result depending on test type (direct, targeted, indirect) ####
    ### If test is direct, positive lesions are detected and removed
    if ("direct" %in% names(l_test_types)) {
      m_lesions[fl_present == 1 & fl_positive == 1 & modality %in% l_test_types[["direct"]], `:=` (
        fl_removed = 1
      )]
    }
    
    ### If test is targeted, apply confirmatory tests only to positive lesions to sample whether they are removed
    if ("targeted" %in% names(l_test_types) & m_lesions[fl_present == 1 & fl_positive == 1 & modality %in% l_test_types[["targeted"]], .N] > 0) {
      m_lesions[fl_present == 1 & fl_positive == 1 & modality %in% l_test_types[["targeted"]], `:=` (
        fl_conf = 1,
        fl_removed = rbinom(
          .N,
          size = 1,
          prob = l_test_chars[[v_mod_conf[screen_type]]][["p_sens"]][["L"]]
        )), by = screen_type]
    }
    
    ### If test is indirect (non-targeted), apply confirmatory tests to all lesions if at least one produces positive result
    if ("indirect" %in% names(l_test_types) & m_lesions[fl_present == 1 & modality %in% l_test_types[["indirect"]], .N] > 0) {
      # Check whether any lesion led to positive screening test
      m_detect <- m_lesions[fl_present == 1 & modality %in% l_test_types[["indirect"]],
                            .(fl_positive_any = max(fl_positive)),
                            by = pt_id]
      
      # Merge positive screening test flag
      m_lesions[m_detect, fl_positive_any := i.fl_positive_any]
      
      # Sample removal with confirmatory test
      m_lesions[fl_present == 1 & fl_positive_any == 1,`:=` (
        fl_conf = 1,
        fl_removed = rbinom(
          .N,
          size = 1,
          prob = l_test_chars[[v_mod_conf[screen_type]]][["p_sens"]][["L"]]
        )), by = screen_type]
    }
    
    # For removed lesions, set time to cancer onset to Inf, and set detection time
    m_lesions[fl_present == 1 & fl_removed == 1, `:=` (
      time_H_Pj = Inf,
      time_detected = screen_age)]
    
    # Check number of eligible lesions and removed lesions among individuals screened at this round
    # and update time to onset of preclinical cancer
    m_lesions_summary <- m_lesions[fl_screen == 1, 
                                   .(screen_type = screen_type[1], # Screening type (routine vs. surveillance)
                                     ct_eligible = sum(fl_present), # Number of lesions present
                                     fl_positive = max(fl_positive), # Whether any lesions produced a positive screen result
                                     fl_conf = max(fl_conf), # Whether confirmatory test was conducted
                                     ct_removed = sum(fl_present == 1 & fl_removed == 1), # Count number of lesions removed during screening round
                                     time_H_P = min(time_H_Pj, na.rm = T)), 
                                   by = pt_id]
    
    # Set primary modality again for summarized data
    m_lesions_summary[, modality := v_mod[screen_type], by = screen_type]
    
    ###### 2.2 Flag false positives among screened individuals with no active lesions
    m_lesions_summary[ct_eligible == 0, fl_positive := rbinom(
      .N,
      size = 1,
      prob = 1 - l_test_chars[[modality]][["p_spec"]]),
      by = modality]
    
    # Add confirmatory test for false positives if applicable
    if ("targeted" %in% names(l_test_types) | "indirect" %in% names(l_test_types)) {
      m_lesions_summary[ct_eligible == 0 & fl_positive == 1 & modality %in% unlist(l_test_types[c("targeted", "indirect")]),
                        `:=` (fl_conf = 1)]
    }
    
    ###### 2.3 Update screening data
    #### Set confirmatory modality if applicable
    m_lesions_summary[fl_conf %in% 1, modality_conf:= v_mod_conf[screen_type], by = screen_type]
    
    #### Set next screen time based on screening type, whether there was confirmatory test applied, and number of lesions detected
    if (!is.null(l_params_surveil)) {
      # For routine screening finding no lesions, set time depending on whether confirmatory test was applied
      if (m_lesions_summary[ct_removed == 0 & screen_type == "screen", .N] > 0) {
        m_lesions_summary[ct_removed == 0 & screen_type == "screen", `:=` ( 
          int_test_next = ifelse(!fl_conf %in% 1, 
                                 l_params_strategy[["int_screen"]],
                                 l_params_strategy[["int_conf"]]))]
      }
      
      # For surveillance finding no lesions, return to routine screening using confirmatory test interval as screening interval
      if (m_lesions_summary[ct_removed == 0 & screen_type == "surveil", .N] > 0) {
        m_lesions_summary[ct_removed == 0 & screen_type == "surveil", `:=` (
          screen_type = "screen",
          int_test_next = ifelse(!fl_conf %in% 1, 
                                 int_surveil_neg["screen"],
                                 int_surveil_neg["conf"]))]
      }
      
      # For any screening finding lesions, set to surveillance and set time based on number of lesions detected
      if (m_lesions_summary[ct_removed > 0, .N] > 0) {
        m_lesions_summary[ct_removed > 0, `:=` (
          screen_type = "surveil",
          int_test_next = eval(parse(text = expr_surveil))
        )]
      }
    } else {
      # For routine screening, set time depending on whether confirmatory test was applied
      m_lesions_summary[, `:=` ( 
        int_test_next = ifelse(!fl_conf %in% 1, 
                               l_params_strategy[["int_screen"]],
                               l_params_strategy[["int_conf"]]))]
    }
    
    # Update time to preclinical cancer and following times, increment number of screens, number of confirmatory tests, and screen age
    m_patients[m_lesions_summary, `:=` (time_H_P = i.time_H_P,
                                        screen_type = i.screen_type,
                                        screen_age = screen_age + i.int_test_next,
                                        modality = i.modality,
                                        modality_conf = i.modality_conf)]
    
    # Reset screening age to NA if after censor date or end age
    m_patients[screen_age >= time_screen_stop, screen_age := NA]
    
    # Add modality counts and delete flags
    for (mod in unique(v_mod)) {
      colname <- paste0("ct_tests_", mod)
      m_patients[modality == mod, (colname) := get(colname) + 1]
    }
    for (mod in unique(v_mod_conf[!is.na(v_mod_conf)])) {
      colname <- paste0("ct_tests_", mod)
      m_patients[modality_conf == mod, (colname) := get(colname) + 1]
    }
    m_patients[, `:=` (modality = NULL, modality_conf = NULL)]
    
    # Reset flags for eligible screeners and lesions
    cols_remove <- c("screen_age", "fl_screen", "fl_present", "fl_positive", 
                     "fl_positive_any", "modality", "modality_conf", "fl_conf")
    m_lesions[, (intersect(cols_remove, names(m_lesions))) := NULL]
    
    # Merge screen age within lesion state and minimum of cancer 
    # onset and death as lesion screening censor date (to be updated after every screening)
    m_lesions[m_patients[screen_age >= time_H_L & 
                           screen_age < pmin(time_H_P, time_screen_stop, na.rm = T)], `:=` (
                             time_lesion_censor = pmin(i.time_H_P, i.time_screen_stop, na.rm = T),
                             screen_type = i.screen_type,
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
  }
  
  # Run while loop
  if (verbose) print("Number of individuals with preclinical cancer remaining to screen:")
  while (m_patients[!is.na(screen_age), .N] > 0) {
    if (verbose) print(m_patients[!is.na(screen_age), .N])
    # Set modality based on routine vs. surveillance screening
    m_patients[!is.na(screen_age), modality := v_mod[screen_type], by = screen_type]
    
    # Sample whether cancer would produce positive test result
    m_patients[!is.na(screen_age), `:=` (
      fl_positive = rbinom(
        .N,
        size = 1,
        prob = l_test_chars[[modality]][["p_sens"]][["P"]]
      )), by = modality]
    
    #### Apply downstream effect of positive result depending on test type (direct, targeted, indirect) ####
    ### If test is direct, positive cancer is detected
    if ("direct" %in% names(l_test_types)) {
      m_patients[!is.na(screen_age) & fl_positive == 1 & modality %in% l_test_types[["direct"]], `:=` (
        fl_conf = 0,
        fl_detected = 1
      )]
    }
    
    ### If test is targeted or indirect, apply confirmatory test to sample whether it is detected
    if ("targeted" %in% names(l_test_types) | "indirect" %in% names(l_test_types)) {
      if (m_patients[!is.na(screen_age) & fl_positive == 1 & modality %in% unlist(l_test_types[c("targeted", "indirect")]), .N] > 0) {
        m_patients[!is.na(screen_age) & fl_positive == 1 & modality %in% unlist(l_test_types[c("targeted", "indirect")]), `:=` (
          fl_conf = 1,
          modality_conf = v_mod_conf[screen_type],
          fl_detected = rbinom(
            .N,
            size = 1,
            prob = l_test_chars[[v_mod_conf[screen_type]]][["p_sens"]][["P"]]
          )), by = screen_type]
      }
    }
    
    # Add primary modality counts and delete flags
    for (mod in unique(v_mod)) {
      colname <- paste0("ct_tests_", mod)
      m_patients[modality == mod, (colname) := get(colname) + 1]
    }
    m_patients[, `:=` (modality = NULL)]
    
    # Add confirmatory modality counts and delete flags
    if ("modality_conf" %in% names(m_patients)) {
      for (mod in unique(v_mod_conf[!is.na(v_mod_conf)])) {
        colname <- paste0("ct_tests_", mod)
        m_patients[modality_conf == mod, (colname) := get(colname) + 1]
      }
      m_patients[, `:=` (modality_conf = NULL)]
    }
    
    # Set next test time for negatives cases that were routinely screened
    m_patients[, int_test_next := NA_integer_]
    m_patients[!is.na(screen_age) & !fl_detected %in% 1 & screen_type == "screen", `:=` (
      int_test_next = ifelse(!fl_conf %in% 1, 
                             l_params_strategy[["int_screen"]],
                             l_params_strategy[["int_conf" ]])
    )]
    
    # Set next test time for negatives cases that were surveillance-screened
    m_patients[!is.na(screen_age) & !fl_detected %in% 1 & screen_type == "surveil", `:=` (
      screen_type = "screen",
      int_test_next = ifelse(!fl_conf %in% 1, 
                             int_surveil_neg["screen"],
                             int_surveil_neg["conf"])
    )]
    
    # Reset clinical cancer age and set next screen age to NA for newly detected cases
    m_patients[!is.na(screen_age) & fl_detected == 1, `:=` (
      time_H_C = screen_age,
      time_P_C = screen_age - time_H_P,
      screen_age = NA)]
    
    # Update screen age, reset to NA if after stop date, and reset temporary variables
    m_patients[!is.na(screen_age), screen_age := screen_age + int_test_next]
    m_patients[screen_age >= time_screen_stop, screen_age := NA]
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
    browser()
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