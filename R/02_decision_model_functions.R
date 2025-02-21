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
  m_cohort_base <- initialize_cohort(l_params_model)
  
  # Simulate baseline characteristics
  simulate_baseline_data(m_cohort_base, l_params_model)
  
  # Simulate disease (lesion and/or cancer) onset
  simulate_disease_onset(m_cohort_base, l_params_model)
  
  # Simulate additional lesion onset and lesion progression
  if("L" %in% l_params_model$v_states) {
    m_lesions <- simulate_additional_lesions(m_cohort_base, l_params_model)
    simulate_lesion_progression(m_cohort_base, m_lesions, l_params_model)
  }
  
  # Simulate cancer progression
  simulate_cancer_progression(m_cohort_base, l_params_model)
  
  # Simulate cancer mortality
  simulate_cancer_mortality(m_cohort_base, l_params_model)
  
  # Compile overall mortality outcomes from background and cancer data
  calc_mortality_outcomes(m_cohort_base)
  
  # Wrap no screening results
  if("L" %in% l_params_model$v_states) {
    res <- list(patient_level = m_cohort_base,
                lesion_level = m_lesions)
  } else {
    res <- m_cohort_base
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
      var_detect <- paste0("time_P", v_cancer[i], "_C")
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
      var_detect <- paste0("time_P", v_cancer[i], "_C")
      
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
    var_detect <- paste0("time_P", tail(v_cancer, 1), "_C")
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
  
  if (length(idx) > 0) {
    with(l_params_model, {
      m_patients[idx, time_C_Dc := query_distr(
        "r", .N, 
        d_time_C_Dc[[stage_dx]]$distr, 
        d_time_C_Dc[[stage_dx]]$params
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


# Rerun version for screening
run_screening_counterfactual <- function(
    m_cohort,
    l_params_model,
    age_screen_start,
    age_screen_end,
    int_screen,
    p_sens,
    p_spec,
    overwrite = FALSE, # If TRUE, allows data in m_cohort to be overwritten; otherwise saves and returns a copy of original data
    verbose = FALSE
) {
  # Save original data if overwrite = FALSE
  if (!overwrite) {
    m_cohort_orig <- copy(m_cohort)
  }
  
  # Separate patient and lesion data as necessary
  if (is.data.table(m_cohort)) {
    m_patients <- m_cohort
  } else {
    m_patients <- m_cohort$patient_level
    m_lesions <- m_cohort$lesion_level
  }
  
  #### Case 1: Simulate screening before disease ####
  # Set variable for time to disease onset
  var_onset <- paste("time", l_params_model$v_states[1], l_params_model$v_states[2], sep = "_")
  
  # Simulate screening during healthy state
  simulate_screening_H(m_patients = m_patients,
                       var_onset = var_onset,
                       age_screen_start = age_screen_start,
                       age_screen_end = age_screen_end,
                       int_screen = int_screen,
                       p_spec = p_spec)
  
  #### Case 2: Lesions developed by screen age, but no preclinical cancer ####
  if ('L' %in% l_params_model$v_states) {
    # Simulate lesion progression to cancer
    simulate_screening_L(m_patients = m_patients, 
                         m_lesions = m_lesions,
                         l_params_model = l_params_model,
                         age_screen_end = age_screen_end,
                         int_screen = int_screen,
                         p_sens = p_sens["L"],
                         p_spec = p_spec,
                         verbose = verbose)
  }
  
  #### Case 3: Preclinical cancer screening and surveillance ####
  simulate_screening_P(m_patients = m_patients,
                       l_params_model = l_params_model,
                       age_screen_end = age_screen_end,
                       int_screen = int_screen,
                       p_sens = p_sens["P"],
                       verbose = verbose)
  
  # Recalculate mortality outcomes
  calc_mortality_outcomes(m_patients)
  
  # Sum tests requiring diagnostic workup
  m_patients[, ct_diag_total := rowSums(.SD, na.rm = T), .SDcols = patterns("ct_diag_")]
  
  # Save original data if overwrite = FALSE
  if (!overwrite) {
    return(list(m_cohort = m_cohort, m_cohort_base = m_cohort_orig))
  } else {
    return(NULL)
  }
}


# Simulate screening occurring during the healthy state
simulate_screening_H <- function(m_patients,
                                 var_onset,
                                 age_screen_start,
                                 age_screen_end,
                                 int_screen,
                                 p_spec
) {
  # Get routine screening ages
  v_screen_ages <- seq(age_screen_start, age_screen_end, int_screen)
  
  # Get number of tests during healthy state before earliest of death or disease onset
  m_patients[, ct_screen := findInterval(pmin(time_H_D, get(var_onset)), v_screen_ages, left.open = T)]
  
  # Simulate number of false positives during healthy state, leading to diagnostic workup
  m_patients[, `:=` (
    ct_diag_FP = rbinom(
      .N,
      size = ct_screen,
      prob = 1 - p_spec
    ))]
  
  # Get first screening age with disease present
  m_patients[, screen_age := v_screen_ages[ct_screen + 1], by = ct_screen]
  
  # Reset first screening age to NA if after censor date
  m_patients[screen_age >= time_screen_censor, screen_age := NA]
  
  return(NULL)
}


# Simulate screening occurring during the precancerous lesion state
# Assumes that screening during healthy state has been completed for population,
# which initializes the screen_age variable
simulate_screening_L <- function(m_patients,
                                 m_lesions,
                                 var_onset,
                                 l_params_model,
                                 age_screen_end,
                                 int_screen,
                                 p_sens,
                                 p_spec,
                                 verbose = verbose
) {
  # Merge first screen age within lesion state and minimum of cancer 
  # onset and death as lesion screening censor date (to be updated after every screening)
  m_lesions[m_patients[screen_age < time_H_P], `:=` (time_lesion_censor = pmin(time_H_P, time_H_Do, na.rm = T),
                                                     screen_age = i.screen_age)]
  
  # Among patients who will receive screening during the lesion state,
  # initialize count for number of screens and flag whether lesion was removed
  m_lesions[!is.na(screen_age), `:=` (ct_screen = 0,
                                      fl_removed = 0)]
  
  # Get number of lesions available for screening
  n_lesion_screen <- m_lesions[!is.na(screen_age), .N]
  
  # Initialize patient-level variables to track number of positive screens and lesions detected
  m_patients[, `:=` (ct_diag_TP_L = 0,
                     ct_lesion_detected = 0)]
  
  # Loop through screening tests until there are no more lesions to screen
  if (verbose) print("Number of lesions remaining to screen:")
  while (n_lesion_screen > 0) {
    if (verbose) print(n_lesion_screen)
    
    ###### 2.1 Flag whether eligible lesions would be detected
    # Initialize flag for eligible screeners and lesions
    m_lesions[, `:=` (fl_screen = NA,
                      fl_eligible = NA)]
    
    # Flag individuals eligible for screening - 
    # Defined as patient having not developed preclinical cancer yet
    m_lesions[time_lesion_censor > screen_age, fl_screen := 1]
    
    # Flag lesions that are eligible to be detected at current screen time - 
    # Onset occurred before screen age and not removed yet
    m_lesions[time_H_L + time_L_Lj <= screen_age & fl_removed == 0,
              fl_eligible := 1]
    
    # Sample whether eligible lesions produce positive test result
    m_lesions[fl_eligible == 1, `:=` (
      fl_positive = rbinom(
        .N,
        size = 1,
        prob = p_sens
      ))]
    
    # Apply downstream effect of positive result - 
    # Assume that lesions with positive result are removed,
    # reset time to cancer onset reset to Inf, and set detection time
    m_lesions[fl_eligible == 1 & fl_positive == 1, `:=` (fl_removed = 1,
                                                         time_H_Pj = Inf,
                                                         time_detected = screen_age)]
    
    # Check number of eligible lesions and removed lesions among individuals screened at this round
    # and update time to onset of preclinical cancer
    m_lesions_summary <- m_lesions[fl_screen == 1, 
                                   .(ct_eligible = sum(fl_eligible, na.rm = T),
                                     ct_removed = sum(fl_eligible == 1 & fl_removed == 1, na.rm = T),
                                     time_H_P = min(time_H_Pj, na.rm = T)), 
                                   by = pt_id]
    
    # Sample false positives among people in lesion state without any active lesions
    m_lesions_summary[ct_eligible == 0, fl_FP := rbinom(
      .N,
      size = 1,
      prob = 1 - p_spec)]
    
    # Update time to preclinical cancer and following times;
    # increment number of screens, positive screens, removed lesions, false positives, and screen age
    m_patients[m_lesions_summary, `:=` (time_H_P = i.time_H_P,
                                        ct_screen = ct_screen + 1,
                                        ct_diag_TP_L = ct_diag_TP_L + (i.ct_removed > 0),
                                        ct_lesion_detected = ct_lesion_detected + i.ct_removed,
                                        ct_diag_FP = ct_diag_FP + pmax(0, i.fl_FP, na.rm = T),
                                        screen_age = screen_age + int_screen)]
    
    # Reset first screening age to NA if after censor date or end age
    m_patients[screen_age >= time_screen_censor | screen_age > age_screen_end, screen_age := NA]
    
    # Merge first screen age within lesion state and minimum of cancer 
    # onset and death as lesion screening censor date (to be updated after every screening)
    m_lesions[, `:=` (time_lesion_censor = NULL,
                      screen_age = NULL)]
    m_lesions[m_patients[screen_age < time_H_P], `:=` (time_H_P = i.time_H_P,
                                                       time_lesion_censor = pmin(i.time_H_P, time_H_Do, na.rm = T),
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
                                 age_screen_end,
                                 int_screen,
                                 p_sens,
                                 verbose = FALSE
) {
  if (verbose) print("Number of individuals with preclinical cancer remaining to screen:")
  while (m_patients[!is.na(screen_age), .N] > 0) {
    if (verbose) print(m_patients[!is.na(screen_age), .N])
    # Increment number of screening tests and sample whether cancer leads to positive screen
    m_patients[time_H_P <= screen_age, `:=` (ct_screen = ct_screen + 1,
                                             fl_detected = rbinom(
                                               .N,
                                               size = 1,
                                               prob = p_sens
                                             ))]
    
    # Reset clinical cancer age and set next screen age to NA for detected cases
    m_patients[!is.na(screen_age) & fl_detected == 1, `:=` (
      time_H_C = screen_age,
      time_P_C = screen_age - time_H_P,
      screen_age = NA)]
    
    # Update screen age and set to NA if beyond censor age or maximum screen age
    m_patients[!is.na(screen_age), screen_age := screen_age + int_screen]
    m_patients[screen_age >= time_screen_censor | screen_age > age_screen_end, screen_age := NA]
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
      d_time_C_Dc[[stage_dx]]$distr,
      d_time_C_Dc[[stage_dx]]$params
    ), by = stage_dx]
  })
  
  # Flag cancers diagnosed with screening vs symptoms
  m_patients[fl_detected == 1, ct_diag_TP_P := 1]
  m_patients[!fl_detected %in% 1 & time_H_C < time_H_D, ct_diag_C := 1]
  
  return(NULL)
}
