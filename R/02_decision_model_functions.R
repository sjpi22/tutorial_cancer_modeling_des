#' Generate time to event data for base case (no screening) decision model
#'
#' @param l_params_all List with all parameters of decision model
#' 
#' @return Data table of time to event data for each individual
#' 
#' @export
run_base_model <- function(l_params_all) {
  # Set seed
  if(!is.null(l_params_all$seed)) {
    set.seed(l_params_all$seed)
  }
  
  # Initialize matrix of patient data
  m_cohort_base <- initialize_cohort(l_params_all)
  
  # Simulate baseline characteristics
  simulate_baseline_data(m_cohort_base, l_params_all)
  
  # Simulate disease (lesion and/or cancer) onset
  simulate_disease_onset(m_cohort_base, l_params_all)
  
  # Simulate additional lesion onset and lesion progression
  if("L" %in% l_params_all$v_states) {
    m_lesion <- simulate_additional_lesions(m_cohort_base, l_params_all)
    simulate_lesion_progression(m_cohort_base, m_lesion, l_params_all)
  }
  
  # Simulate cancer progression
  simulate_cancer_progression(m_cohort_base, l_params_all)
  
  # Simulate cancer mortality
  simulate_cancer_mortality(m_cohort_base, l_params_all)
  
  # Compile overall mortality outcomes from background and cancer data
  calc_mortality_outcomes(m_cohort_base)
  
  # Wrap no screening results
  if("L" %in% l_params_all$v_states) {
    res <- list(patient_level = m_cohort_base,
                lesion_level = m_lesion)
  } else {
    res <- m_cohort_base
  }
  
  return(res)
}


# Initialize time to event matrix for patient cohort with patient IDs and baseline strategy
initialize_cohort <- function(l_params_all) {
  m_times_init <- with(as.list(l_params_all), {
    m_times <- data.table(
      pt_id = 1:n_cohort
    )
    setkey(m_times, pt_id)
    return(m_times)
  }
  )
  return(m_times_init)
}


# Simulate baseline characteristics
simulate_baseline_data <- function(m_times, l_params_all) {
  with(l_params_all, {
    # Sample time to death from other causes for multiple sexes
    if (length(sex) > 1) {
      # Sample male / female
      m_times[, male := query_distr(
        "r", .N, 
        d_male$distr, 
        d_male$params
      )]
      
      # Sample time to death from other causes if male
      m_times[male == 1, time_H_Do := query_distr(
        "r", .N, 
        d_time_H_Do[["male"]]$distr, 
        d_time_H_Do[["male"]]$params
      )]
      
      # Sample time to death from other causes if female
      m_times[male == 0, time_H_Do := query_distr(
        "r", .N, 
        d_time_H_Do[["female"]]$distr, 
        d_time_H_Do[["female"]]$params
      )]
    } else { # Sample time to death from other causes for single sex case
      # Assign male / female
      m_times[, male := (sex == "male")]
      
      # Sample time to death from other causes if female
      m_times[, time_H_Do := query_distr(
        "r", .N, 
        d_time_H_Do[[sex]]$distr, 
        d_time_H_Do[[sex]]$params
      )]
    }
  })
  return(NULL)
}


# Simulate disease (lesion or cancer) onset
simulate_disease_onset <- function(m_times, l_params_all) {
  with(l_params_all, {
    # If including lesion state, simulate time to lesion and cancer onset
    if("L" %in% v_states) {
      # Simulate time to first lesion onset
      m_times[, time_H_L := query_distr(
        "r", .N, 
        d_time_H_L$distr, 
        d_time_H_L$params
      )]
      
      # Simulate number of additional lesions (only for individuals with lesion onset before death)
      # Scale rate down by mean screening duration, then stretch by duration until death from other causes
      m_times[time_H_Do > time_H_L, n_L_add := query_distr(
        "r", .N, 
        d_n_L$distr, 
        lapply(d_n_L$params, '*', (time_H_Do - time_H_L))
      )]
      
      # Cap number of additional lesions
      m_times[n_L_add >= n_lesions_max, n_L_add := n_lesions_max - 1]
    } else {
      # Simulate time to de novo cancer onset
      m_times[, time_H_P := query_distr("r", .N, d_time_H_P$distr, d_time_H_P$params)]
    }
  })
  return(NULL)
}


# Simulate development of additional lesions and output lesion-level data table
simulate_additional_lesions <- function(m_times, l_params_all) {
  m_lesion <- with(as.list(l_params_all), {
    # Subset to individuals with lesion onset before death
    m_lesion <- m_times[time_H_L < time_H_Do]
    
    # Account for case with no lesions
    if (nrow(m_lesion) > 0) {
      # Set initial lesion ID and individual onset time variable
      m_lesion[, `:=` (lesion_id = 1, time_L_Lj = 0)]
      
      # Create new row for each additional lesion (account for case where no one has >1 lesion)
      if (nrow(m_lesion[n_L_add > 0]) > 0) {
        m_other_lesion <- m_lesion[n_L_add > 0, .(
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
        m_lesion <- rbindlist(list(
          m_lesion[, .SD, .SDcols = intersect(names(m_lesion), names(m_other_lesion))],
          m_other_lesion
        ))
      }
      
      # Set pt_id as key
      setkey(m_lesion, pt_id)
    }
    return(m_lesion)
  })
  return(m_lesion)
}


# Simulate lesion progression to cancer
simulate_lesion_progression <- function(m_times, m_lesion, l_params_all) {
  with(as.list(l_params_all), {
    # Account for case of no lesions
    if (!is.null(m_lesion)) {
      # Sample time to cancer
      m_lesion[, time_Lj_Pj := query_distr(
        "r", .N, 
        d_time_L_P$distr, 
        d_time_L_P$params)]
      
      # Calculate time from first lesion to lesion-specific cancer onset
      m_lesion[, time_L_Pj := time_L_Lj + time_Lj_Pj]
      
      # Calculate time from healthy to lesion-specific cancer onset
      m_lesion[, time_H_Pj := time_H_L + time_L_Pj]
      
      # Get earliest time to cancer onset
      m_times_cancer <- m_lesion[, .(time_H_P = min(time_H_Pj)), by = pt_id]
      
      # Merge to patient data
      m_times[m_times_cancer, time_H_P := i.time_H_P]
    } else {
      m_times[, time_H_P := Inf]
    }
  })
  return(NULL)
}


# Simulate cancer stage progression by stage
simulate_cancer_progression <- function(m_times, l_params_all) {
  with(l_params_all, {
    # Populate time to progression to next preclinical stage
    for (i in 1:(length(v_cancer)-1)) {
      var_progress <- paste0("time_P", v_cancer[i], "_P", v_cancer[i+1])
      m_times[time_H_P < time_H_Do, (var_progress) := query_distr(
        "r", .N, 
        get(paste0("d_", var_progress))$distr, 
        get(paste0("d_", var_progress))$params
      )]
    }
    
    # Populate time to detection (clinical cancer)
    for (i in 1:length(v_cancer)) {
      var_detect <- paste0("time_P", v_cancer[i], "_C")
      m_times[time_H_P < time_H_Do, (var_detect) := query_distr(
        "r", .N, 
        get(paste0("d_", var_detect))$distr, 
        get(paste0("d_", var_detect))$params
      )]
    }
    
    # Initialize stage of diagnosis variable
    m_times[, stage_dx := eval(parse(text = paste0("as.", class(v_cancer), "(NA)")))]
    
    # Initialize time from cancer onset to detection (should be NA for patients without cancer onset before death from other causes)
    m_times[time_H_P < time_H_Do, time_P_C := 0]
    
    # Loop through stages until detection occurs
    for (i in 1:(length(v_cancer)-1)) {
      # Get progression and detection variables
      var_progress <- paste0("time_P", v_cancer[i], "_P", v_cancer[i+1])
      var_detect <- paste0("time_P", v_cancer[i], "_C")
      
      # If stage at diagnosis has not been set yet and detection occurs before 
      # progression, add time to detection to running total of time from cancer 
      # onset to detection and set stage at diagnosis
      m_times[time_H_P < time_H_Do & get(var_detect) < get(var_progress) & is.na(stage_dx), `:=` (
        stage_dx = v_cancer[i],
        time_P_C = time_P_C + get(var_detect))]
      
      # Otherwise if not detected yet, add time to detection to running total 
      # of time from cancer onset to detection
      m_times[time_H_P < time_H_Do & is.na(stage_dx), 
              time_P_C := time_P_C + get(var_progress)]
    }
    
    # If not detected yet by second to last stage, set last stage as stage at detection
    var_detect <- paste0("time_P", tail(v_cancer, 1), "_C")
    m_times[time_H_P < time_H_Do & is.na(stage_dx), `:=` (
      stage_dx = tail(v_cancer, 1),
      time_P_C = time_P_C + get(var_detect))]
  })
  
  # Set time to cancer diagnosis
  m_times[!is.na(time_P_C), time_H_C := time_H_P + time_P_C]
  return(NULL)
}


# Simulate cancer mortality by stage at diagnosis among people with cancer onset before death from other causes
simulate_cancer_mortality <- function(m_times, l_params_all) {
  # Get indices of patients with cancer onset before death from other causes
  # Not cancer diagnosis before death from other causes to be more conservative, in case screening may convert from preclinical to clinical before death
  idx <- m_times[time_H_P < time_H_Do, which = TRUE]
  
  if (length(idx) > 0) {
    with(l_params_all, {
      m_times[idx, time_C_Dc := query_distr(
        "r", .N, 
        d_time_C_Dc[[stage_dx]]$distr, 
        d_time_C_Dc[[stage_dx]]$params
      ), by = stage_dx]
    })
    
    # Calculate death from cancer
    m_times[idx, time_H_Dc := time_H_C + time_C_Dc]
  } else {
    m_times[, `:=`(time_C_Dc = NA,
                   time_H_Dc = NA)]
  }
  return(NULL)
}


# Generate mortality outcomes
calc_mortality_outcomes <- function(m_times) {
  # Join cancer data to patient-level data
  # Calculate all-cause death and cause of death
  m_times[, time_H_D := pmin(time_H_Do, time_H_Dc, na.rm = TRUE)]
  m_times[, fl_Dc := (time_H_Do > pmin(time_H_Dc, Inf, na.rm = TRUE))]
  
  # Calculate survival from cancer diagnosis
  m_times[time_H_C <= time_H_D, time_C_D := time_H_D - time_H_C]
  
  # Set censor age for screening (earliest of cancer diagnosis or death)
  m_times[, time_screen_censor := pmin(time_H_C, time_H_D, na.rm = T)]
  
  return(NULL)
}


# Rerun version for screening (note: overwrites data by updating by reference)
run_screening_counterfactual <- function(
    age_screen_start,
    age_screen_end,
    int_screen,
    p_sens,
    p_spec,
    l_params_all,
    m_times,
    m_lesion = NULL,
    verbose = FALSE
) {
  
  #### Case 1: Simulate screening before disease ####
  # Set variable for time to disease onset
  var_onset <- paste("time", l_params_all$v_states[1], l_params_all$v_states[2], sep = "_")
  
  # Simulate screening during healthy state
  simulate_screening_H(m_times = m_times,
                       var_onset = var_onset,
                       age_screen_start = age_screen_start,
                       age_screen_end = age_screen_end,
                       int_screen = int_screen,
                       p_spec = p_spec)
  
  
  #### Case 2: Lesions developed by screen age, but no preclinical cancer ####
  if ('L' %in% l_params_all$v_states) {
    # Simulate lesion progression to cancer
    simulate_screening_L(m_times = m_times, 
                         m_lesion = m_lesion,
                         l_params_all = l_params_all,
                         age_screen_end = age_screen_end,
                         int_screen = int_screen,
                         p_sens = p_sens[["L"]],
                         p_spec = p_spec,
                         verbose = verbose)
  }
  
  #### Case 3: Preclinical cancer screening and surveillance ####
  simulate_screening_P(m_times = m_times,
                       l_params_all = l_params_all,
                       age_screen_end = age_screen_end,
                       int_screen = int_screen,
                       p_sens = p_sens[["P"]],
                       verbose = verbose)

  # Recalculate mortality outcomes
  calc_mortality_outcomes(m_times)
  
  return(NULL)
}


# Simulate screening occurring during the healthy state
simulate_screening_H <- function(m_times,
                                 var_onset,
                                 age_screen_start,
                                 age_screen_end,
                                 int_screen,
                                 p_spec
) {
  # Get routine screening ages
  v_screen_ages <- seq(age_screen_start, age_screen_end, int_screen)
  
  # Get number of tests during healthy state before earliest of death or disease onset
  m_times[, ct_screen := findInterval(pmin(time_H_D, get(var_onset)), v_screen_ages, left.open = T)]
  
  # Simulate number of false positives during healthy state
  m_times[, `:=` (
    ct_FP = rbinom(
      .N,
      size = ct_screen,
      prob = 1 - p_spec
    ))]
  
  # Get first screening age with disease present
  m_times[, screen_age := v_screen_ages[ct_screen + 1], by = ct_screen]
  
  # Reset first screening age to NA if after censor date
  m_times[screen_age >= time_screen_censor, screen_age := NA]
  
  return(NULL)
}


# Simulate screening occurring during the precancerous lesion state
# Assumes that screening during healthy state has been completed for population,
# which initializes the screen_age variable
simulate_screening_L <- function(m_times,
                                 m_lesion,
                                 var_onset,
                                 l_params_all,
                                 age_screen_end,
                                 int_screen,
                                 p_sens,
                                 p_spec,
                                 verbose = verbose
) {
  browser()
  # Merge first screen age within lesion state and minimum of cancer 
  # onset and death as lesion screening censor date (to be updated after every screening)
  m_lesion[m_times[screen_age < time_H_P], `:=` (time_lesion_censor = pmin(i.time_H_P, time_H_Do, na.rm = T),
                                                 screen_age = i.screen_age)]
  
  # Among patients who will receive screening during the lesion state,
  # initialize count for number of screens and flag whether lesion was removed
  m_lesion[!is.na(screen_age), `:=` (ct_screen = 0,
                                     fl_removed = 0)]
  
  # Get number of lesions available for screening
  n_lesion_screen <- m_lesion[!is.na(screen_age), .N]
  
  # Initialize patient-level variables to track number of positive screens and lesions detected
  m_times[, `:=` (ct_TP_lesion = 0,
                  ct_lesion_detected = 0)]
  
  # Loop through screening tests until there are no more lesions to screen
  if (verbose) print("Number of lesions remaining to screen:")
  while (n_lesion_screen > 0) {
    if (verbose) print(n_lesion_screen)
    
    ###### 2.1 Flag whether eligible lesions would be detected
    # Initialize flag for eligible screeners and lesions
    m_lesion[, `:=` (fl_screen = NA,
                     fl_eligible = NA)]
    
    # Flag individuals eligible for screening - 
    # Defined as patient having not developed preclinical cancer yet
    m_lesion[time_lesion_censor > screen_age, fl_screen := 1]
    
    # Flag lesions that are eligible to be detected at current screen time - 
    # Onset occurred before screen age and not removed yet
    m_lesion[time_H_L + time_L_Lj <= screen_age & fl_removed == 0,
             fl_eligible := 1]
    
    # Sample whether eligible lesions produce positive test result
    m_lesion[fl_eligible == 1, `:=` (
      fl_positive = rbinom(
        .N,
        size = 1,
        prob = p_sens
      ))]
    
    # Apply downstream effect of positive result - 
    # Assume that lesions with positive result are removed,
    # reset time to cancer onset reset to Inf, and set detection time
    m_lesion[fl_eligible == 1 & fl_positive == 1, `:=` (fl_removed = 1,
                                                        time_H_Pj = Inf,
                                                        time_detected = screen_age)]
    
    # Check number of eligible lesions and removed lesions among individuals screened at this round
    # and update time to onset of preclinical cancer
    m_lesion_summary <- m_lesion[fl_screen == 1, 
                                 .(ct_eligible = sum(fl_eligible, na.rm = T),
                                   ct_removed = sum(fl_eligible == 1 & fl_removed == 1, na.rm = T),
                                   time_H_P = min(time_H_Pj, na.rm = T)), 
                                 by = pt_id]
    
    # Sample false positives among people in lesion state without any active lesions
    m_lesion_summary[ct_eligible == 0, fl_FP := rbinom(
      .N,
      size = 1,
      prob = 1 - p_spec)]
    
    # Update time to preclinical cancer and following times;
    # increment number of screens, positive screens, removed lesions, false positives, and screen age
    m_times[m_lesion_summary, `:=` (time_H_P = i.time_H_P,
                                    ct_screen = ct_screen + 1,
                                    ct_TP_lesion = ct_TP_lesion + (i.ct_removed > 0),
                                    ct_lesion_detected = ct_lesion_detected + i.ct_removed,
                                    ct_FP = ct_FP + i.fl_FP,
                                    screen_age = screen_age + int_screen)]
    
    # Merge first screen age within lesion state and minimum of cancer 
    # onset and death as lesion screening censor date (to be updated after every screening)
    m_lesion[, `:=` (time_lesion_censor = NULL,
                     screen_age = NULL)]
    m_lesion[m_times[screen_age < time_H_P], `:=` (time_lesion_censor = pmin(i.time_H_P, time_H_Do, na.rm = T),
                                                   screen_age = i.screen_age)]
    
    # Update number of lesions available for screening
    n_lesion_screen <- m_lesion[!is.na(screen_age), .N]
  }
  
  # Update mortality outcomes
  m_times[, `:=` (time_H_C = time_H_P + time_P_C,
                  time_H_Dc = time_H_P + time_P_C + time_C_Dc)]
  calc_mortality_outcomes(m_times)
  return(NULL)
}


# Simulate screening occurring during the preclinical cancer state
# Assumes that screening during healthy state has been completed for population,
# which initializes the screen_age variable
simulate_screening_P <- function(m_times,
                                 l_params_all,
                                 age_screen_end,
                                 int_screen,
                                 p_sens,
                                 verbose = FALSE
) {
  
  if (verbose) print("Number of individuals with preclinical cancer remaining to screen:")
  while (m_times[!is.na(screen_age), .N] > 0) {
    if (verbose) print(m_times[!is.na(screen_age), .N])
    
    # Increment number of screening tests and sample whether cancer leads to positive screen
    m_times[time_H_P < screen_age, `:=` (ct_screen = ct_screen + 1,
                                         fl_detected = rbinom(
                                           .N,
                                           size = 1,
                                           prob = p_sens
                                         ))]
    
    # Convert detected cases to clinical cancer, get stage at diagnosis, set next screen age to NA
    # m_times[!is.na(screen_age) & fl_detected == 1, `:=` (
    #   time_H_C = screen_age,
    #   stage_dx = fcase(
    #     time_H_C - get(paste0("time_P", stage_dx, "_C")) <= screen_age, as.double(stage_dx),
    #     time_H_C - get(paste0("time_P", stage_dx, "_C")) -
    #       get(paste0("time_P", stage_dx-1, "_P", stage_dx)) <= screen_age, as.double(stage_dx-1),
    #     time_H_C - get(paste0("time_P", stage_dx, "_C")) -
    #       get(paste0("time_P", stage_dx-1, "_P", stage_dx)) -
    #       get(paste0("time_P", stage_dx-2, "_P", stage_dx-1)) <= screen_age, as.double(stage_dx-2),
    #     time_H_C - get(paste0("time_P", stage_dx, "_C")) -
    #       get(paste0("time_P", stage_dx-1, "_P", stage_dx)) -
    #       get(paste0("time_P", stage_dx-2, "_P", stage_dx-1)) -
    #       get(paste0("time_P", stage_dx-3, "_P", stage_dx-2)) <= screen_age, as.double(stage_dx-3)
    #   ),
    #   screen_age = NA),
    #   by = stage_dx]
    m_times[!is.na(screen_age) & fl_detected == 1, `:=` (
      time_H_C = screen_age,
      time_P_C = screen_age - time_H_P,
      screen_age = NA)]
    
    # Update screen age
    m_times[!is.na(screen_age), screen_age := screen_age + int_screen]
    m_times[screen_age >= time_screen_censor | screen_age > age_screen_end, screen_age := NA]
  }
  
  ##### For patients with detected cancer, recalculate stage at diagnosis and resimulate progression to clinical cancer and mortality outcomes
  # Note: Assumes at least 2 stages of cancer
  m_times[fl_detected == 1, `:=` (stage_dx = 1,
                                  time_P_C_running = time_P1_P2)]
  
  for (stg in 2:length(l_params_all$v_cancer)) {
    m_times[fl_detected == 1 & time_P_C > time_P_C_running, `:=` (stage_dx = stage_dx + 1,
                                                                  time_P_C_running = time_P_C_running + get(paste0("time_P", stg, "_P", stg + 1)))]
  }
  
  # Recalculate time to death from cancer among people with screen-detected cancer
  with(l_params_all, {
    m_times[fl_detected == 1, time_C_Dc := query_distr(
      "r", .N,
      d_time_C_Dc[[v_surv[stage_dx]]]$distr,
      d_time_C_Dc[[v_surv[stage_dx]]]$params,
      d_time_C_Dc[[v_surv[stage_dx]]]$src
    ), by = stage_dx]
  })
  
  return(NULL)
}
