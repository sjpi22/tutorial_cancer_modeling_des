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
    l_df_lesion <- simulate_additional_lesions(m_cohort_base, l_params_all)
    simulate_lesion_progression(m_cohort_base, l_df_lesion, l_params_all)
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
                lesion_level = l_df_lesion)
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
    # If including a lesion state, simulate time to lesion and cancer onset
    if("L" %in% v_states) {
      # By lesion type
      for (i in v_lesion) {
        # Simulate time to first lesion onset
        m_times[, paste0("time_H_L", i) := query_distr(
          "r", .N, 
          get(paste0("d_time_H_L", i))$distr, 
          get(paste0("d_time_H_L", i))$params
        )]
        
        # Simulate number of additional lesions (only for individuals with lesion onset before death)
        # Scale rate down by mean screening duration, then stretch by duration until death from other causes
        m_times[time_H_Do > get(paste0("time_H_L", i)), paste0("n_L", i) := query_distr(
          "r", .N, 
          get(paste0("d_n_L", i))$distr, 
          lapply(get(paste0("d_n_L", i))$params, '*', (time_H_Do - get(paste0("time_H_L", i))))
        )]
        
        # Cap number of lesions
        m_times[get(paste0("n_L", i)) > max_lesions, paste0("n_L", i) := max_lesions]
      }
      
      # Take earliest lesion onset time as time_H_L
      m_times[, time_H_L := do.call(pmin, c(.SD, na.rm = TRUE)), .SDcols = patterns("time_H_L")]
    } else {
      # Simulate time to de novo cancer onset
      m_times[, time_H_P := query_distr("r", .N, d_time_H_P$distr, d_time_H_P$params)]
    }
  })
  return(NULL)
}


# Simulate development of additional lesions and output lesion-level data table
simulate_additional_lesions <- function(m_times, l_params_all) {
  l_df_lesion <- with(as.list(l_params_all), {
    # Initialize data structure to save lesion data
    l_df_lesion <- list()
    
    # By lesion type
    for (i in v_lesion) {
      # Subset to individuals with lesion onset before death
      m_lesion <- m_times[time_H_Do > get(paste0("time_H_L", i))]
      
      # Account for case with no lesions
      if (nrow(m_lesion) > 0) {
        # Set initial lesion ID and individual onset time variable
        m_lesion[, `:=` (lesion_id = 1, time_L_Lj = 0)]
        
        # Reset time_H_L variable to be lesion-specific
        m_lesion[, time_H_L := NULL]
        setnames(m_lesion, paste0("time_H_L", i), "time_H_L")
        
        # Create new row for each additional lesion (account for case where no one has >1 lesion)
        if (nrow(m_lesion[get(paste0("n_L", i)) > 0]) > 0) {
          m_other_lesion <- m_lesion[get(paste0("n_L", i)) > 0, .(
            time_H_Do,
            time_H_L = time_H_L,
            lesion_id = seq(2, get(paste0("n_L", i)) + 1)
          ), 
          pt_id]
          
          # Simulate onset time of each additional lesion
          m_other_lesion[, time_L_Lj := query_distr(
            "r", .N, 
            get(paste0("d_time_L", i, "_Lj"))$distr, 
            lapply(get(paste0("d_time_L", i, "_Lj"))$params, '*', (time_H_Do - time_H_L))
          )]
          
          # Append original lesion data
          m_lesion <- rbindlist(list(
            m_lesion[, .SD, .SDcols = intersect(names(m_lesion), names(m_other_lesion))],
            m_other_lesion
          ))
        }
        
        # Save lesion data frame
        setkey(m_lesion, pt_id)
        l_df_lesion[[i]] <- m_lesion
      }
    }
    
    # Return lesion data
    return(l_df_lesion)
  })
  
  return(l_df_lesion)
}

# Simulate lesion progression to cancer
simulate_lesion_progression <- function(m_times, l_df_lesion, l_params_all) {
  with(as.list(l_params_all), {
    # By lesion type
    for (i in v_lesion) {
      # Extract lesion data
      m_lesion <- l_df_lesion[[i]]
      
      # Account for case of no lesions
      if (!is.null(m_lesion)) {
        # Sample time to cancer
        var_convert <- paste0("time_L", i, "_P")
        m_lesion[, time_Lj_P := query_distr(
          "r", .N, 
          get(paste0("d_", var_convert))$distr, 
          get(paste0("d_", var_convert))$params)]
        
        # Calculate time from first lesion to lesion-specific cancer onset
        m_lesion[, time_L_P := time_L_Lj + time_Lj_P]
        
        # Calculate time from healthy to lesion-specific cancer onset (kept for summarization of lesion multiplicity by age)
        m_lesion[, time_H_P := time_H_L + time_L_P]
        
        # Get earliest time to cancer onset and rename to be lesion-type-specific
        m_times_cancer <- m_lesion[, .(time_L_P = min(time_L_P)), by = pt_id]
        
        # Merge to patient data
        m_times[m_times_cancer, (var_convert) := i.time_L_P]
        
        # Calculate time from healthy to lesion-type-specific cancer onset
        m_times[, paste0("time_H_P_", i) := get(paste0("time_H_L", i)) + get(var_convert)]
        
        # Reset NA to Inf since max.col cannot handle NA
        m_times[is.na(get(paste0("time_H_P_", i))), paste0("time_H_P_", i) := Inf]
      } else {
        m_times[, paste0("time_H_P_", i) := Inf]
      }
    }
    
    # Take earliest cancer onset time as time_H_P
    m_times[, time_H_P := do.call(pmin, c(.SD, na.rm = TRUE)),
            .SDcols = paste0("time_H_P_", v_lesion)]
    
    # Get source lesion for preclinical cancer
    m_times[is.finite(time_H_P),
            source_lesion := l_params_all$v_lesion[max.col(-1*.SD)], 
            .SDcols = paste0("time_H_P_", v_lesion)]
    
    # Create censor variable for assessing lesion characteristics and merge to lesion data
    m_times[, time_lesion_censor := pmin(time_H_Do, time_H_P, na.rm = TRUE)]
    for (i in names(l_df_lesion)) {
      l_df_lesion[[i]][m_times, time_lesion_censor := i.time_lesion_censor]
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
    l_m_lesion = NULL,
    verbose = FALSE
) {
  
  #### Initialize screening variables ####
  
  # Initialize number of screening tests and false positives
  m_times[, `:=` (ct_screen = 0,
                  ct_FP = 0)]
  
  # Set first screening age and screening end age
  m_times[time_screen_censor > age_screen_start, `:=` (
    screen_age = age_screen_start,
    time_screen_end = age_screen_end)]
  
  
  #### Case 1: Simulate screening before disease  ####
  
  # Set variable for time to disease onset
  if ('L' %in% l_params_all$v_states) {
    var_onset <- "time_H_L"
  } else {
    var_onset <- "time_H_P"
  }
  
  # Get routine screening ages
  v_screen_ages <- seq(age_screen_start, age_screen_end, int_screen)
  
  # Get number of tests before earliest of death or disease onset
  m_times[, ct_screen := findInterval(pmin(time_H_D, get(var_onset)), v_screen_ages, left.open=T)]
  
  # Simulate number of false positives
  m_times[, `:=` (
    ct_FP := rbinom(
      .N,
      size = ct_screen,
      prob = 1 - p_spec
    ))]
  
  # Get first screening age with lesion present and set to NA if after censor date
  m_times[, screen_age := v_screen_ages[ct_screen + 1], by = ct_screen]
  m_times[time_screen_censor <= screen_age, screen_age := NA]
  
  ##### Case 2: Lesions developed by screen age, but no preclinical cancer ####
  if ('L' %in% l_params_all$v_states) {
    # Prepare data for screening
    v_n_lesion_screen <- c() # For number of lesions eligible for screen detection
    
    for (i in names(l_m_lesion)) {
      # Get lesion data
      m_lesion <- l_m_lesion[[i]]
      
      # Merge screen age, as well as minimum of cancer onset and death as lesion 
      # screening censor date (to be updated after every screening)
      m_lesion[m_times, `:=` (time_screen_censor = pmin(i.time_H_P, time_H_Do, na.rm=T),
                              time_screen_end = i.time_screen_end,
                              screen_age = i.screen_age)]
      
      # Initialize number of tests among patients who are not censored and flag whether lesion was removed
      m_lesion[!is.na(screen_age), `:=` (ct_screen = 0,
                                         ct_FP = 0,
                                         fl_removed = 0)]
      
      # Count number of lesions available for screen detection
      v_n_lesion_screen <- c(v_n_lesion_screen, m_lesion[!is.na(screen_age), .N])
    }
    
    # Get number of lesions available for screening across all lesion types
    n_lesion_screen <- sum(v_n_lesion_screen)
    
    # Loop through screening tests until there are no more lesions to screen
    if (verbose) print("Number of lesions remaining to screen:")
    while (n_lesion_screen > 0) {
      if (verbose) print(n_lesion_screen)
      
      # Reset vector to store counts for number of screening-eligible lesions
      v_n_lesion_screen <- c()
      
      ###### 2.1 Flag whether eligible lesions would produce positive screen
      for (i in names(l_m_lesion)) {
        # Get lesion data
        m_lesion <- l_m_lesion[[i]]
        
        # Flag lesions that are eligible to be detected at current screen time
        m_lesion[fl_removed == 0 & time_H_L + time_L_Lj <= screen_age & time_screen_censor > screen_age,
                 fl_eligible := 1]
        
        # Sample whether lesion produces positive screen based on sensitivity for lesion size at diagnosis
        # and increment screen test count
        m_lesion[fl_eligible == 1, `:=` (
          fl_positive = rbinom(
            .N,
            size = 1,
            prob = p_sens
          ))]
      }
      
      # Check positive screens for any lesion for possible downstream diagnostic test detection
      m_detect <- rbindlist(l_m_lesion)[
        fl_eligible == 1,
        .(fl_positive_any = max(fl_positive)),
        by = pt_id]
      
      # Flag detection with diagnostic test
      for (i in names(l_m_lesion)) {
        # Get lesion data
        m_lesion <- l_m_lesion[[i]]
        
        # Merge positive screening test flag
        m_lesion[m_detect, fl_positive_any := i.fl_positive_any]
        
        # Sample detection with diagnostic test follow up
        m_lesion[fl_eligible == 1 & fl_positive_any == 1,`:=` (
          fl_detected = rbinom(
            .N,
            size = 1,
            prob = p_sens[['COL']][[i]]
          ))]
      }
      
      ####### 2.1.1 Detected lesions
      # To determine detection and follow up surveillance - combine lesion datasets, 
      # save screening age
      m_detected <- rbindlist(l_m_lesion)[
        fl_detected == 1,
        .(screen_age = unique(screen_age)),
        by = pt_id]
      
      setkey(m_detected, pt_id)
      
      # Increment test count
      m_times[m_detected, `:=` (
        ct_screen = ct_screen + 1,
        ct_diag = ct_diag + 1)]
      
      # Update surveillance data from detected lesions
      m_times[m_detected, `:=` (prev_screen_age = i.screen_age)]
      
      
      ####### 2.1.2 Undetected lesions
      # Get patients with no detected eligible lesions
      m_undetected <- rbindlist(l_m_lesion)[
        fl_eligible == 1 & !(pt_id %in% m_detected$pt_id),
        .(screen_age = unique(screen_age),
          fl_positive_any = max(fl_positive_any, 0, na.rm=T)), by = pt_id]
      setkey(m_undetected, pt_id)
      
      # Increment test count
      m_times[m_undetected, `:=` (
        ct_screen = ct_screen + 1,
        ct_diag = ct_diag + (fl_positive_any == 1))]
      
      # For patients with no detected eligible lesions, set surveillance risk to 0
      m_times[m_undetected, 
              `:=` (prev_screen_age = i.screen_age)]
      
      ###### 2.2 Flag patients with no active lesions
      # Find patients with no eligible lesions at screening age but who have not transitioned to preclinical cancer
      m_inactive <- rbindlist(l_m_lesion)[!is.na(screen_age), .(
        screen_age = unique(screen_age),
        fl_eligible_any = max(fl_eligible, 0, na.rm=T)), by = pt_id][fl_eligible_any == 0]
      
      # Sample false positives among people without any active lesions
      m_inactive[, fl_FP := rbinom(
        .N,
        size = 1,
        prob = 1 - p_spec[[modality]])]
      
      # Increment test count
      m_times[m_inactive[fl_FP==1], `:=` (
        ct_screen = ct_screen + 1,
        ct_diag = ct_diag + 1)]
      
      m_times[m_inactive[fl_FP==0], `:=` (
        ct_screen = ct_screen + 1)]
      
      # Update screen age and reset false positive flag
      m_times[m_inactive, `:=` (prev_screen_age = i.screen_age)]
      
      
      ###### 2.3 Update screening info and outcomes
      for (i in names(l_m_lesion)) {
        # Get lesion data
        m_lesion <- l_m_lesion[[i]]
        
        # Reset time to cancer onset for detected lesions to Inf and set detection time
        m_lesion[fl_detected == 1, `:=` (fl_removed = 1,
                                         time_H_P = Inf,
                                         time_detected = screen_age)]
        
        # Get earliest time to cancer onset and rename to be lesion-type-specific
        # Note: could potentially optimize by only updating for patients with lesions detected
        m_times_cancer <- m_lesion[, .(time_H_P = min(time_H_P)), by = pt_id]
        
        # Merge to patient data
        m_times[m_times_cancer, paste0("time_H_P_", i) := i.time_H_P]
      }
      
      # Update progression to cancer for people who had lesions detected
      # NOTE: save_orig = FALSE since columns were saved at the beginning, with originally time_H_P which is already changed here
      run_lesion_counterfactual(m_times, l_params_all, v_lesion_subset = names(l_m_lesion), save_orig = FALSE)
      
      # Update screening info for next round
      for (i in names(l_m_lesion)) {
        # Get lesion data
        m_lesion <- l_m_lesion[[i]]
        
        # Merge updates for any with detected lesions and set screen type to surveillance
        m_lesion[m_detected, `:=` (int_test_next = i.int_test_next,
                                   surveil_risk = i.surveil_risk_next,
                                   screen_type = "surveillance",
                                   time_screen_end = age_surveil_end)]
        
        # Merge updates for any with eligible lesions but none detected (screen type stays the same)
        m_lesion[m_undetected, `:=` (int_test_next = i.int_test_next,
                                     surveil_risk = 0)]
        
        # Merge updates for any with no eligible lesions
        m_lesion[m_inactive, `:=` (int_test_next = i.int_test_next,
                                   surveil_risk = 0)]
        
        # Merge new minimum of cancer onset and death as lesion screening censor date
        m_lesion[m_times, `:=` (time_screen_censor = pmin(i.time_H_P, time_H_Do, na.rm=T))]
        
        # Update age for next screening and set screen age to NA for patients who have been censored or have no lesions left
        m_lesion[, screen_age := screen_age + int_test_next]
        m_lesion[screen_age >= time_screen_censor | screen_age > time_screen_end, screen_age := NA]
        
        # Reset helper variables for updated screening test
        m_lesion[, `:=` (fl_eligible = NULL,
                         fl_positive = NULL,
                         fl_positive_any = NULL,
                         fl_detected = NULL,
                         int_test_next = NULL)]
        
        # Count number of lesions available for screening
        v_n_lesion_screen <- c(v_n_lesion_screen, m_lesion[!is.na(screen_age), .N])
      }
      
      # Get number of lesions available for screening across all lesion types
      n_lesion_screen <- sum(v_n_lesion_screen)
    }
  }
  
  ##### Case 3: Preclinical cancer screening and surveillance
  # Note: screen censor and mortality outcomes already updated using run_lesion_counterfactual()
  
  # Update screen age for people who received screenings at the lesion stage
  m_times[!is.na(prev_screen_age), `:=` (screen_age = prev_screen_age + int_test_next)]
  m_times[screen_age >= time_screen_censor | screen_age > time_screen_end, screen_age := NA]
  
  if (verbose) print("Number of individuals with preclinical cancer remaining to screen:")
  while (m_times[!is.na(screen_age), .N] > 0) {
    if (verbose) print(m_times[!is.na(screen_age), .N])
    
    # Sample whether cancer leads to positive screen
    # If colonoscopy, positive screen --> detection and transition to clinical cancer
    if (modality == 'COL') {
      m_times[time_H_P < screen_age, fl_detected := rbinom(
        .N,
        size = 1,
        prob = p_sens[[modality]][["cancer"]]
      )]
    } else {
      # For routine screening, use stool test
      m_times[time_H_P < screen_age & screen_type == "routine", fl_positive := rbinom(
        .N,
        size = 1,
        prob = p_sens[[modality]][["cancer"]]
      )]
      
      # For routine screening with positive test or surveillance, sample detection with colonoscopy follow up
      m_times[time_H_P < screen_age & ((screen_type == "routine" & fl_positive == 1) | (screen_type == "surveillance")), 
              fl_detected := rbinom(
                .N,
                size = 1,
                prob = p_sens[['COL']][["cancer"]]
              )]
    }
    
    # Convert detected cases to clinical cancer, get stage at diagnosis, set next screen age to NA
    # Note: Assumes 4 stages
    m_times[!is.na(screen_age) & fl_detected==1, `:=` (
      time_H_C = screen_age,
      stage_dx = fcase(
        time_H_C - get(paste0("time_P", stage_dx, "_C")) <= screen_age, as.double(stage_dx),
        time_H_C - get(paste0("time_P", stage_dx, "_C")) -
          get(paste0("time_P", stage_dx-1, "_P", stage_dx)) <= screen_age, as.double(stage_dx-1),
        time_H_C - get(paste0("time_P", stage_dx, "_C")) -
          get(paste0("time_P", stage_dx-1, "_P", stage_dx)) -
          get(paste0("time_P", stage_dx-2, "_P", stage_dx-1)) <= screen_age, as.double(stage_dx-2),
        time_H_C - get(paste0("time_P", stage_dx, "_C")) -
          get(paste0("time_P", stage_dx-1, "_P", stage_dx)) -
          get(paste0("time_P", stage_dx-2, "_P", stage_dx-1)) -
          get(paste0("time_P", stage_dx-3, "_P", stage_dx-2)) <= screen_age, as.double(stage_dx-3)
      ),
      screen_age = NA),
      by = stage_dx]
    
    # If not detected, set surveillance risk to zero and calculate next screening interval
    m_times[time_H_P < screen_age & !(fl_detected %in% 1), `:=` (
      int_test_next = v_screen_int[modality])]
    
    # Update screen age
    m_times[!is.na(screen_age), screen_age := screen_age + int_test_next]
    m_times[screen_age >= time_screen_censor | screen_age > time_screen_end, screen_age := NA]
  }
  
  ##### For patients with detected cancer, resimulate progression to clinical cancer and mortality outcomes
  
  # Recalculate time to death from cancer among people with screen-detected cancer
  with(l_params_all, {
    m_times[fl_detected == 1, time_C_Dc := query_distr(
      "r", .N,
      d_time_C_Dc[[v_surv[stage_dx]]]$distr,
      d_time_C_Dc[[v_surv[stage_dx]]]$params,
      d_time_C_Dc[[v_surv[stage_dx]]]$src
    ), by = stage_dx]
  })
  
  # Calculate death from cancer - constrain to be at least the former time to death from cancer
  m_times[fl_detected == 1, time_H_Dc := pmax(time_H_Dc, time_H_C + time_C_Dc)]
  
  # Calculate mortality outcomes
  calc_mortality_outcomes(m_times)
  
  
  #### 5. Calculate outcomes ===========================================
  
  # Calculate life years in population
  l_outcomes <- list()
  l_outcomes$lifeyears <- m_times[, sum(time_H_D)] / m_times[, .N] * rate_unit_outcome
  
  # Add number of colonoscopies by summing colonoscopies and adding patients with symptom-detected cancer, and number of other screening tests if applicable
  if (modality == 'COL') {
    l_outcomes$n_screen_col <- (m_times[, sum(ct_screen)] + m_times[time_H_C < time_H_D & !(fl_detected %in% 1), .N]) / m_times[, .N] * rate_unit_outcome
    l_outcomes$n_screen_other <- 0
  } else {
    l_outcomes$n_screen_col <- (m_times[, sum(ct_diag)] + m_times[time_H_C < time_H_D & !(fl_detected %in% 1), .N]) / m_times[, .N] * rate_unit_outcome
    l_outcomes$n_screen_other <- m_times[, sum(ct_screen)] / m_times[, .N] * rate_unit_outcome
  }
  
  return(l_outcomes)
}
