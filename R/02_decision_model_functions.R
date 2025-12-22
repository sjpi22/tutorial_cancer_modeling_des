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

