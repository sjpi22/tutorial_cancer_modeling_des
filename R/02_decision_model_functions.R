#' Generate time to event data for decision model
#'
#' @param l_params_all List with all parameters of decision model
#' 
#' @return Data table of time to event data for each individual
#' 
#' @export
run_model <- function(l_params_all, verbose = FALSE) {
  ##############################################################################
  #### Initialization
  ##############################################################################
  # Set seed
  if(!is.null(l_params_all$seed)) {
    set.seed(l_params_all$seed)
  }
  
  # Initialize list of results for each strategy
  l_results <- list()
  ##############################################################################
  #### Generate data under base case strategy
  ##############################################################################
  # Initialize matrix of patient data
  m_cohort_base <- initialize_cohort(l_params_all, verbose = verbose)
  
  # Simulate baseline characteristics
  simulate_baseline_data(m_cohort_base, l_params_all, verbose = verbose)
  
  # Simulate disease (lesion and/or cancer) onset
  simulate_disease_onset(m_cohort_base, l_params_all, verbose = verbose)
    
  # Simulate cancer progression
  simulate_cancer_progression(m_cohort_base, l_params_all, verbose = verbose)
  
  # Simulate cancer mortality
  simulate_cancer_mortality(m_cohort_base, l_params_all, verbose = verbose)
  
  # Compile overall mortality outcomes from background and cancer data
  calc_mortality_outcomes(m_cohort_base, verbose = verbose)
  
  # Add to results
  l_results[[l_params_all$v_strats[1]]] <- m_cohort_base
  
  return(l_results)
}

# Initialize time to event matrix for patient cohort with patient IDs and baseline strategy
initialize_cohort <- function(l_params_all, verbose = FALSE) {
  if(verbose) print('Initializing cohort')
  m_times_init <- with(as.list(l_params_all), {
    m_times <- data.table(
      pt_id = 1:n_cohort,
      modality = v_strats[1]
    )
    setkey(m_times, pt_id)
    return(m_times)
  }
  )
  return(m_times_init)
}

# Simulate baseline characteristics
simulate_baseline_data <- function(m_times, l_params_all, verbose = FALSE) {
  if(verbose) print('Simulating baseline data')
  with(as.list(l_params_all), {
    # Sample male / female
    m_times[, b_male := query_distr("r", .N, b_male$distr, b_male$params)]
    
    # Sample time to death from other causes if male
    m_times[b_male == 1, time_H_Do := query_distr("r", .N, time_H_Do_male$distr, time_H_Do_male$params)]
    
    # Sample time to death from other causes if female
    m_times[b_male == 0, time_H_Do := query_distr("r", .N, time_H_Do_female$distr, time_H_Do_female$params)]
  }
  )
}


# Simulate disease (lesion or cancer) onset
simulate_disease_onset <- function(m_times, l_params_all, verbose = FALSE) {
  if(verbose) print('Simulating disease onset')
  with(as.list(l_params_all), {
    # If including a lesion state, simulate time to lesion and cancer onset
    if("L" %in% v_states) {
      # Simulate time to lesion onset
      m_times[, time_H_L := query_distr("r", .N, time_H_L$distr, time_H_L$params)]
      
      # Simulate time from lesion to cancer onset
      m_times[, time_L_P := query_distr("r", .N, time_L_P$distr, time_L_P$params)]
      
      # Calculate time from birth to cancer onset
      m_times[, time_H_P := time_H_L + time_L_P]
    } else {
      # Simulate time to cancer onset
      m_times[, time_H_P := query_distr("r", .N, time_H_P$distr, time_H_P$params)]
    }
  })
}

# Simulate cancer stage progression
simulate_cancer_progression <- function(m_times, l_params_all, verbose = FALSE) {
  if(verbose) print('Simulating cancer progression')
  with(l_params_all, {
    # Initialize stage of diagnosis variable
    m_times[, stage_dx := eval(parse(text = paste0("as.", class(v_cancer), "(NA)")))]
      
    # Initialize time from cancer onset to detection
    m_times[, time_P_C := 0]
    
    for (i in 1:length(v_cancer)) {
      # Time from preclinical cancer stage to symptomatic detection
      var_detect <- paste0("time_P", v_cancer[i], "_C")
      m_times[is.na(stage_dx), (var_detect) := query_distr("r", .N, get(var_detect)$distr, get(var_detect)$params)]
      
      if (i < length(v_cancer)) {
        # Time from one stage to next in preclinical cancer
        var_progress <- paste0("time_P", v_cancer[i], "_P", v_cancer[i+1])
        m_times[is.na(stage_dx), (var_progress) := query_distr("r", .N, get(var_progress)$distr, get(var_progress)$params)]
        
        # Stage at cancer diagnosis - earliest stage at which time to detection is less than time to next stage
        m_times[is.na(stage_dx), stage_dx := fifelse(get(var_detect) < get(var_progress), v_cancer[i], NA)]
        
        # Add time to detection to running total of time from cancer onset to detection
        m_times[is.na(stage_dx) | (stage_dx == v_cancer[i]), time_P_C := fifelse(get(var_detect) < get(var_progress), 
                                                                                 time_P_C + get(var_detect), 
                                                                                 time_P_C + get(var_progress))]
      } else {
        # If not detected yet, set last stage at cancer diagnosis
        m_times[is.na(stage_dx), stage_dx := v_cancer[i]]
        
        # Add time to detection to running total of time from cancer onset to detection
        m_times[stage_dx == v_cancer[i], time_P_C := time_P_C + get(var_detect)]
      }
    }
    
  }
  )
  
  # Time to cancer diagnosis
  m_times[, time_H_C := time_H_P + time_P_C]
}

# Simulate cancer mortality
simulate_cancer_mortality <- function(m_times, l_params_all, verbose = FALSE) {
  if(verbose) print('Simulating cancer progression')
  with(l_params_all, {
    # Cancer mortality after diagnosis
    for (i in 1:length(v_cancer)) {
      var <- paste0("time_C", v_cancer[i], "_Dc")
      m_times[stage_dx == v_cancer[i], time_C_Dc := query_distr("r", .N, get(var)$distr, get(var)$params)]
    }
  }
  )
  
  # Calculate death from cancer
  m_times[, time_H_Dc := time_H_C + time_C_Dc]
}

# Generate mortality outcomes
calc_mortality_outcomes <- function(m_times, verbose = FALSE) {
  if(verbose) print('Calculating mortality outcomes')
  
  # Join cancer data to patient-level data
  # Calculate death all-cause death and cause of death
  m_times[, time_H_D := pmin(time_H_Do, time_H_Dc, na.rm = TRUE)]
  m_times[, fl_Dc := (time_H_Do > pmin(time_H_Dc, Inf, na.rm = TRUE))]
  
  # Calculate survival from cancer diagnosis
  m_times[time_H_C <= time_H_D, time_C_D := time_H_D - time_H_C]
}
