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
  
  # Generate baseline characteristics
  simulate_baseline_data(m_cohort_base, l_params_all, verbose = verbose)
  
  # Generate cancer onset
  simulate_cancer_onset(m_cohort_base, l_params_all, verbose = verbose)
  
  # Cancer progression
  simulate_cancer_progression(m_cohort_base, l_params_all, verbose = verbose)
  
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

# Generate baseline data - updates by reference in-place, so does not return anything
simulate_baseline_data <- function(m_times, l_params_all, verbose = FALSE) {
  if(verbose) print('Simulating baseline data')
  with(as.list(l_params_all), {
    # Sample male / female
    m_times[, b_male := query_distr("r", .N, b_male$distr, b_male$params)]
    
    # Sample time to death from other causes if male
    m_times[b_male == 1, time_0_Do := query_distr("r", .N, time_0_Do_male$distr, time_0_Do_male$params)]
    
    # Sample time to death from other causes if female
    m_times[b_male == 0, time_0_Do := query_distr("r", .N, time_0_Do_female$distr, time_0_Do_female$params)]
  }
  )
}


# Generate precancerous lesions - updates by reference in-place, so does not return anything
simulate_cancer_onset <- function(m_times, l_params_all, verbose = FALSE) {
  if(verbose) print('Simulating lesion data')
  with(as.list(l_params_all), {
    
    # Simulate time to cancer onset
    m_times[, time_0_1 := query_distr("r", .N, time_0_1$distr, time_0_1$params)]
  })
}

# Simulate cancer stage progression and mortality
simulate_cancer_progression <- function(m_times, l_params_all, verbose = FALSE) {
  if(verbose) print('Simulating cancer progression')
  with(l_params_all, {
    # Time from first to second stage
    m_times[, time_1i_1ii := query_distr("r", .N, time_1i_1ii$distr, time_1i_1ii$params)]
    
    # Time to symptomatic detection in first stage
    m_times[, time_1i_2 := query_distr("r", .N, time_1i_2$distr, time_1i_2$params)]
    
    # Time to symptomatic detection in second stage
    m_times[, time_1ii_2 := query_distr("r", .N, time_1ii_2$distr, time_1ii_2$params)]
    
    # Stage at cancer diagnosis - earliest stage at which time to detection is less than time to next stage
    m_times[, stage_dx := fifelse(time_1i_2 < time_1i_1ii, 1, 2)]
    m_times[, time_1_2 := fifelse(time_1i_2 < time_1i_1ii, time_1i_2, time_1i_1ii + time_1ii_2)]
    
    # Time to cancer diagnosis
    m_times[, time_0_2 := time_0_1 + time_1_2]
    
    # Cancer mortality after diagnosis
    for (stg in 1:length(v_cancer)) {
      var <- paste0("time_2", v_cancer[stg], "_Dc")
      m_times[stage_dx == stg, time_2_Dc := query_distr("r", .N, get(var)$distr, get(var)$params)]
    }
  }
  )
  
  # Calculate death from cancer
  m_times[, time_0_Dc := time_0_2 + time_2_Dc]
}


# Generate mortality outcomes
calc_mortality_outcomes <- function(m_times, verbose = FALSE) {
  if(verbose) print('Calculating mortality outcomes')
  
  # Join cancer data to patient-level data
  # Calculate death all-cause death and cause of death
  m_times[, time_0_D := pmin(time_0_Do, time_0_Dc, na.rm = TRUE)]
  m_times[, fl_death_cancer := (time_0_Do > pmin(time_0_Dc, Inf, na.rm = TRUE))]
  
  # Calculate survival from cancer diagnosis
  m_times[time_0_2 <= time_0_D, time_2_D := time_0_D - time_0_2]
}
