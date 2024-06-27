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
  
  # Generate precancerous lesions
  simulate_lesion_data(m_cohort_base, l_params_all, verbose = verbose)
  
  # Consolidate lesion-level data to patient-level data for cancer onset and progression
  simulate_cancer_progression(m_cohort_base, l_params_all, verbose = verbose)
  
  # Compile overall mortality outcomes from background and cancer data
  calc_mortality_outcomes(m_cohort_base, verbose = verbose)
  
  # Add to results
  l_results[[l_params_all$v_strats[1]]] <- m_cohort_base
  
  ##############################################################################
  #### Regenerate disease process under other intervention scenarios
  ##############################################################################
  if(length(l_params_all$v_strats) > 1) {
    for (str in l_params_all$v_strats[-1]) {
      
      # Add to results
      l_results[[str]] <- list(m_cohort = m_cohort_base, 
                               m_lesions = m_lesions)
    }
  }
  
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
simulate_lesion_data <- function(m_times, l_params_all, verbose = FALSE) {
  if(verbose) print('Simulating lesion data')
  with(as.list(l_params_all), {
    
    # Simulate time to lesion onset
    m_times[, time_0_1 := query_distr("r", .N, time_0_1$distr, time_0_1$params)]
    
    # Simulate time from lesion onset to cancer onset
    m_times[, time_1_2 := query_distr("r", .N, time_1_2$distr, time_1_2$params)]
    
    # Calculate time from birth to cancer onset
    m_times[, time_0_2 := time_0_1 + time_1_2]
    
  })
}

# Simulate cancer stage progression and mortality
simulate_cancer_progression <- function(m_times, l_params_all, verbose = FALSE) {
  if(verbose) print('Simulating cancer progression')
  with(l_params_all, {
    # Time from first to second stage
    m_times[, time_2i_2ii := query_distr("r", .N, time_2i_2ii$distr, time_2i_2ii$params)]
    
    # Time to symptomatic detection in first stage
    m_times[, time_2i_3 := query_distr("r", .N, time_2i_3$distr, time_2i_3$params)]
    
    # Time to symptomatic detection in second stage
    m_times[, time_2ii_3 := query_distr("r", .N, time_2ii_3$distr, time_2ii_3$params)]
    
    # Stage at cancer diagnosis - earliest stage at which time to detection is less than time to next stage
    m_times[, stage_dx := fifelse(time_2i_3 < time_2i_2ii, 1, 2)]
    m_times[, time_2_3 := fifelse(time_2i_3 < time_2i_2ii, time_2i_3, time_2i_2ii + time_2ii_3)]
    
    # Time to cancer diagnosis
    m_times[, time_0_3 := time_0_2 + time_2_3]
    
    # Cancer mortality after diagnosis
    for (stg in 1:length(v_cancer)) {
      var <- paste0("time_3", v_cancer[stg], "_Dc")
      m_times[stage_dx == stg, time_3_Dc := query_distr("r", .N, get(var)$distr, get(var)$params)]
    }
  }
  )
  
  # Calculate death from cancer
  m_times[, time_0_Dc := time_0_3 + time_3_Dc]
}


# Generate mortality outcomes
calc_mortality_outcomes <- function(m_times, verbose = FALSE) {
  if(verbose) print('Calculating mortality outcomes')
  
  # Join cancer data to patient-level data
  # Calculate death all-cause death and cause of death
  m_times[, time_0_D := pmin(time_0_Do, time_0_Dc, na.rm = TRUE)]
  m_times[, fl_death_cancer := (time_0_Do > pmin(time_0_Dc, Inf, na.rm = TRUE))]
  
  # Calculate screening censor date
  m_times[, time_screen_censor := pmin(time_0_D, time_0_3, na.rm = TRUE)]
  
  # Calculate survival from cancer diagnosis
  m_times[time_0_3 <= time_0_D, time_3_D := time_0_D - time_0_3]
}
