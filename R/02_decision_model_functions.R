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
  m_cohort_init <- initialize_cohort(l_params_all, verbose = verbose)
  
  # Generate baseline characteristics
  simulate_baseline_data(m_cohort_init, l_params_all, verbose = verbose)
  
  # Generate precancerous lesions
  m_lesions <- simulate_lesion_data(m_cohort_init, l_params_all, verbose = verbose)
  
  # Consolidate lesion-level data to patient-level data for cancer onset and progression
  m_cohort_cancer <- simulate_cancer_progression(m_cohort_init, m_lesions, l_params_all, 
                                                 verbose = verbose)
  
  # Compile overall mortality outcomes from backgrount and cancer data
  m_cohort_base <- calc_mortality_outcomes(m_cohort_init, m_cohort_cancer, 
                                           verbose = verbose)
  
  # Add to results
  l_results[[l_params_all$v_strats[1]]] <- list(m_cohort = m_cohort_base, 
                                                m_lesions = m_lesions)
  
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


# Generate precancerous lesions
simulate_lesion_data <- function(m_times, l_params_all, verbose = FALSE) {
  if(verbose) print('Simulating lesion data')
  
  # Create lesion-level dataframe for each lesion type
  m_times_lesion <- with(as.list(l_params_all), {
    # Initialize lesion dataframe
    m_times_lesion <- data.table()
    
    for(lesiontype in v_lesions) {
      # Record lesion type and simulate time to first lesion 
      m_times[, `:=` (time_0_1 = query_distr("r", .N, get(paste0("time_0_1", lesiontype))$distr, get(paste0("time_0_1", lesiontype))$params),
                      lesion_type = lesiontype)]
      
      # Simulate number of lesions before death if first lesion occurs before death from other causes
      m_times[time_0_1 < time_0_Do, n_add := query_distr("r", .N, get(paste0("n_add", lesiontype))$distr, lapply(get(paste0("n_add", lesiontype))$params, function(x) x*(time_0_Do - time_0_1)))]
      
      # Create row for every lesion
      temp_m_times_lesion <- m_times[time_0_1 < time_0_Do, .(lesionID = seq(n_add + 1)), by = names(m_times)] 
      
      # Get lesion-specific onset time
      temp_m_times_lesion[, time_1_1j := query_distr("r", .N, get(paste0("time_1_1j", lesiontype))$distr, eval(parse(text = get(paste0("time_1_1j", lesiontype))$params)))]
      
      # Reset time of first lesion
      temp_m_times_lesion[lesionID == 1, time_1_1j := 0]
      
      # Calculate time from birth to lesion-specific onset 
      temp_m_times_lesion[, time_0_1j := time_0_1 + time_1_1j]
      
      # Filter to lesions that occur before death time for efficiency
      temp_m_times_lesion <- temp_m_times_lesion[time_0_1j < time_0_Do] 
      
      # Order by patient ID and time of lesion onset and reset lesion IDs
      setorder(temp_m_times_lesion, pt_id, time_0_1j)
      temp_m_times_lesion[, lesionID := seq(.N),  by = pt_id]
      
      # Get lesion-specific onset time and cancer conversion time
      temp_m_times_lesion[, time_1j_2i := query_distr("r", .N, get(paste0("time_1j_2i", lesiontype))$distr, get(paste0("time_1j_2i", lesiontype))$params)]
      
      # Calculate time from birth to lesion cancer conversion
      temp_m_times_lesion[, time_0_2i := time_0_1j + time_1j_2i]
      
      # Bind to lesion-level data table
      m_times_lesion <- rbind(m_times_lesion, temp_m_times_lesion)
    }
    
    return(m_times_lesion)
  })
  

  # Delete new rows from m_times
  m_times[, c("time_0_1", "n_add", "lesion_type"):=NULL]
  
  return(m_times_lesion)
}

# Simulate cancer stage progression and mortality
simulate_cancer_progression <- function(m_times, m_times_lesion, l_params_all, verbose = FALSE) {
  if(verbose) print('Simulating cancer progression')
  
  # Calculate cancer onset from lesions
  m_times_cancer <- calc_cancer_onset(m_times, m_times_lesion, l_params_all)
  
  with(l_params_all, {
    # Time to next stage of preclinical cancer (variable name = time_2{start stage}_2{next stage})
    m_times_cancer[, (vars_preclin_nextstage) := lapply(vars_preclin_nextstage, function(var) query_distr("r", .N, get(var)$distr, get(var)$params))]
    
    # Time to detection from preclinical stage (variable name = time_2{stage}_3)
    m_times_cancer[, (vars_preclin_detect) := lapply(vars_preclin_detect, function(var) query_distr("r", .N, get(var)$distr, get(var)$params))]
    
    # Stage at cancer diagnosis - earliest stage at which time to detection is less than time to next stage
    m_times_cancer[, stage_dx := 0]
    m_times_cancer[, time_2_3 := 0]
    for (i in 1:(length(v_cancer)-1)) # Iteratively set stage at diagnosis
      m_times_cancer[stage_dx == 0, stage_dx := fifelse(get(vars_preclin_detect[i]) < get(vars_preclin_nextstage[i]), i, 0)]
      m_times_cancer[, time_2_3 := time_2_3 + (stage_dx == 0) * get(vars_preclin_nextstage[i]) + (stage_dx == i) * get(vars_preclin_detect[i])]
    m_times_cancer[stage_dx == 0, stage_dx := length(v_cancer)] # If not set yet, detected at last stage
    
    # Time to cancer diagnosis
    m_times_cancer[, time_0_3 := time_0_2 + time_2_3]
    
    # Cancer mortality after diagnosis
    for (stg in sort(unique(m_times_cancer$stage_dx))) {
      var <- paste0("time_3", v_cancer[stg], "_Dc")
      m_times_cancer[stage_dx == stg, time_3_Dc := query_distr("r", .N, get(var)$distr, get(var)$params)]
    }
  }
  )
  
  # Calculate death from cancer
  m_times_cancer[, time_0_Dc := time_0_3 + time_3_Dc]
  
  return(m_times_cancer)
}

# Calculate cancer onset from lesions
calc_cancer_onset <- function(m_times, m_times_lesion, l_params_all) {
  # Find time to first lesion that converts to cancer and flag lesion
  m_times_lesion[, time_0_2 := min(time_0_2i, na.rm = TRUE), by = pt_id]
  m_times_lesion[, first_cancerous_lesion := (time_0_2 == time_0_2i)]
  
  # Filter to earliest lesion and calculate time from first lesion to cancer onset
  m_times_cancer <- m_times_lesion[first_cancerous_lesion == 1, ]
  m_times_cancer[, time_1_2i := time_0_2 - time_0_1]
  
  return(m_times_cancer)
}

# Generate mortality outcomes
calc_mortality_outcomes <- function(m_times, m_times_cancer, verbose = FALSE) {
  if(verbose) print('Calculating mortality outcomes')
  
  # Join cancer data to patient-level data
  dupe_names <- intersect(names(m_times), names(m_times_cancer))[-1]
  m_times_updated <- merge(m_times, m_times_cancer[, (dupe_names) := NULL], by = "pt_id", all.x = TRUE)
  
  # Calculate death all-cause death and cause of death
  m_times_updated[, time_0_D := pmin(time_0_Do, time_0_Dc, na.rm = TRUE)]
  m_times_updated[, fl_death_cancer := (time_0_Do > pmin(time_0_Dc, Inf, na.rm = TRUE))]
  
  # Calculate screening censor date
  m_times_updated[, time_screen_censor := pmin(time_0_D, time_0_3, na.rm = TRUE)]
  
  # Calculate survival from cancer diagnosis
  m_times_updated[time_0_3 <= time_0_D, time_3_D := time_0_D - time_0_3]
  assertthat::assert_that(nrow(m_times_updated[(fl_death_cancer == 1) & is.na(time_3_D), ]) == 0) # Check that all patients with death due to cancer have time from diagnosis to death populated
  
  return(m_times_updated)
}
