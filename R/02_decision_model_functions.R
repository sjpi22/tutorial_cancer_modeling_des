#' Generate time to event data for decision model
#'
#' @param l_params_all List with all parameters of decision model
#' 
#' @return Data table of time to event data for each individual
#' 
#' @export
run_model <- function(l_params_all) {
  ##############################################################################
  #### Initialization
  ##############################################################################
  # Set seed
  if(!is.null(l_params_all$seed)) {
    set.seed(l_params_all$seed)
  }
  
  # Initialize list of results for each strategy
  l_results <- list()
  
  # Convert dataframe of variables to list of lists for easier use
  l_distr <- lapply(split(l_params_all$df_vars, l_params_all$df_vars$varID), function(x) as.list(x))
  
  # Clean up params lists
  l_distr <- lapply(l_distr, function(x) {
    x$params <- x$params[[1]]
    return(x)
  }
  )
  
  # Add list version of variables
  l_params_all$l_distr <- l_distr
  
  ##############################################################################
  #### Generate data under base case strategy
  ##############################################################################
  # Initialize matrix of patient data
  m_cohort_init <- initialize_cohort(l_params_all)
  
  # Generate baseline characteristics (demographics, risk factors, and time to death from other causes)
  m_cohort_init <- simulate_pt_data(m_cohort_init, l_params_all, "baseline")
  
  # Generate precancerous lesions
  m_lesions <- simulate_lesion_data(m_cohort_init, l_params_all)
  
  # Consolidate lesion-level data to patient-level data
  m_cohort_base <- simulate_cancer_mort(m_cohort_init, m_cohort_lesions, l_params_all)
  
  # Generate overall survival
  l_results[[l_params_all$v_states[1]]] <- list(m_cohort = m_cohort_base, 
                                                m_lesions = m_lesions)
  
  ##############################################################################
  #### Regenerate disease process under other intervention scenarios
  ##############################################################################
  if(length(l_params_all$v_strats) > 1) {
    for (str in l_params_all$v_strats[-1]) {
    }
  }
  
  # Set labels for results list
  names(l_results) <- l_params_all$v_strats
  
  return(l_results)
}

# Initialize time to event matrix for patient cohort with patient IDs and baseline strategy
initialize_cohort <- function(l_params_all) {
  m_times_init <- with(as.list(l_params_all), {
    m_times <- data.table(
      pt_id = 1:n_cohort,
      modality = v_strats[1]
    )
    return(m_times)
  }
  )
  return(m_times_init)
}

# Generate patient-level data
simulate_pt_data <- function(m_times, l_params_all, vargroups = NULL) {
  m_times_updated <- with(as.list(l_params_all), {
    # If vargroups is populated, filter to variables
    if(!is.null(vargroups)) df_vars_calc <- df_vars[df_vars$vargroup %in% vargroups, ]
    else df_vars_calc <- df_vars
    
    # Sample variables from corresponding distributions
    m_times[, (df_vars_calc$varname) := lapply(
      # Convert column names to corresponding distribution variables
      l_distr[df_vars_calc$varID],
      
      # Apply random sampling from distribution
      function(distr) query_distr("r", nrow(m_times), distr$distr, distr$params))
    ]
      # m_times[, unlist(lapply(.SD, my.summary), recursive = FALSE), 
         # .SDcols = ColChoice, by = category]
    return(m_times)
  }
  )
  return(m_times_updated)
}


# Generate precancerous lesions
simulate_lesion_data <- function(m_times, l_params_all, vargroup_term = "lesion") {
  # Initialize lesion dataframe
  m_times_lesion <- data.table()
  
  # Create lesion-level dataframe for each lesion type
  for(lesiontype in l_params_all$v_lesions) {
    # Simulate patient level lesion variables (time to first lesion and number of lesions)
    temp_m_times <- simulate_pt_data(m_times, l_params_all, vargroup = paste(vargroup_term, "pt-level", lesiontype)) %>%
      # Filter to lesions that occur before death time for efficiency
      filter(time_0_1 < time_0_Do)
    
    # Create row for every lesion
    temp_m_times_lesion <- temp_m_times %>%
      group_by(pt_id) %>%
      uncount(n_add + 1, .remove = FALSE) %>%
      mutate(lesiontype = lesiontype,
             lesionID = 1:n()) %>%
      setDT()
    
    # Get lesion-specific onset time and cancer conversion time
    temp_m_times_lesion <- simulate_pt_data(temp_m_times_lesion, l_params_all, vargroup = paste(vargroup_term, "lesion-level", lesiontype))
    
    # Reset time of first lesion
    temp_m_times_lesion[lesionID == 1, time_1_1i := 0]
    
    # Calculate time from birth to lesion-specific onset 
    temp_m_times_lesion[, time_0_1i := time_0_1 + time_1_1i]
    
    # Filter to lesions that occur before death time for efficiency
    temp_m_times_lesion <- temp_m_times_lesion[time_0_1i < time_0_Do, ] %>%
      # Reset lesion ID in order of lesion onset
      group_by(pt_id) %>%
      arrange(time_0_1i, .by_group = TRUE) %>%
      mutate(lesionID = 1:n()) %>%
      ungroup() %>%
      setDT()
    
    # Calculate time from birth to lesion cancer conversion
    temp_m_times_lesion[, time_0_2i := time_0_1i + time_1i_2]
    
    # Bind to lesion-level data table
    m_times_lesion <- rbind(m_times_lesion, temp_m_times_lesion)
  }
  
  return(m_times_lesion)
}

# Simulate cancer mortality
simulate_cancer_mort <- function(m_times, m_times_lesion, l_params_all) {
  # Find lesion that first converts to cancer
  m_times_cancer <- m_times_lesion %>%
    group_by(pt_id) %>%
    mutate(time_0_2 = min(time_0_2i, na.rm = TRUE)) %>%
    ungroup %>%
    filter(time_0_2 == time_0_2i) %>%
    select(-time_0_2i) %>%
    mutate(time_1_2i = time_0_2 - time_0_1) %>%
    setDT()
  
  # Simulate cancer-related survival
  m_times_cancer <- simulate_pt_data(m_times_cancer, l_params_all, vargroup = "cancer")
  m_times_cancer[, time_2_3 := time_2_Du * prog_3s]
  m_times_cancer[, time_0_3 := time_0_2 + time_2_3]
  
  # Join lesion data to patient-level data
  m_times_updated <- m_times %>%
    left_join(m_times_cancer, by = "pt_id", suffix=c("",".y")) %>%
              select(-ends_with(".y"))
}

# Generate overall survival

