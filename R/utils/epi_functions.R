

# Calculate prevalence in the intervals spanned by [v_ages[i], v_ages[i+1])
calc_prevalence <- function(m_time, start_var, end_var, censor_var, v_ages = NULL) {
  # Number in cohort
  n_cohort <- nrow(m_time)
  
  # If v_ages is null, set to find prevalence by year
  if(is.null(v_ages)) {
    v_ages <- seq(0, ceiling(max(m_time[, get(censor_var)])))
  } else {
    # Sort v_ages
    v_ages <- sort(v_ages)
  }
  
  # Create age range data frame
  res <- data.frame(
    age_start = v_ages[-length(v_ages)],
    age_end = v_ages[-1]
  )
  
  # Calculate prevalence in the given age ranges
  res <- cbind(res, t(mapply(
    function(age_start, age_end) {
      denom = sum(pmax(m_time[ , pmin(get(censor_var), age_end)] - age_start, 0))
      num = sum(pmax(m_time[ , pmin(get(end_var), get(censor_var), age_end)] - m_time[ , pmax(get(start_var), age_start)], 0, na.rm = TRUE))
      return(c(person_years_cases = num, person_years_total = denom, prevalence = num/denom))
    }, res$age_start, res$age_end))) %>%
    mutate(n_total = round(person_years_total / (age_end - age_start))) %>%
    relocate(prevalence, .after = last_col())
  
  return(res)
}

# Calculate precancerous lesion (PCL) prevalence from lesion-level and individual-level data
calc_pcl_prevalence <- function(l_params_all, m_patients, m_lesions, v_ages) {
  # Initialize containers for results
  output_prev <- c()
  names_output_prev <- c()
  
  # Loop over lesion types
  for (lesiontype in l_params_all$v_lesions) {
    # Get times for specific lesion type
    temp_lesion_type <- unique(m_lesions[lesion_type == lesiontype, c("pt_id", "time_0_1")], by = c("pt_id", "time_0_1"))
    
    # Rename variable to avoid conflict when merging with original data
    setnames(temp_lesion_type, "time_0_1", "time_0_1_lesion")
    
    # Merge to final data
    temp_m_cohort <- merge(m_patients, temp_lesion_type, by = "pt_id", all.x = TRUE)
    temp_prev <- calc_prevalence(temp_m_cohort, "time_0_1_lesion", "time_0_2", "time_screen_censor", v_ages[[lesiontype]]) %>%
      select(-c("person_years_cases", "person_years_total")) %>%
      mutate(confint_lb = qbinom((1-l_params_all$conf_level)/2, size = n_total, prob = prevalence) / n_total,
             confint_ub = qbinom((1+l_params_all$conf_level)/2, size = n_total, prob = prevalence) / n_total)
    
    # Save outputs
    output_prev <- c(output_prev, list(temp_prev))
    names_output_prev <- c(names_output_prev, lesiontype)
  }
  
  # Save results
  names(output_prev) <- names_output_prev
  
  # Return prevalence output
  return(output_prev)
}

# Calculate annual incidence
calc_incidence <- function(m_time, time_var, censor_var, age_lower_bounds,
                           strat_var = NULL,
                           rate_unit = 100000) {
  
  # Number in cohort
  n_cohort <- nrow(m_time)
  
  # Create age category labels
  age_lower_bounds_orig <- age_lower_bounds # Save original vector of ages
  age_lower_bounds <- unique(c(0, age_lower_bounds)) # Add 0 as lower bound if necessary
  age_ranges <- paste(age_lower_bounds[-length(age_lower_bounds)], 
                      age_lower_bounds[-1], sep = "-")
  age_ranges <- c(age_ranges, paste0(max(age_lower_bounds), "+"))
  age_df <- data.frame(age_range = age_ranges, 
                       age_start = age_lower_bounds,
                       age_end = lead(age_lower_bounds, default = ceiling(max(m_time[, get(censor_var)]))))
  
  # Exposure time by age range
  person_years_at_risk <- data.frame(age_df, 
                                     total_atrisk = mapply(
                                       function(age_start, age_end) {
                                         total_atrisk = sum(pmax(m_time[ , pmin(get(censor_var), age_end)] - age_start, 0))
                                       }, age_df$age_start, age_df$age_end),
                                     age_diff = age_df$age_end - age_df$age_start) %>%
    mutate(n_population = round(total_atrisk / age_diff)) %>%
    filter(age_range %in% age_df$age_range[age_df$age_start %in% age_lower_bounds_orig])
  
  # Set grouping variables
  group_vars <- "age_range"
  if(!is.null(strat_var)) group_vars <- c(group_vars, strat_var)
  
  # Set rate unit
  if(is.null(rate_unit)) rate_unit <- 1
  
  # Events by category and stage
  event_counts <- m_time %>%
    filter(get(time_var) <= get(censor_var)) %>%
    mutate(age_range = sapply(get(time_var),
                              function(x) age_ranges[which.max(age_lower_bounds[x >= age_lower_bounds])])) %>%
    group_by_at(group_vars) %>%
    summarise(n_events = n(),
              .groups = "drop")
  
  # Make sure all values of stratifying variable are represented
  if(!is.null(strat_var)) {
    # Get sorted values of stratifying variable
    strat_var_vals <- sort(unique(m_time[[strat_var]]))
    
    # Initialize dataframe with all combinations of age range and stratifying variables
    all_var_vals <- data.frame(
      age_range = rep(sort(unique(event_counts$age_range)), each = length(strat_var_vals))
    )
    
    # Add stratifying variables
    all_var_vals[, as.character(strat_var)] <- rep(strat_var_vals, length(unique(event_counts$age_range)))
    
    # Merge event counts
    event_counts <- all_var_vals %>%
      left_join(event_counts, by = group_vars)
  }
  
  if(nrow(event_counts) > 0) {
    event_counts <- person_years_at_risk %>%
      left_join(event_counts, by = "age_range") %>%
      mutate(n_events = replace_na(n_events, 0)) %>%
      mutate(incidence_raw = n_events / total_atrisk) %>%
      mutate(incidence = incidence_raw * rate_unit)
  } else {
    event_counts <- person_years_at_risk %>%
      mutate(n_events = 0,
             incidence = 0)
  }
  
  return(event_counts)
}

