

# Calculate prevalence in the intervals spanned by [v_ages[i], v_ages[i+1])
calc_prevalence <- function(m_time, start_var, end_var, censor_var, v_ages = NULL, conf_level = 0.95) {
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
    relocate(prevalence, .after = last_col()) %>%
    mutate(ci_lb = qbinom((1-conf_level)/2, size = n_total, prob = prevalence) / n_total,
           ci_ub = qbinom((1+conf_level)/2, size = n_total, prob = prevalence) / n_total) %>%
    mutate_all(~replace(., is.na(.), 0))
  
  return(res)
}

# Calculate annual incidence
calc_incidence <- function(m_time, time_var, censor_var, v_ages,
                           strat_var = NULL,
                           rate_unit = 100000) {
  
  # Number in cohort
  n_cohort <- nrow(m_time)
  
  # Create age category labels
  age_bounds <- unique(c(0, v_ages)) # Add 0 as lower bound if necessary
  age_ranges <- paste(age_bounds[-length(age_bounds)], 
                      age_bounds[-1], sep = "-")
  age_df <- data.frame(age_range = age_ranges, 
                       age_start = age_bounds[-length(age_bounds)],
                       age_end = age_bounds[-1])
  
  # Exposure time by age range
  person_years_at_risk <- data.frame(age_df, 
                                     total_atrisk = mapply(
                                       function(age_start, age_end) {
                                         total_atrisk = sum(pmax(m_time[ , pmin(get(censor_var), age_end)] - age_start, 0))
                                       }, age_df$age_start, age_df$age_end),
                                     age_diff = age_df$age_end - age_df$age_start) %>%
    mutate(n_population = round(total_atrisk / age_diff)) %>%
    filter(age_range %in% age_df$age_range[age_df$age_start %in% v_ages[-length(v_ages)]])
  
  # Set grouping variables
  grouping_vars <- "age_range"
  if(!is.null(strat_var)) grouping_vars <- c(grouping_vars, strat_var)
  
  # Set rate unit
  if(is.null(rate_unit)) rate_unit <- 1
  
  # Number of events by category
  event_counts <- m_time %>%
    filter(get(time_var) <= pmin(get(censor_var), max(age_bounds))) %>%
    mutate(age_range = sapply(get(time_var),
                              function(x) age_ranges[which.max(age_bounds[x >= age_bounds])])) %>%
    group_by_at(grouping_vars) %>%
    summarise(n_events = n(),
              .groups = "drop")
  
  # Make sure all values of stratifying variables are represented
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
      left_join(event_counts, by = grouping_vars)
  }
  
  if(nrow(event_counts) > 0) {
    event_counts <- person_years_at_risk %>%
      left_join(event_counts, by = "age_range") %>%
      mutate(n_events = replace_na(n_events, 0)) %>%
      mutate(unit = rate_unit) %>%
      mutate(incidence = n_events / total_atrisk * rate_unit,
             se = sqrt(n_events / age_diff) / (total_atrisk / age_diff) * rate_unit)
  } else {
    event_counts <- person_years_at_risk %>%
      mutate(n_events = 0,
             unit = rate_unit,
             incidence = 0,
             se = 0)
  }
  
  return(event_counts)
}

# Calculate stage distribution from patient-level data
calc_stage_distr <- function(m_time, grouping_var, event_var, censor_var, groups_expected = c(1,2,3,4), conf_level = 0.95) {
  # Count patients diagnosed at each stage
  distr <- m_time %>%
    filter(get(event_var) <= get(censor_var)) %>%
    group_by_at(grouping_var) %>%
    summarize(count = n(), .groups = 'drop') 
  
  # Create dataframe of expected stages
  df_groups_expected <- setNames(data.frame(groups_expected), grouping_var)
  
  # If any stages not represented, convert to 0
  distr <- df_groups_expected %>%
    left_join(distr, by = grouping_var) %>%
    mutate_all(~replace(., is.na(.), 0)) 
  
  # If no patients diagnosed at any stage, create dataframe of zeros
  if (nrow(distr) == 0) {
    distr <- df_groups_expected %>%
      mutate(count = 0)
  }
  
  # Add percentages and confidence intervals
  distr <- distr %>%
    mutate(pct = count / sum(count)) %>%
    mutate(ci_lb = qbinom((1-conf_level)/2, size = sum(count), prob = pct) / sum(count),
           ci_ub = qbinom((1+conf_level)/2, size = sum(count), prob = pct) / sum(count)) %>%
    mutate_all(~replace(., is.na(.), 0)) 
  
  return(distr)
}

# Get prevalence age range for input into calc_prevalence from list of prevalence dataframes
get_age_range_from_list <- function(l_df_inputs) {
  l_v_ages <- list()
  for (lesiontype in names(l_df_inputs)) {
    df_inputs <- l_df_inputs[[lesiontype]]
    l_v_ages[[lesiontype]] <- get_age_range(df_inputs)
  }
  
  return(l_v_ages)
}

# Get prevalence age range for input into calc_prevalence from prevalence dataframe
get_age_range <- function(df_inputs) {
  v_ages <- c(df_inputs$age_start[1], df_inputs$age_end)
  return(v_ages)
}