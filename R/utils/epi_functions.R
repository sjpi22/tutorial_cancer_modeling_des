
# Calculate prevalence
calc_prevalence <- function(m_time, start_var, end_var, censor_var) {
  # Number in cohort
  n_cohort <- nrow(m_time)
  
  # Flag start times, end times, and death times
  m_time_ordered <- copy(m_time)
  
  start_events <- m_time_ordered[get(start_var) < get(censor_var), ] %>%
    select(all_of(start_var)) %>%
    mutate(delta_condition = 1) %>%
    rename(time = start_var)
  
  end_events <- m_time[get(start_var) < get(censor_var), ] %>%
    mutate(time = pmin(get(end_var), get(censor_var)),
           delta_condition = -1) %>%
    select(time, delta_condition)
  
  censor_events <- m_time %>%
    select(all_of(censor_var)) %>%
    mutate(delta_total = -1) %>%
    rename(time = censor_var)
  
  # Combine and calculate population totals over time
  all_events <- rbind(start_events, end_events, censor_events, fill=TRUE) %>%
    arrange(by = time) %>%
    mutate_all(~replace(., is.na(.), 0)) %>%
    # Collapse by time in case of multiple events at the same time
    group_by(time) %>%
    summarise(delta_condition = sum(delta_condition),
              delta_total = sum(delta_total)) %>%
    ungroup() %>%
    mutate(n_condition = cumsum(delta_condition),
           n_total = n_cohort + cumsum(delta_total)) %>%
    mutate(prevalence = n_condition / n_total) 
  
  # Calculate prevalence closest to each discrete age
  v_ages <- seq(ceiling(min(all_events$time)),
                max(all_events$time[!is.infinite(all_events$time)]))
  indices <- with(all_events, {
    indices <- sapply(v_ages,
                      function(x) which.max(time[time <= x]))
  }
  )
  
  df_prevalence <- cbind(age = v_ages, 
                         prevalence = all_events$prevalence[indices],
                         n_population = all_events$n_total[indices])
  
  return(df_prevalence)
}

# Calculate annual incidence
calc_incidence <- function(m_time, time_var, censor_var, age_upper_bounds,
                           strat_var = NULL,
                           rate_unit = 100000) {
  
  # Number in cohort
  n_cohort <- nrow(m_time)
  
  # Create age labels
  m_time_ordered <- copy(m_time)
  m_time_ordered[, `:=`(year_censor = floor(get(censor_var)),
                        year_event = floor(get(time_var)))]
  
  # Create age category labels
  age_lower_bounds <- lag(age_upper_bounds, default = 0)
  age_ranges <- paste(age_lower_bounds, 
                      age_upper_bounds-1, sep = "-")
  age_ranges <- c(age_ranges, paste0(max(age_upper_bounds), "+"))
  age_bounds <- c(age_lower_bounds, max(age_upper_bounds))
  
  # Deaths by age
  censor_counts_year <- m_time_ordered %>%
    group_by(year_censor) %>%
    summarise(n_died = n()) %>%
    ungroup() %>%
    mutate(age_range = sapply(year_censor,
                              function(x) age_ranges[which.max(age_bounds[x >= age_bounds])]),
           n_start = n_cohort - cumsum(c(0, n_died[-length(n_died)])))
  
  # Deaths by category
  censor_counts_range <- censor_counts_year %>%
    group_by(age_range) %>%
    summarise(person_years_at_risk = sum(n_start))
  
  # Set grouping variables
  group_vars <- "age_range"
  if(!is.null(strat_var)) group_vars <- c(group_vars, strat_var)
  
  # Set rate unit
  if(is.null(rate_unit)) rate_unit <- 1
  
  # Events by category and stage
  event_counts <- m_time_ordered %>%
    filter(get(time_var) <= get(censor_var)) %>%
    mutate(age_range = sapply(year_event,
                              function(x) age_ranges[which.max(age_bounds[x >= age_bounds])])) %>%
    group_by_at(group_vars) %>%
    summarise(n_events = n()) %>%
    ungroup() 
  
  if(nrow(event_counts) > 0) {
    event_counts <- censor_counts_range %>%
      left_join(event_counts, by = "age_range") %>%
      mutate(n_events = replace_na(n_events, 0)) %>%
      mutate(incidence = n_events / person_years_at_risk * rate_unit)
  } else {
    event_counts <- censor_counts_range %>%
      mutate(n_events = 0,
             incidence = 0)
  }
  
  return(event_counts)
}

