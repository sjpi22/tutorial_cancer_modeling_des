# Sample patient cohort for screening if 'idx_screen_sample' or 
# 'n_screen_sample" is populated
screen_sample <- function(m_patients, n_screen_sample = NULL, idx_screen_sample = NULL) {
  if (!is.null(idx_screen_sample)) {
    return(m_patients[idx_screen_sample])
  } else if (!is.null(n_screen_sample)) {
    if (n_screen_sample < nrow(m_patients)) return(m_patients[sample(.N, n_screen_sample)])
  }
  return(m_patients)
}

# Calculate age-specific prevalence
calc_prevalence <- function(m_patients, start_var, end_var, censor_var, v_ages, 
                            n_screen_sample = NULL,
                            idx_screen_sample = NULL,
                            conf_level = 0.95) {
  # Create age range data frame
  res <- data.frame(
    age_start = v_ages[-length(v_ages)],
    age_end = v_ages[-1]
  )
  
  # Sample patient cohort for screening if applicable
  m_screen_sample <- screen_sample(m_patients, n_screen_sample, idx_screen_sample)
  
  # Calculate prevalence in the given age ranges
  res <- cbind(res, t(mapply(
    function(age_start, age_end) {
      denom <- unname(unlist(m_screen_sample[, sum(pmax(pmin(get(censor_var), age_end) - age_start, 0))]))
      num <- unname(unlist(m_screen_sample[, sum(pmax(pmin(get(end_var), get(censor_var), age_end, na.rm = TRUE) - pmax(pmin(get(start_var), Inf, na.rm = T), age_start, na.rm = TRUE), 0))]))
      return(c(
        person_years_cases = num, 
        person_years_total = denom, 
        value = num/denom)
      )
    }, res$age_start, res$age_end))) %>%
    mutate(
      se = sqrt(value * (1-value) / person_years_total),
      ci_lb = qbinom((1-conf_level)/2, size = round(person_years_total), prob = value) / round(person_years_total),
      ci_ub = qbinom((1+conf_level)/2, size = round(person_years_total), prob = value) / round(person_years_total)
    ) %>%
    mutate_all(~replace(., is.na(.), 0))
  
  return(res)
}


# Calculate number of lesions
calc_nlesions <- function(m_lesion, start_var, end_var, censor_var, age, n_max = 3) {
  # Account for case of null data
  if (!is.null(m_lesion)) {
    # Filter to lesions existing at age and uncensored, and get distribution of lesion frequency
    lesion_cts <- m_lesion[
      get(start_var) <= age & get(end_var) > age & get(censor_var) > age, 
      .(num_lesion = .N), by = pt_id
    ][, .(ct_per_num = .N), by = num_lesion]
    
    # Collapse together any number of lesions larger than n_max (e.g., if n_max is 3, 3+ is grouped together)
    lesion_cts[num_lesion > n_max, num_lesion := n_max]
    
    # Sum number of patients for each category for number of lesions
    lesion_cts <- lesion_cts[, .(ct = sum(ct_per_num)), by = num_lesion]
    
    # Get percentage of patients with each lesion count
    lesion_cts[, value := ct / sum(ct)]
    
    # Calculate standard error
    lesion_cts[, se := sqrt(value * (1-value) / sum(ct))]
    
    # Ensure that all numbers from 0 to n_max are represented with count of 0 if necessary
    if (nrow(lesion_cts) < n_max) {
      lesion_cts <- data.frame(num_lesion = seq(n_max)) %>%
        left_join(lesion_cts, by = "num_lesion") %>%
        mutate_all(~replace(., is.na(.), 0))
    } else {
      # Sort
      lesion_cts <- lesion_cts[order(num_lesion)]
    }
  } else {
    lesion_cts <- data.frame(num_lesion = seq(n_max),
                             value = 0)
  }
  return(lesion_cts)
}


# Calculate annual incidence
calc_incidence <- function(m_patients, time_var, censor_var, v_ages,
                           strat_var = NULL,
                           rate_unit = 100000,
                           n_screen_sample = NULL,
                           idx_screen_sample = NULL) {
  
  # Sample patient cohort for screening if applicable
  m_screen_sample <- screen_sample(m_patients, n_screen_sample, idx_screen_sample)
  
  # Create age category labels
  age_bounds <- unique(c(0, v_ages)) # Add 0 as lower bound if necessary
  age_ranges <- paste(age_bounds[-length(age_bounds)], 
                      age_bounds[-1], sep = "-")
  age_df <- data.frame(age_range = age_ranges, 
                       age_start = age_bounds[-length(age_bounds)],
                       age_end = age_bounds[-1])
  
  # Exposure time by age range
  person_years_at_risk <- data.frame(
    age_df, 
    total_atrisk = mapply(
      function(age_start, age_end) {
        total_atrisk = m_screen_sample[, sum(pmax(pmin(get(censor_var), age_end) - age_start, 0))]
      }, age_df$age_start, age_df$age_end),
    age_diff = age_df$age_end - age_df$age_start) %>%
    filter(age_range %in% age_df$age_range[age_df$age_start %in% v_ages[-length(v_ages)]])
  
  # Set grouping variables
  grouping_vars <- "age_range"
  if(!is.null(strat_var)) grouping_vars <- c(grouping_vars, strat_var)
  
  # Set rate unit
  if(is.null(rate_unit)) rate_unit <- 1
  
  # Number of events by category
  event_counts <- m_screen_sample[get(time_var) <= pmin(get(censor_var), max(age_bounds))][
    , age_range := age_ranges[findInterval(get(time_var), age_bounds)]
  ][, .N, by = grouping_vars]
  setnames(event_counts, "N", "n_events")
  
  # Make sure all values of stratifying variables are represented
  if(!is.null(strat_var)) {
    # Get sorted values of stratifying variable
    strat_var_vals <- sort(unique(m_screen_sample[[strat_var]]))
    
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
      mutate(
        unit = rate_unit,
        value = n_events / total_atrisk * rate_unit,
        se = sqrt(n_events) / total_atrisk * rate_unit
      )
  } else {
    event_counts <- person_years_at_risk %>%
      mutate(n_events = 0,
             unit = rate_unit,
             value = 0,
             se = 0)
  }
  
  return(event_counts)
}


# Calculate distribution (such as cancer stage or lesion type) from patient-level data
calc_distr <- function(m_patients, grouping_var, event_var, censor_var, 
                       groups_expected, min_age = 0, n_screen_sample = NULL, 
                       idx_screen_sample = NULL,
                       conf_level = 0.95) {
  
  # Sample patient cohort for screening if applicable
  m_screen_sample <- screen_sample(m_patients, n_screen_sample, idx_screen_sample)
  
  # Count patients diagnosed at each stage
  cts <- m_screen_sample[get(event_var) >= min_age & get(event_var) < get(censor_var), 
                         .(ct = .N), by = grouping_var]
  
  # Create dataframe of expected stages if necessary
  if (nrow(cts) < length(groups_expected)) {
    df_groups_expected <- setNames(data.frame(groups_expected), grouping_var)
    
    # If no patients diagnosed at any stage, create dataframe of zeros
    if (nrow(cts) == 0) {
      cts <- df_groups_expected %>%
        mutate(ct = 0) %>%
        setDT()
    } else {
      # Ensure that all stages are represented with count of 0 if necessary
      cts <- df_groups_expected %>%
        left_join(cts, by = grouping_var) %>%
        mutate_all(~replace(., is.na(.), 0)) %>%
        setDT()
    }
  } else {
    # Order
    cts <- cts[order(get(grouping_var))]
  }
  
  if (sum(cts$ct) == 0) {
    cts[, value := 0]
  } else {
    # Get percentage of patients at each stage
    cts[, value := ct / sum(ct)]
  }
  
  # Add SE and confidence intervals
  cts <- cts %>%
    mutate(
      se = sqrt(value * (1-value) / sum(ct)),
      ci_lb = qbinom((1-conf_level)/2, size = sum(ct), prob = value) / sum(ct),
      ci_ub = qbinom((1+conf_level)/2, size = sum(ct), prob = value) / sum(ct)
    ) %>%
    mutate_all(~replace(., is.na(.), 0))
  
  return(cts)
}


# Life years in intervention scenario vs. baseline
calc_lifeyears <- function(
    m_patients, 
    sum_var,
    censor_var = NULL,
    age_min = NULL,
    unit = 1
) {
  # Sum life years among people above the age cutoff
  if (!is.null(age_min)) {
    res <- m_patients[get(censor_var) >= age_min, .(sum(get(sum_var)))]
  } else {
    res <- m_patients[, .(sum(get(sum_var)))]
  }
  return(res / unit)
}