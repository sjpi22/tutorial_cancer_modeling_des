# Calculate age-specific prevalence from longitudinal data
calc_prevalence <- function(m_patients, start_var, end_var, censor_var, v_ages) {
  # Create age range data frame
  res <- data.frame(
    age_start = v_ages[-length(v_ages)],
    age_end = v_ages[-1]
  )
  
  # Calculate prevalence in the given age ranges
  res <- cbind(res, t(mapply(
    function(age_start, age_end) {
      denom <- unname(unlist(m_patients[, sum(pmax(pmin(get(censor_var), age_end) - age_start, 0))]))
      num <- unname(unlist(m_patients[, sum(pmax(pmin(get(end_var), get(censor_var), age_end, na.rm = TRUE) - pmax(pmin(get(start_var), Inf, na.rm = T), age_start, na.rm = TRUE), 0))]))
      return(c(
        person_years_cases = num, 
        person_years_total = denom, 
        value = num/denom)
      )
    }, res$age_start, res$age_end))) %>%
    mutate_all(~replace(., is.na(.), 0)) %>%
    setDT()
  
  return(res)
}

# Calculate age-specific prevalence cross-sectionally
# Caution: Overwrites and deletes columns named 'sample_age' or 'fl_case' in m_patients
# If populated, dt_sample_ages should be a data.table with columns 'pt_id' and 'sample_age'
calc_prevalence_cs <- function(m_patients, start_var, end_var, censor_var, v_ages, 
                               sample_distr = NULL,
                               dt_sample_ages = NULL,
                               conf_level = 0.95) {
  # Create age range data frame
  dt_ages <- data.table(
    age_start = v_ages[-length(v_ages)],
    age_end = v_ages[-1]
  )
  
  # Sample age of study
  if (is.null(dt_sample_ages)) {
    if (is.null(sample_distr)) {
      # If no sample distribution is provided, sample uniformly from the age range
      sample_distr <- list(
        distr = "unif",
        params = list(min = min(v_ages), max = max(v_ages))
      )
    }
    dt_sample_ages <- data.table(
      pt_id = m_patients$pt_id,
      sample_age = query_distr("r", nrow(m_patients), sample_distr$distr, sample_distr$params)
    )
  }
  
  # Get age range of sample age
  dt_sample_ages[, age_idx := findInterval(sample_age, v_ages)]
  dt_sample_ages[, age_start := v_ages[age_idx], by = age_idx]
  
  # Set patient ID as key of sample age data table
  setkey(dt_sample_ages, pt_id)
  
  # Merge sample age to patient data table
  m_patients[dt_sample_ages, `:=` (sample_age = i.sample_age,
                                   age_start = i.age_start)]
  
  # Calculate cross-sectional prevalence by age group among people not censored by sample age
  m_patients[get(censor_var) > sample_age, fl_case := (get(start_var) <= sample_age & get(end_var) > sample_age)]
  m_patients[get(censor_var) > sample_age & is.na(fl_case), fl_case := F] # Reset NA to FALSE
  summ_prevalence <- m_patients[get(censor_var) > sample_age, .(
    n_total = .N,
    n_cases = sum(fl_case),
    value = mean(fl_case)), by = age_start]
  summ_prevalence <- merge(dt_ages, summ_prevalence, by = "age_start", all.x = T)
  
  # Calculate confidence intervals and merge to summary table
  df_confint <- data.frame(t(mapply(function(x, y) prop.test(x, y, conf.level = conf_level)$conf.int, summ_prevalence$n_cases, summ_prevalence$n_total)))
  summ_prevalence[, c("ci_lb", "ci_ub") := df_confint]
  
  # Estimate SE from CI
  summ_prevalence[, se := (ci_ub - ci_lb) / (2 * qnorm((1 + conf_level)/2))]
  
  # Remove added variables from m_patients
  m_patients[, c("sample_age", "age_start", "fl_case") := NULL]
  return(summ_prevalence)
}


# Calculate number of lesions from longitudinal data (assuming uniform weighting of times during which lesion is present)
calc_nlesions <- function(m_lesions, start_var, end_var, censor_var, 
                          start_age, end_age, n_max = 3) {
  # Account for case of null data
  if (!is.null(m_lesions)) {
    # Convert data to lesion start and end events
    dt_events <- rbindlist(list(
      m_lesions[get(start_var) < pmin(get(censor_var), end_age), .(pt_id, event_time = pmax(get(start_var), start_age), delta = 1)],  # Start of lesion or eligible screening period
      m_lesions[get(start_var) < pmin(get(censor_var), end_age), .(pt_id, event_time = pmin(get(end_var), get(censor_var), end_age), delta = -1)] # End of lesion or eligible screening period
    ))
    
    # Sort by patient and event time
    setorder(dt_events, pt_id, event_time)
    
    # Compute cumulative lesion count and time intervals
    dt_events[, `:=` (lesion_count = cumsum(delta),
                      next_time = shift(event_time, type = "lead")),
              by = pt_id]
    
    # Filter to time intervals with > 0 lesions
    dt_intervals <- dt_events[!is.na(next_time) & lesion_count > 0, .(pt_id, start_time = event_time, end_time = next_time, lesion_count)]
    
    # Get duration of time intervals
    dt_intervals[, duration := end_time - start_time]
    
    # Collapse together any number of lesions larger than n_max (e.g., if n_max is 3, 3+ is grouped together)
    dt_intervals[, n_lesions := ifelse(lesion_count > n_max, n_max, lesion_count)]
    
    # Sum amount of time per lesion count
    lesion_cts <- dt_intervals[, .(person_years_total = sum(duration)), by = n_lesions]
    
    # Get percentage of time with each lesion count
    lesion_cts[, person_years_cases := sum(person_years_total)]
    lesion_cts[, value := person_years_total / person_years_cases]
    
    # Ensure that all numbers from 0 to n_max are represented with count of 0 if necessary
    if (nrow(lesion_cts) < n_max) {
      lesion_cts <- data.frame(n_lesions = seq(n_max)) %>%
        left_join(lesion_cts, by = "n_lesions") %>%
        mutate_all(~replace(., is.na(.), 0)) %>%
        setDT()
    } else {
      # Sort
      lesion_cts <- lesion_cts[order(n_lesions)]
    }
  } else {
    lesion_cts <- data.table(n_lesions = seq(n_max),
                             value = 0)
  }
  
  # Make categorical variable for lesion count
  lesion_cts[, n_lesions_cat := as.factor(ifelse(n_lesions == n_max, paste0(n_max, "+"), as.character(n_lesions)))]
  return(lesion_cts)
}


# Calculate number of lesions cross-sectionally
calc_nlesions_cs <- function(m_lesions, start_var, end_var, censor_var, 
                             start_age, end_age, n_max = 3, 
                             sample_distr = NULL, dt_sample_ages = NULL, 
                             conf_level = 0.95) {
  # Account for case of null data
  if (!is.null(m_lesions)) {
    # Sample age of study
    if (is.null(dt_sample_ages)) {
      if (is.null(sample_distr)) {
        # If no sample distribution is provided, sample uniformly from the age range
        sample_distr <- list(
          distr = "unif",
          params = list(min = start_age, max = end_age)
        )
      }
      dt_sample_ages <- data.table(
        pt_id = unique(m_lesions$pt_id)
      )
      dt_sample_ages[, sample_age := query_distr("r", nrow(dt_sample_ages), sample_distr$distr, sample_distr$params)]
    }
    
    # Set patient ID as key of sample age data table
    setkey(dt_sample_ages, pt_id)
    
    # Merge sample age to lesion data table
    m_lesions[dt_sample_ages, `:=` (sample_age = i.sample_age)]
    
    # Filter to lesions existing at age and uncensored, and get distribution of lesion frequency
    lesion_cts <- m_lesions[
      get(start_var) <= sample_age & get(end_var) > sample_age & get(censor_var) > sample_age, 
      .(lesion_count = .N), by = pt_id
    ][, .(n_cases = .N), by = lesion_count]
    
    # Collapse together any number of lesions larger than n_max (e.g., if n_max is 3, 3+ is grouped together)
    lesion_cts[, n_lesions := ifelse(lesion_count > n_max, n_max, lesion_count)]
    
    # Sum number of patients for each category for number of lesions
    lesion_cts <- lesion_cts[, .(n_cases = sum(n_cases)), by = n_lesions]
    
    # Get percentage of patients with each lesion count
    lesion_cts[, n_total := sum(n_cases)]
    lesion_cts[, value := n_cases / n_total]
    
    # Calculate confidence intervals and merge to summary table
    df_confint <- data.frame(t(mapply(function(x, y) prop.test(x, y, conf.level = conf_level)$conf.int, lesion_cts$n_cases, lesion_cts$n_total)))
    lesion_cts[, c("ci_lb", "ci_ub") := df_confint]
    
    # Estimate SE from CI
    lesion_cts[, se := (ci_ub - ci_lb) / (2 * qnorm((1 + conf_level)/2))]
    
    # Ensure that all numbers from 0 to n_max are represented with count of 0 if necessary
    if (nrow(lesion_cts) < n_max) {
      lesion_cts <- data.frame(n_lesions = seq(n_max)) %>%
        left_join(lesion_cts, by = "n_lesions") %>%
        mutate_all(~replace(., is.na(.), 0)) %>%
        setDT()
    } else {
      # Sort
      lesion_cts <- lesion_cts[order(n_lesions)]
    }
  } else {
    lesion_cts <- data.table(n_lesions = seq(n_max),
                             value = 0)
  }
  
  # Make categorical variable for lesion count
  lesion_cts[, n_lesions_cat := as.factor(ifelse(n_lesions == n_max, paste0(n_max, "+"), as.character(n_lesions)))]
  return(lesion_cts)
}


# Calculate annual incidence from longitudinal data
calc_incidence <- function(m_patients, time_var, censor_var, v_ages,
                           strat_var = NULL,
                           rate_unit = 100000) {
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
    person_years_total = mapply(
      function(age_start, age_end) {
        person_years_total = m_patients[, sum(pmax(pmin(get(censor_var), age_end) - age_start, 0))]
      }, age_df$age_start, age_df$age_end)
    ) %>%
    filter(age_range %in% age_df$age_range[age_df$age_start %in% v_ages[-length(v_ages)]])
  
  # Set grouping variables
  grouping_vars <- "age_range"
  if(!is.null(strat_var)) grouping_vars <- c(grouping_vars, strat_var)
  
  # Set rate unit
  if(is.null(rate_unit)) rate_unit <- 1
  
  # Number of events by category
  event_counts <- m_patients[get(time_var) <= pmin(get(censor_var), max(age_bounds))][
    , age_range := age_ranges[findInterval(get(time_var), age_bounds)]
  ][, .N, by = grouping_vars]
  setnames(event_counts, "N", "n_events")
  
  # Make sure all values of stratifying variables are represented
  if(!is.null(strat_var)) {
    # Get sorted values of stratifying variable
    strat_var_vals <- sort(unique(m_patients[[strat_var]]))
    
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
        value = n_events / person_years_total * rate_unit,
        se = sqrt(n_events) / person_years_total * rate_unit
      ) %>%
      setDT()
  } else {
    event_counts <- person_years_at_risk %>%
      mutate(n_events = 0,
             unit = rate_unit,
             value = 0,
             se = 0) %>%
      setDT()
  }
  
  return(event_counts)
}


# Calculate distribution (such as cancer stage or lesion type) from patient-level data
calc_distr <- function(m_patients, grouping_var, event_var, censor_var, 
                       groups_expected, min_age = 0, conf_level = 0.95) {
  # Count patients diagnosed at each stage
  cts <- m_patients[get(event_var) >= min_age & get(event_var) < get(censor_var), 
                    .(n_cases = .N), by = grouping_var]
  
  # Create dataframe of expected stages if necessary
  if (nrow(cts) < length(groups_expected)) {
    df_groups_expected <- setNames(data.frame(groups_expected), grouping_var)
    
    # If no patients diagnosed at any stage, create dataframe of zeros
    if (nrow(cts) == 0) {
      cts <- df_groups_expected %>%
        mutate(n_cases = 0) %>%
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
  
  if (sum(cts$n_cases) == 0) {
    cts[, c("n_total", "value", "ci_lb", "ci_ub", "se") := 0]
  } else {
    # Get percentage of patients at each stage
    cts[, n_total := sum(n_cases)]
    cts[, value := n_cases / n_total]
    
    # Calculate confidence intervals and merge to summary table
    df_confint <- data.frame(t(mapply(function(x, y) prop.test(x, y, conf.level = conf_level)$conf.int, cts$n_cases, cts$n_total)))
    cts[, c("ci_lb", "ci_ub") := df_confint]
    
    # Estimate SE from CI
    cts[, se := (ci_ub - ci_lb) / (2 * qnorm((1 + conf_level)/2))]
  }
  
  return(cts)
}


# Calculate total life years in cohort among those who are alive and disease-free above a minimum age
calc_lifeyears <- function(
    m_patients, 
    sum_var = "time_H_D",
    censor_var = "time_screen_censor",
    age_min = 0,
    unit = 1
) {
  # Sum life years among people above the age cutoff
  if (age_min > 0 | !is.null(age_min)) {
    res <- m_patients[get(censor_var) >= age_min, .(time_total = sum(get(sum_var)),
                                                    N = .N)]
  } else {
    res <- m_patients[, .(time_total = sum(get(sum_var)),
                          N = .N)]
  }
  
  # Scale to unit if necessary
  if (unit == 1) {
    return(unlist(res))
  } else {
    res[, time_total := time_total / N * unit]
    return(unlist(res))
  }
}


# Calculate mean dwell time (can also be used to calculate mean sojourn time, which is done by default)
calc_dwell_time <- function(
    m_patients, 
    start_var = "time_H_P",
    end_var = "time_H_C",
    event_var = "time_H_C",
    censor_var = "time_H_D"
) {
  # Get mean sojourn time among people diagnosed with cancer in lifetime
  res <- m_patients[get(event_var) < get(censor_var), mean(get(end_var) - get(start_var))]
  return(res)
}