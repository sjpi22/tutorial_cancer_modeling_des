#' Calculate Age-Specific Prevalence from Longitudinal Data
#'
#' Estimates the prevalence of a condition in specified age ranges from a matrix or 
#' data frame of patient-level time-to-event data.
#'
#' @param m_patients A data frame or matrix where each row represents a patient, 
#'   and columns include age at event onset, transition, and censoring.
#' @param start_var The name (string) of the variable in \code{m_patients} 
#'   indicating age at condition onset.
#' @param end_var The name (string) of the variable indicating age at which 
#'   condition ends (e.g., progression to cancer).
#' @param censor_var The name (string) of the variable indicating age at censoring 
#'   (e.g., death or loss to follow-up). Patients who are censored before a given 
#'   age are excluded from the denominator at that age.
#' @param id_var The name (string) of the variable indicating the patient ID.
#' @param v_ages A numeric vector of ages at which to estimate prevalence. If \code{NULL}, 
#'   prevalence will be calculated across the age range of the population.
#' @param method Method for calculating prevalence. \code{"cs"} for cross- 
#'   sectional, \code{"long"} for longitudinal, and \code{"rcs"} for repeated 
#'   cross-sectional.
#' @param sample_var Optional string with the variable indicating age at which 
#'   patients are sampled for cross-sectional calculation. If \code{sample_var}
#'   and \code{dt_sample_ages} are both populated, \code{sample_var} is used.
#' @param dt_sample_ages Optional data frame with patient ID variable and 
#'   sampled patient ages that must be under variable name \code{sample_age} for
#'   cross-sectional calculation. If \code{sample_var} and \code{dt_sample_ages}
#'   are both populated, \code{sample_var} is used.
#' @param output_uncertainty Binary indicator for whether to output standard
#'   errors and confidence intervals. Formulas are still in development mode
#'   for longitudinal and repeated cross-sectional formulations.
#' @param conf_level Confidence level for binomial confidence intervals, default is 0.95.
#'
#' @return A data table with estimated prevalence and confidence intervals at each age.
#'
#' @examples
#' # Example usage with simulated data:
#' m_patients <- data.table(
#'   pt_id = 1:5,
#'   time_S = c(45, 50, 60, 55, 47), # Time to disease onset
#'   time_R = c(52, 65, 70, 58, 55), # Time to recovery
#'   time_D = c(50, 75, 68, 60, 90)  # Time to death
#' )
#' calc_prevalence(
#'   m_patients = m_patients,
#'   start_var = "time_S",
#'   end_var = "time_R",
#'   censor_var = "time_D",
#'   v_ages = seq(50, 60, by = 5)
#' )
#' 
#' @import data.table
#' @importFrom dplyr %>%
#' @export
calc_prevalence <- function(m_patients, 
                            start_var, 
                            end_var, 
                            censor_var, 
                            id_var = "pt_id",
                            v_ages = NULL,
                            method = "cs",
                            sample_var = NULL,
                            dt_sample_ages = NULL,
                            output_uncertainty = FALSE,
                            conf_level = 0.95) {
  # If v_ages is NULL, calculate prevalence across entire age range of population
  if (is.null(v_ages)) {
    v_ages <- c(0, m_patients[, max(get(censor_var))])
  }
  
  # Create age range data frame
  dt_ages <- data.table(
    age_start = v_ages[-length(v_ages)],
    age_end = v_ages[-1]
  )
  
  # Ensure that key is set for cohort data
  if (is.null(key(m_patients))) setkeyv(m_patients, id_var)
  
  # Cross-sectional vs. longitudinal formulation
  if (method == "cs") {
    # Save variables that will be overwritten
    l_overwrite <- list()
    for (varname in c("sample_age", "age_start", "fl_case")) {
      # Save sample_age variable if it will be overwritten
      if (varname %in% colnames(m_patients)) {
        # Exception for sample_age when sample_var is "sample_age"
        if (!(varname == "sample_age" & ifelse(is.null(sample_var), "", sample_var) == "sample_age")) {
          l_overwrite[[varname]] <- copy(m_patients[[varname]])
        }
      }
    }
    
    # Assign sample ages if not given
    if (is.null(sample_var)) {
      # Sample age of study
      if (is.null(dt_sample_ages)) {
        dt_sample_ages <- data.table(pt_id = m_patients[[id_var]])
        dt_sample_ages[, sample_age := runif(.N, min(v_ages), max(v_ages))]
      }
      
      # Get age range of sample age
      dt_sample_ages[, age_idx := findInterval(sample_age, v_ages)]
      dt_sample_ages[, age_start := v_ages[age_idx], by = age_idx]
      
      # Set patient ID as key of sample age data table
      setkeyv(dt_sample_ages, id_var)
      
      # Rename ID if necessary
      if (id_var != "pt_id") setnames(dt_sample_ages, "pt_id", id_var)
      
      # Merge sample age to patient data table
      m_patients[dt_sample_ages, `:=` (sample_age = i.sample_age,
                                       age_start = i.age_start)]
    } else {
      # Rename given sample variable
      if (sample_var != "sample_age") {
        if (!is.null(m_patients[["sample_age"]])) m_patients[, sample_age := NULL] # Remove sample_age variable from data table if already present
        setnames(m_patients, sample_var, "sample_age") # Replace name of sample age variable to sample_age
      }
      
      # Save age_idx variable if it will be overwritten
      if ("age_idx" %in% colnames(m_patients)) {
        age_idx_saved <- m_patients$age_idx
      }
      
      # Get age range of sample age
      m_patients[, age_idx := findInterval(sample_age, v_ages)]
      m_patients[, age_start := v_ages[age_idx], by = age_idx]
      
      # Replace age_idx variable
      if (exists("age_idx_saved")) {
        m_patients[, age_idx := age_idx_saved]
      } else {
        m_patients[, age_idx := NULL]
      }
    }
    
    # Calculate cross-sectional prevalence by age group among people not censored by sample age
    m_patients[get(censor_var) > sample_age, fl_case := (get(start_var) <= sample_age & get(end_var) > sample_age)]
    m_patients[get(censor_var) > sample_age & is.na(fl_case), fl_case := F] # Reset NA to FALSE
    summ_prevalence <- m_patients[get(censor_var) > sample_age, .(
      n_total = .N,
      n_cases = sum(fl_case),
      value = mean(fl_case)), by = age_start]
    summ_prevalence <- merge(dt_ages, summ_prevalence, by = "age_start", all.x = T)
    
    # If required to output uncertainty estimates
    if (output_uncertainty) {
      # Calculate confidence intervals and merge to summary table
      df_confint <- data.frame(t(mapply(function(x, y) prop.test(x, y, conf.level = conf_level)$conf.int, 
                                        summ_prevalence$n_cases, 
                                        summ_prevalence$n_total)))
      summ_prevalence[, c("ci_lb", "ci_ub") := df_confint]
      
      # Estimate SE from CI
      summ_prevalence[, se := (ci_ub - ci_lb) / (2 * qnorm((1 + conf_level)/2))]
    }
    
    # Reset or remove sample_age variable
    vars_remove <- c("age_start", "fl_case")
    if (!is.null(sample_var)) { # If sample age variable was already included in patient matrix
      # Reset variable name
      if (sample_var != "sample_age") { # If sample age variable was not the standardized name
        setnames(m_patients, "sample_age", sample_var) # Revert to unstandardized name
      }
    } else { # If sample age variable was not included in patient matrix
      # Add "sample_age" to list of variables to remove
      vars_remove <- c(vars_remove, "sample_age")
    }
    
    # Remove variables that were added as side effects
    vars_remove <- setdiff(vars_remove, names(l_overwrite)) # Remove variables that were added and not replaced
    vars_remove <- intersect(vars_remove, names(m_patients))  # Only remove columns that exist
    if (length(vars_remove) > 0) {
      m_patients[, c(vars_remove) := NULL]
    }
    
    # Replace saved variables that were overwritten
    for (varname in names(l_overwrite)) {
      m_patients[, (varname) := l_overwrite[[varname]]]
    }
    
  } else if (method == "long") {
    # Calculate prevalence in the given age ranges
    summ_prevalence <- cbind(dt_ages, t(mapply(
      function(age_start, age_end) {
        denom <- unname(unlist(m_patients[, sum(pmax(pmin(get(censor_var), age_end) - age_start, 0))]))
        num <- unname(unlist(m_patients[, sum(pmax(pmin(get(end_var), get(censor_var), age_end, na.rm = TRUE) - pmax(pmin(get(start_var), Inf, na.rm = T), age_start, na.rm = TRUE), 0))]))
        return(c(
          person_years_cases = num, 
          person_years_total = denom, 
          value = num/denom)
        )
      }, dt_ages$age_start, dt_ages$age_end))) %>%
      mutate_all(~replace(., is.na(.), 0)) %>%
      setDT()
    
    # If required to output uncertainty estimates
    if (output_uncertainty) {
      # Generate confidence interval
      df_confint <- data.frame(t(mapply(function(x, y) prop.test(x, y, conf.level = conf_level)$conf.int, 
                                        summ_prevalence$person_years_cases / (summ_prevalence$age_end - summ_prevalence$age_start) * 2, 
                                        summ_prevalence$person_years_total / (summ_prevalence$age_end - summ_prevalence$age_start) * 2)))
      summ_prevalence[, c("ci_lb", "ci_ub") := df_confint]
      
      # Estimate SE from CI
      summ_prevalence[, se := (ci_ub - ci_lb) / (2 * qnorm((1 + conf_level)/2))]
    }
  } else if (method == "rcs") {
    # Create sequence of age ranges including (-Inf, 0] with 1-year intervals from min to max of v_ages
    v_ages_tabulate <- unique(c(-Inf, seq(min(v_ages), max(v_ages))))
    
    # Tabulate counts in age intervals for start age, end age, and censor age
    counts_censor <- tabulate(cut(m_patients[, get(censor_var)], v_ages_tabulate))
    counts_start <- tabulate(cut(m_patients[, fifelse(get(start_var) < get(censor_var), get(start_var), NA)], v_ages_tabulate))
    counts_end <- tabulate(cut(m_patients[, fifelse(get(start_var) < get(censor_var), pmin(get(end_var), get(censor_var), na.rm = T), NA)], v_ages_tabulate))
    
    # Pad start and end counts
    counts_start <- c(counts_start, rep(0, length(v_ages_tabulate) - 1 - length(counts_start)))
    counts_end <- c(counts_end, rep(0, length(v_ages_tabulate) - 1 - length(counts_end)))
    
    # Create data table of prevalence-related counts in each age range
    dt_counts <- data.table(
      age_year = tail(v_ages_tabulate, -1), # Lower bound for all but -Inf
      age_idx = findInterval(tail(v_ages_tabulate, -1), v_ages), # Index of age range vector (left closed)
      age_idx_ub = findInterval(tail(v_ages_tabulate, -1), v_ages, left.open = T), # Index of age range vector (left open)
      person_years_cases = cumsum(counts_start) - cumsum(counts_end), # Count number of prevalence cases at beginning of each year starting from v_ages_tabulate[2]
      person_years_total = nrow(m_patients) - cumsum(counts_censor) # Count number of uncensored patients at beginning of each year starting from v_ages_tabulate[2]
    )
    
    # Map age range groups
    dt_counts[, `:=` (age_start = v_ages[age_idx]), by = age_idx]
    
    # Set key as age group followed by age_start to ensure order
    setkey(dt_counts, age_start, age_year)
    
    # Map boundary age range groups where age_idx != age_idx_ub
    dt_counts[age_idx != age_idx_ub, `:=` (age_grp_ub = v_ages[age_idx_ub]), by = age_idx_ub]
    
    # Sum counts by age group except upper boundary
    summ_prevalence <- dt_counts[age_start < max(v_ages), 
                                 .(person_years_cases = sum(person_years_cases),
                                   person_years_total = sum(person_years_total)),
                                 by = age_start]
    
    # Merge age ranges
    summ_prevalence <- merge(dt_ages, summ_prevalence, by = "age_start")
    
    # Extract boundary counts
    summ_prevalence_ub <- dt_counts[!is.na(age_grp_ub), ]
    
    # Add boundary counts to each age group
    summ_prevalence[, `:=` (person_years_cases = person_years_cases + summ_prevalence_ub$person_years_cases,
                            person_years_total = person_years_total + summ_prevalence_ub$person_years_total)]
    
    # Calculate prevalence
    summ_prevalence[, value := person_years_cases/person_years_total]
    
    # If required to output uncertainty estimates
    if (output_uncertainty) {
      # Generate confidence interval
      df_confint <- data.frame(t(mapply(function(x, y) prop.test(x, y, conf.level = conf_level)$conf.int, 
                                        summ_prevalence$person_years_cases / (summ_prevalence$age_end - summ_prevalence$age_start) * 2, 
                                        summ_prevalence$person_years_total / (summ_prevalence$age_end - summ_prevalence$age_start) * 2)))
      summ_prevalence[, c("ci_lb", "ci_ub") := df_confint]
      
      # Estimate SE from CI
      summ_prevalence[, se := (ci_ub - ci_lb) / (2 * qnorm((1 + conf_level)/2))]
    }
  }
  
  return(summ_prevalence)
}


#' Calculate Number of Lesions from Longitudinal Data (Lesion Multiplicity)
#'
#' Estimates the distribution for the number of lesions among individuals with 
#' at least one lesion.
#'
#' @param m_lesions A data frame or matrix where each row represents a lesion, 
#'   and columns include the patient ID, age at lesion onset, transition, and censoring.
#' @param start_var The name (string) of the variable in \code{m_patients} 
#'   indicating age at condition onset.
#' @param end_var The name (string) of the variable indicating age at which 
#'   condition ends (e.g., progression to cancer).
#' @param censor_var The name (string) of the variable indicating age at censoring 
#'   (e.g., death or loss to follow-up). Patients who are censored before a given 
#'   age are excluded from the denominator at that age.
#' @param id_var The name (string) of the variable indicating the patient ID.
#' @param start_age The start of the age range to evaluate lesion distribution.
#' @param end_age The end of the age range to evaluate lesion distribution.
#' @param n_max Maximum number of lesions to reportz
#' @param method Method for calculating number of lesions. \code{"cs"} for cross-sectional 
#'   and \code{"long"} for longitudinal.
#' @param dt_sample_ages Optional data frame with patient ID variable and 
#'   sampled patient ages that must be under variable name \code{sample_age} for
#'   cross-sectional calculation.
#' @param output_uncertainty Binary indicator for whether to output standard
#'   errors and confidence intervals.
#' @param conf_level Confidence level for binomial confidence intervals, default is 0.95.
#'
#' @return A data table with the proportion of individuals or time with each 
#' number of lesions.
#'
#' @examples
#' # Example usage with simulated data:
#' m_lesions <- data.table(
#'   pt_id = c(1, 1, 2, 3, 3),
#'   lesion_id = c(1, 2, 1, 1, 2),
#'   time_S = c(45, 50, 60, 55, 47), # Time to lesion onset
#'   time_R = c(52, 65, 70, 58, 55), # Time to lesion end
#'   time_D = c(50, 75, 68, 60, 90)  # Time to death
#' )
#' calc_nlesions(
#'   m_lesions = m_lesions,
#'   start_var = "time_S",
#'   end_var = "time_R",
#'   censor_var = "time_D",
#'   start_age = 45,
#'   end_age = 65
#' )
#' 
#' @import data.table
#' @importFrom dplyr %>%
#' @export
calc_nlesions <- function(m_lesions, 
                          start_var, 
                          end_var, 
                          censor_var,
                          id_var = "pt_id",
                          start_age = 0, 
                          end_age = NULL, 
                          n_max = 3,
                          method = "cs",
                          dt_sample_ages = NULL, 
                          output_uncertainty = FALSE,
                          conf_level = 0.95) {
  # Ensure that key is set for cohort data
  if (is.null(key(m_lesions))) setkeyv(m_lesions, id_var)
  
  # If end_age is NULL, set to max age in data
  if (is.null(end_age)) {
    end_age <- m_lesions[, max(get(censor_var))]
  }
  
  if (method == "cs") {
    # Account for case of null data
    if (!is.null(m_lesions)) {
      # Sample age of study
      if (is.null(dt_sample_ages)) {
        dt_sample_ages <- data.table(
          pt_id = unique(m_lesions[[id_var]])
        )
        dt_sample_ages[, sample_age := runif(.N, start_age, end_age)]
      }
      
      # Set patient ID as key of sample age data table
      setkey(dt_sample_ages, pt_id)
      
      # Rename ID if necessary
      if (id_var != "pt_id") setnames(dt_sample_ages, "pt_id", id_var)
      
      # Merge sample age to lesion data table
      m_lesions[dt_sample_ages, `:=` (sample_age = i.sample_age)]
      
      # Filter to lesions existing at age and uncensored, and get distribution of lesion frequency
      lesion_cts <- m_lesions[
        get(start_var) <= sample_age & get(end_var) > sample_age & get(censor_var) > sample_age, 
        .(lesion_count = .N), by = get(id_var)
      ][, .(n_cases = .N), by = lesion_count]
      
      # Collapse together any number of lesions larger than n_max (e.g., if n_max is 3, 3+ is grouped together)
      lesion_cts[, n_lesions := ifelse(lesion_count > n_max, n_max, lesion_count)]
      
      # Sum number of patients for each category for number of lesions
      lesion_cts <- lesion_cts[, .(n_cases = sum(n_cases)), by = n_lesions]
      
      # Get percentage of patients with each lesion count
      lesion_cts[, n_total := sum(n_cases)]
      lesion_cts[, value := n_cases / n_total]
      
      # Calculate confidence intervals and merge to summary table
      if (output_uncertainty) {
        df_confint <- data.frame(t(mapply(function(x, y) prop.test(x, y, conf.level = conf_level)$conf.int, lesion_cts$n_cases, lesion_cts$n_total)))
        lesion_cts[, c("ci_lb", "ci_ub") := df_confint]
        
        # Estimate SE from CI
        lesion_cts[, se := (ci_ub - ci_lb) / (2 * qnorm((1 + conf_level)/2))]
      }
      
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
  } else {
    # Account for case of null data
    if (!is.null(m_lesions)) {
      # Convert data to lesion start and end events
      dt_events <- rbindlist(list(
        m_lesions[get(start_var) < pmin(get(censor_var), end_age), .(get(id_var), event_time = pmax(get(start_var), start_age), delta = 1)],  # Start of lesion or eligible screening period
        m_lesions[get(start_var) < pmin(get(censor_var), end_age), .(get(id_var), event_time = pmin(get(end_var), get(censor_var), end_age), delta = -1)] # End of lesion or eligible screening period
      ))
      
      # Sort by patient and event time
      setorderv(dt_events, c(id_var, "event_time"))
      
      # Compute cumulative lesion count and time intervals
      dt_events[, `:=` (lesion_count = cumsum(delta),
                        next_time = shift(event_time, type = "lead")),
                by = get(id_var)]
      
      # Filter to time intervals with > 0 lesions
      dt_intervals <- dt_events[!is.na(next_time) & lesion_count > 0, .(get(id_var), start_time = event_time, end_time = next_time, lesion_count)]
      
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
}


#' Calculate Age-Specific Incidence from Longitudinal Data
#'
#' Estimates the incidence of a condition at specific ages from a matrix or 
#' data frame of patient-level time-to-event data.
#'
#' @param m_patients A data frame or matrix where each row represents a patient, 
#'   and columns include age at event onset, transition, and censoring.
#' @param time_var The name (string) of the variable in \code{m_patients} 
#'   indicating age at condition onset.
#' @param censor_var The name (string) of the variable indicating age at censoring 
#'   (e.g., death or loss to follow-up). Patients who are censored before a given 
#'   age are excluded from the denominator at that age.
#' @param id_var The name (string) of the variable indicating the patient ID.
#' @param strat_var Optional name (string) of variable for further stratification.
#' @param v_ages A numeric vector of ages at which to estimate incidence. If \code{NULL}, 
#'   incidence will be calculated across the age range of the population.
#' @param method Method for calculating incidence. Currently \code{"long"} for 
#'   longitudinal is implemented.
#' @param rate_unit Quantity to divide incidence rate by for reporting.
#' @param output_uncertainty Binary indicator for whether to output standard
#'   errors and confidence intervals.
#' @param conf_level Confidence level for binomial confidence intervals, default is 0.95.
#'
#' @return A data table with estimated incidence in each age range.
#'
#' @examples
#' # Example usage with simulated data:
#' m_patients <- data.table(
#'   pt_id = 1:5,
#'   time_S = c(45, 50, 60, 55, 47), # Time to disease onset
#'   time_R = c(52, 65, 70, 58, 55), # Time to recovery
#'   time_D = c(50, 75, 68, 60, 90)  # Time to death
#' )
#' calc_incidence(
#'   m_patients = m_patients,
#'   start_var = "time_S",
#'   censor_var = "time_D",
#'   v_ages = seq(50, 60, by = 5)
#' )
#' 
#' @import data.table
#' @importFrom dplyr %>%
#' @export
calc_incidence <- function(m_patients, 
                           time_var, 
                           censor_var, 
                           id_var = "pt_id",
                           strat_var = NULL,
                           v_ages = NULL,
                           method = "long",
                           rate_unit = 100000,
                           output_uncertainty = FALSE,
                           conf_level = 0.95) {
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
    # Merge event counts and person-years at risk
    event_counts <- person_years_at_risk %>%
      left_join(event_counts, by = "age_range") %>%
      mutate(n_events = replace_na(n_events, 0)) %>%
      mutate(
        unit = rate_unit,
        value = n_events / person_years_total * rate_unit,
        se = sqrt(n_events) / person_years_total * rate_unit
      ) %>%
      setDT()
    
    if (output_uncertainty) {
      # Calculate chi-squared critical values
      event_counts[, `:=` (chi2_lb = qchisq((1-conf_level)/2, 2*n_events),
                           chi2_ub = qchisq((1+conf_level)/2, 2*(n_events + 1)))]
      event_counts[, `:=` (ci_lb = chi2_lb/2/person_years_total,
                           ci_ub = chi2_ub/2/person_years_total)]
      
      # Adjust for rate unit
      if (rate_unit != 1) {
        event_counts[, `:=` (ci_lb = ci_lb * rate_unit,
                             ci_ub = ci_ub * rate_unit)]
      }
    }
  } else {
    event_counts <- person_years_at_risk %>%
      mutate(n_events = 0,
             unit = rate_unit,
             value = 0,
             se = 0) %>%
      setDT()
    
    if (output_uncertainty) {
      event_counts[, c("ci_lb", "ci_ub") := NA]
    }
  }
  
  return(event_counts)
}


#' Calculate Categorical Distribution
#'
#' Calculates a categorical distribution (such as cancer stage or primary lesion 
#' type) from patient-level data at a point in time.
#'
#' @param m_patients A data frame or matrix where each row represents a patient, 
#'   and columns include age at event onset, transition, and censoring.
#' @param grouping_var The name (string) of the variable indicating categories.
#' @param event_var The name (string) of the variable in \code{m_patients} 
#'   indicating age of event.
#' @param censor_var The name (string) of the variable indicating age at censoring 
#'   (e.g., death or loss to follow-up). Patients who are censored before the event
#'   age are excluded from the denominator.
#' @param groups_expected Vector of expected categories in case not all are 
#'   represented in the patient data.
#' @param min_age Minimum age at which to begin evaluating distribution. 
#'   Default 0.
#' @param output_uncertainty Binary indicator for whether to output standard
#'   errors and confidence intervals.
#' @param conf_level Confidence level for binomial confidence intervals, default is 0.95.
#'
#' @return A data table with proportion in each category.
#'
#' @examples
#' # Example usage with simulated data:
#' m_patients <- data.table(
#'   pt_id = 1:5,
#'   time_cancer = c(45, 50, 60, 55, 47), # Time to cancer diagnosis
#'   stage = c(1, 3, 4, 3, 1), # Stage at diagnosis
#'   time_D = c(50, 75, 68, 60, 90)  # Time to death
#' )
#' calc_distr(
#'   m_patients = m_patients,
#'   grouping_var = "stage",
#'   event_var = "time_cancer",
#'   censor_var = "time_D",
#'   groups_expected = 1:4
#' )
#' 
#' @import data.table
#' @importFrom dplyr %>%
#' @export
calc_distr <- function(m_patients, 
                       grouping_var, 
                       event_var, 
                       censor_var, 
                       groups_expected, 
                       min_age = 0, 
                       output_uncertainty = FALSE,
                       conf_level = 0.95) {
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
    cts[, c("n_total", "value") := 0]
    
    # Add CI and SE if needed
    if (output_uncertainty) {
      cts[, c("ci_lb", "ci_ub", "se") := 0]
    }
  } else {
    # Get percentage of patients at each stage
    cts[, n_total := sum(n_cases)]
    cts[, value := n_cases / n_total]
    
    # Calculate confidence intervals and merge to summary table if needed
    if (output_uncertainty) {
      df_confint <- data.frame(t(mapply(function(x, y) prop.test(x, y, conf.level = conf_level)$conf.int, cts$n_cases, cts$n_total)))
      cts[, c("ci_lb", "ci_ub") := df_confint]
      
      # Estimate SE from CI
      cts[, se := (ci_ub - ci_lb) / (2 * qnorm((1 + conf_level)/2))]
    }
  }
  
  return(cts)
}


#' Calculate Total Life Years
#'
#' Calculate total life years in cohort among those who are alive and disease-
#' free above a minimum age
#'
#' @param m_patients A data frame or matrix where each row represents a patient, 
#'   and columns include age at event onset, transition, and censoring.
#' @param sum_var The name (string) of the variable in \code{m_patients} 
#'   that will be summed.
#' @param censor_var The name (string) of the variable indicating age at censoring 
#'   (e.g., death or loss to follow-up). Patients who are censored before a given 
#'   age are excluded from the sum. If NULL, \code{sum_var} will be used.
#' @param min_age Minimum age at which to begin evaluating distribution. 
#'   Default 0.
#'
#' @return A data frame with estimated prevalence and confidence intervals at each age.
#'
#' @examples
#' # Example usage with simulated data:
#' m_patients <- data.table(
#'   pt_id = 1:5,
#'   time_D = c(50, 75, 68, 60, 90)  # Time to death
#' )
#' calc_lifeyears(
#'   m_patients = m_patients,
#'   sum_var = "time_D",
#'   censor_var = "time_D",
#'   min_age = 25
#' )
#' 
#' @import data.table
#' @importFrom dplyr %>%
#' @export
calc_lifeyears <- function(
    m_patients, 
    sum_var,
    censor_var = NULL,
    min_age = 0
) {
  # Sum life years among people above the age cutoff
  if (min_age > 0 | !is.null(min_age)) {
    if (is.null(censor_var)) censor_var <- sum_var
    res <- m_patients[get(censor_var) >= min_age, .(time_total = sum(get(sum_var)),
                                                    N = .N)]
  } else {
    res <- m_patients[, .(time_total = sum(get(sum_var)),
                          N = .N)]
  }
  
  return(unlist(res))
}


#' Calculate Mean Duration (Dwell or Sojourn Time)
#'
#' Calculate mean duration in a state, such as dwell time or sojourn time.
#'
#' @param m_patients A data frame or matrix where each row represents a patient, 
#'   and columns include age at event onset, transition, and censoring.
#' @param start_var The name (string) of the variable in \code{m_patients} 
#'   indicating age at condition onset.
#' @param end_var The name (string) of the variable indicating age at which 
#'   condition ends.
#' @param censor_var The name (string) of the variable indicating age at censoring 
#'   (e.g., death or loss to follow-up). Patients who are censored before a given 
#'   age are excluded from the denominator at that age.
#' @param event_var The name (string) of the variable indicating the event that 
#'   must occur before the censor time for an individual to be included. If 
#'   \code{NULL}, \code{start_var} will be used.
#'
#' @return A data table with the mean duration.
#'
#' @examples
#' # Example usage with simulated data:
#' m_patients <- data.table(
#'   pt_id = 1:5,
#'   time_S = c(45, 50, 60, 55, 47), # Time to disease onset
#'   time_R = c(52, 65, 70, 58, 55), # Time to recovery
#'   time_D = c(50, 75, 68, 60, 90)  # Time to death
#' )
#' calc_dwell_time(
#'   m_patients = m_patients,
#'   start_var = "time_S",
#'   end_var = "time_R",
#'   censor_var = "time_D"
#' )
#' 
#' @import data.table
#' @importFrom dplyr %>%
#' @export
calc_dwell_time <- function(
    m_patients, 
    start_var,
    end_var,
    censor_var,
    event_var = NULL
) {
  # Get mean sojourn time among people diagnosed with cancer in lifetime
  if (is.null(event_var)) event_var <- start_var
  res <- m_patients[get(event_var) < get(censor_var), mean(get(end_var) - get(start_var))]
  return(res)
}

#' Calculate Risk
#' 
#' Calculate age-conditional or lifetime risk of a condition
#' 
#' @param m_patients A data frame or matrix where each row represents a patient, 
#'   and columns include age at event onset, transition, and censoring.
#' @param start_var The name (string) of the variable in \code{m_patients} 
#'   indicating the start age of the condition.
#' @param censor_var The name (string) of the variable indicating age at censoring 
#'   (e.g., death or loss to follow-up). Patients who are censored before a given 
#'   age are excluded from the sum. If NULL, \code{sum_var} will be used.
#' @param min_age Minimum age at which to begin evaluating risk. Default 0.
#' @param max_age Maximum age at which to evaluate risk.
#' @param output_uncertainty Binary indicator for whether to output standard
#'   errors and confidence intervals.
#' @param conf_level Confidence level for binomial confidence intervals, default is 0.95.
#'
calc_risk <- function(
    m_patients, 
    start_var,
    censor_var,
    min_age = 0,
    max_age = NULL,
    output_uncertainty = FALSE,
    conf_level = 0.95
) {
  # Set minimum age as 0 and maximum age as maximum observed age in data if not provided
  if (is.null(min_age)) min_age <- 0
  if (is.null(max_age)) max_age <- m_patients[, max(get(censor_var))]
  
  # Among individuals who are uncensored without condition by min_age,
  # calculate number that develop condition before earliest of death and max_age
  res <- m_patients[pmin(get(start_var), get(censor_var), na.rm = T) >= min_age, 
                    .(n_cases = sum(get(start_var) < pmin(max_age, get(censor_var))),
                      n_total = .N)]
  
  # Calculate risk
  res[, value := n_cases/n_total]
  
  # Calculate confidence intervals and merge to summary table if needed
  if (output_uncertainty) {
    df_confint <- data.frame(t(mapply(function(x, y) prop.test(x, y, conf.level = conf_level)$conf.int, res$n_cases, res$n_total)))
    res[, c("ci_lb", "ci_ub") := df_confint]
    
    # Estimate SE from CI
    res[, se := (ci_ub - ci_lb) / (2 * qnorm((1 + conf_level)/2))]
  }
}
