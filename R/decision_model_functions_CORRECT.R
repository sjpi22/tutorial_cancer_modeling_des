# Initialize time to event matrix for patient cohort with patient IDs and baseline strategy
initialize_cohort <- function(n_cohort) {
  m_times <- data.table(
    pt_id = 1:n_cohort
  )
  setkey(m_times, pt_id)
  return(m_times)
}

# Generate baseline data
simulate_baseline_data <- function(m_times,
                                   ages_H_Do,
                                   prob_H_Do) {
  
  # Sample time to death from other causes
  m_times[, time_H_Do := sample(x = ages_H_Do, 
                                size = .N, 
                                replace = TRUE, 
                                prob = prob_H_Do) + runif(.N)]
}


# Simulate time to cancer onset
simulate_cancer_onset <- function(m_times, shape_H_P, scale_H_P) {
  m_times[, time_H_P := rweibull(.N, shape = shape_H_P, scale = scale_H_P)]
}

# Simulate cancer stage progression and diagnosis
simulate_cancer_progression <- function(m_times,
                                        rate_P1_P2,
                                        rate_P1_C1,
                                        rate_P2_C2) {
  # Time from first to second stage
  m_times[, time_P1_P2 := rexp(.N, rate_P1_P2)]
  
  # Time to symptomatic detection in first stage
  m_times[, time_P1_C1 := rexp(.N, rate_P1_C1)]
  
  # Time to symptomatic detection in second stage
  m_times[, time_P2_C2 := rexp(.N, rate_P2_C2)]
  
  # Stage at cancer diagnosis - earliest stage at which time to detection is less than time to next stage
  m_times[, stage_dx := fifelse(time_P1_C1 < time_P1_P2, 1, 2)]
  m_times[, time_P_C := fifelse(time_P1_C1 < time_P1_P2, time_P1_C1, time_P1_P2 + time_P2_C2)]
  
  # Time to cancer diagnosis # @@@ ERROR - time_P_C was switched to time_H_P
  m_times[, time_H_C := time_H_P + time_P_C]
}

# Simulate time from cancer diagnosis to death
simulate_cancer_mortality <- function(m_times,
                                      ages_C1_Dc,
                                      prob_C1_Dc,
                                      ages_C2_Dc,
                                      prob_C2_Dc) {
  # Cancer mortality after stage 1 diagnosis
  m_times[stage_dx == 1, time_C_Dc := sample(x = ages_C1_Dc, 
                                             size = .N, 
                                             replace = TRUE, 
                                             prob = prob_C1_Dc) + runif(.N)]
  
  # Reset to cured if maximum age was sampled
  m_times[stage_dx == 1 & time_C_Dc > max(ages_C1_Dc), time_C_Dc := Inf]
  
  # Cancer mortality after stage 2 diagnosis
  m_times[stage_dx == 2, time_C_Dc := sample(x = ages_C2_Dc, 
                                             size = .N, 
                                             replace = TRUE, 
                                             prob = prob_C2_Dc) + runif(.N)]
  
  # Reset to cured if maximum age was sampled
  m_times[stage_dx == 2 & time_C_Dc > max(ages_C2_Dc), time_C_Dc := Inf]
  
  # Calculate time to death from cancer
  m_times[, time_H_Dc := time_H_C + time_C_Dc]
}


# Calculate mortality outcomes
calc_mortality_outcomes <- function(m_times) {
  # Calculate age at all-cause death and cause of death
  m_times[, time_H_D := pmin(time_H_Do, time_H_Dc, na.rm = TRUE)]
  m_times[, fl_death_cancer := (time_H_Do > pmin(time_H_Dc, Inf, na.rm = TRUE))]
  
  # Calculate survival from cancer diagnosis
  m_times[time_H_C <= time_H_Dc, time_C_D := time_H_Dc - time_H_C]
}

# Combine parameters and steps to simulate one run of decision model
run_no_screening_model <- function(n_cohort,
                                   ages_H_Do,
                                   prob_H_Do, 
                                   shape_H_P, 
                                   scale_H_P,
                                   rate_P1_P2,
                                   rate_P1_C1,
                                   rate_P2_C2,
                                   ages_C1_Dc,
                                   prob_C1_Dc,
                                   ages_C2_Dc,
                                   prob_C2_Dc,
                                   seed = NULL) {
  
  # Set seed
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  # Initialize matrix of patient data
  m_cohort_base <- initialize_cohort(n_cohort)
  
  # Generate baseline characteristics
  simulate_baseline_data(m_cohort_base,
                         ages_H_Do,
                         prob_H_Do)
  
  # Simulate time to cancer onset
  simulate_cancer_onset(m_cohort_base, 
                        shape_H_P, 
                        scale_H_P)
  
  # Simulate cancer stage progression and diagnosis
  simulate_cancer_progression(m_cohort_base,
                              rate_P1_P2,
                              rate_P1_C1,
                              rate_P2_C2)
  
  # Simulate time from cancer diagnosis to death
  simulate_cancer_mortality(m_cohort_base,
                            ages_C1_Dc,
                            prob_C1_Dc,
                            ages_C2_Dc,
                            prob_C2_Dc)
  
  # Calculate mortality outcomes
  calc_mortality_outcomes(m_cohort_base)
  
  return(m_cohort_base)
}
