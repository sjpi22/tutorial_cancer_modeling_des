# Dependencies:
# calc_pcl_prevalence from utils/epi_functions.R

# Recursively read CSV files from list of lists of file paths
recursive_read_csv <- function(l_files) {
  file_data <- list()
  for (item in names(l_files)) {
    if (is.list(l_files[[item]])) {
      file_data[[item]] <- recursive_read_csv(l_files[[item]])
    } else
      file_data[[item]] <- read.csv(l_files[[item]])
  }
  return(file_data)
}

# Update list of model parameters given in the specified parameter vector and mapping
update_calib_params <- function(l_params_all, v_params_update, param_map) {
  # Make copy of parameter list
  l_params_update <- copy(l_params_all)
  
  # Update parameters based on mapping
  assertthat::are_equal(length(v_params_update), nrow(param_map))
  for (i in seq(length(v_params_update))) {
    # Get value to update variable to
    val <- v_params_update[i]
    
    # Get parameter to update
    var_info <- unlist(param_map[i,])
    l_params_update[[var_info['var_name']]]$params[[var_info['param_name']]][as.integer(var_info['param_index'])] <- val
  }
  
  # After resetting all tunable parameters, for empirical distributions, rescale so that probability sums to 1
  for (var in unique(param_map$var_name)) {
    if (l_params_update[[var]]$distr == 'empirical') {
      l_params_update[[var]]$params$probs <- l_params_update[[var]]$params$probs / sum(l_params_update[[var]]$params$probs)
    }
  }
  
  return(l_params_update)
}

# Run model and generate calibration target outputs
calc_calib_targets <- function(l_params_all, m_patients, m_lesions, v_ages_prevalence, v_ages_incidence,
                               n_screen_sample = NULL) {
  
  # Sample patient cohort for screening if 'n_screen_sample" is populated
  if (!is.null(n_screen_sample)) {
    if (n_screen_sample < nrow(m_patients)) {
      m_screen_sample <- m_patients[sample(.N, n_screen_sample)]
    } else {
      m_screen_sample <- m_patients
    }
  } else {
    m_screen_sample <- m_patients
  }
  
  # Create censor variable for precancerous lesion prevalence screening
  m_screen_sample[, time_screen_censor := pmin(time_0_D, time_0_3, na.rm = TRUE)]
  
  # Calculate precancerous lesion prevalence with default model
  prevalence <- calc_pcl_prevalence(l_params_all, m_screen_sample, m_lesions, v_ages_prevalence)
  
  # Calculate cancer age-specific incidence
  incidence <- calc_incidence(m_patients, "time_0_3", "time_0_D", v_ages_incidence) %>%
    dplyr::select(-c("total_atrisk", "age_diff"))
  
  # Calculate cancer stage distribution
  stage_distr <- calc_stage_distr(m_patients, "stage_dx", "time_0_3", "time_0_D")
  
  # Return outputs
  return(list(prevalence = prevalence,
              incidence = incidence,
              stage_distr = stage_distr))
}

# Reshape calibration targets to single vector
# Inputs: list with items list[['prevalence']][['lesiontype']] and list[['incidence']]
reshape_calib_targets <- function(l_calib_targets, output_se = FALSE, conf_level = 0.95) {
  # Initialize vector of outputs
  v_targets <- c()
  if (output_se) {
    v_se <- c()
    z <- qnorm(conf_level + (1-conf_level)/2, 0, 1)
  }
  
  # Prevalence by lesion type
  for (lesiontype in names(l_calib_targets[['prevalence']])) {
    prev_df <- l_calib_targets[['prevalence']][[lesiontype]] %>%
      mutate(varname = paste('prev', lesiontype, age_start, age_end, sep = '_'))
    
    # Extract prevalence and name as varname
    v_prev <- prev_df$prevalence
    names(v_prev) <- prev_df$varname
    
    # Append to target vector
    v_targets = c(v_targets, v_prev)
    
    # Get standard error
    if (output_se) {
      v_temp_se <- (prev_df$ci_ub - prev_df$ci_lb) / (2 * z)
      names(v_temp_se) <- prev_df$varname
      v_se <- c(v_se, v_temp_se)
    }
  }
  
  # Extract incidence and name as varname and append to target vector
  inc_df <- l_calib_targets[['incidence']] %>%
    mutate(varname = paste('inc', 'cancer', age_start, age_end, sep = '_'))
  v_inc <- inc_df$incidence
  names(v_inc) <- inc_df$varname
  v_targets = c(v_targets, v_inc)
  
  # Get standard error
  if (output_se) {
    v_temp_se <- inc_df$se
    names(v_temp_se) <- inc_df$varname
    v_se <- c(v_se, v_temp_se)
  }
  
  # Stage distribution
  stage_df <- l_calib_targets[['stage_distr']] %>%
    mutate(varname = paste('stage', stage_dx, sep = '_'))
  v_stage <- stage_df$pct
  names(v_stage) <- stage_df$varname
  v_targets = c(v_targets, v_stage)
  
  # Get standard error
  if (output_se) {
    v_temp_se <- (stage_df$ci_ub - stage_df$ci_lb) / (2 * z)
    names(v_temp_se) <- stage_df$stage_dx
    v_se <- c(v_se, v_temp_se)
  }
  
  # Return outputs
  if (output_se) {
    return(list(v_targets = v_targets, v_se = v_se))
  } else return(v_targets)
}

# Wrapper for running calibration and outputting vector of outputs
params_to_calib_targets <- function(l_params_all, v_params_update, param_map, v_ages_prevalence, v_ages_incidence) {
  # Update parameters
  l_params_update <- update_calib_params(l_params_all, v_params_update, param_map)
  
  # Run decision model
  results <- run_model(l_params_update)
  results_noscreening <- results[['None']]
  
  # Get calibration targets and reshape to vector
  l_calib_targets <- calc_calib_targets(l_params_all, 
                                        results_noscreening$m_cohort, 
                                        results_noscreening$m_lesions, 
                                        v_ages_prevalence, 
                                        v_ages_incidence)
  v_calib_targets <- reshape_calib_targets(l_calib_targets)
  
  return(v_calib_targets)
}