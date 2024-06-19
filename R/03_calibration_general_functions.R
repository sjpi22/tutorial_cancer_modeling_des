# Dependencies:
# calc_pcl_prevalence from utils/epi_functions.R

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
calc_calib_targets <- function(l_params_all, m_patients, m_lesions, v_ages_prevalence, v_ages_incidence) {
  
  # Create censor variable for precancerous lesion prevalence screening
  m_patients[, time_screen_censor := pmin(time_0_D, time_0_3, na.rm = TRUE)]
  
  # Calculate precancerous lesion prevalence with default model
  model_prevalence <- calc_pcl_prevalence(l_params_all, m_patients, m_lesions, v_ages_prevalence)
  
  # Calculate cancer incidence
  model_incidence <- calc_incidence(m_patients, "time_0_3", "time_0_D", v_ages_incidence)
  
  # Return outputs
  return(list(model_prevalence = model_prevalence,
              model_incidence = model_incidence))
}

# Reshape calibration targets to single vector
reshape_calib_targets <- function(l_calib_targets) {
  # Initialize vector of outputs
  v_targets <- c()
  
  # Prevalence by lesion type
  for (lesiontype in names(l_calib_targets$model_prevalence)) {
    prev_df <- l_calib_targets$model_prevalence[[lesiontype]] %>%
      mutate(varname = paste('prev', lesiontype, age_start, age_stop, sep = '_'))
    
    # Extract prevalence and name as varname
    v_prev <- prev_df$prevalence
    names(v_prev) <- prev_df$varname
    
    # Append to target vector
    v_targets = c(v_targets, v_prev)
  }
  
  # Incidence
  inc_df <- l_calib_targets$model_incidence %>%
    mutate(varname = paste('inc', 'cancer', age_min, age_max, sep = '_'))
  
  # Extract incidence and name as varname
  v_inc <- inc_df$incidence
  names(v_inc) <- inc_df$varname
  
  # Append to target vector
  v_targets = c(v_targets, v_inc)
  
  return(v_targets)
}

# Wrapper for running calibration and outputting vector of outputs
params_to_calib_targets <- function(l_params_all, v_params_update, param_map, v_ages_prevalence, v_ages_incidence) {
  # Update parameters
  l_params_update <- update_calib_params(l_params_all, v_params_update, param_map)
  
  # Run decision model
  results <- run_model(l_params_update)
  results_noscreening <- results[['None']]
  
  # Get calibration targets and reshape to vector
  l_calib_targets <- calc_calib_targets(l_params_all, results_noscreening$m_cohort, results_noscreening$m_lesions, v_ages_prevalence, v_ages_incidence)
  v_calib_targets <- reshape_calib_targets(l_calib_targets)
  
  return(v_calib_targets)
}