# Load general parameters for calibration
load_calib_params <- function(l_params_model, # Model parameters to update
                              l_outcome_params, # List of outcome parameters
                              l_censor_vars, # List of variables to combine for censor variables
                              file_priors, # Path to parameter mapping as .rds file
                              outpath = 'output', # Path for output
                              n_cores_reserved_local = 2, # Number of cores to remove when running in parallel locally; all other detected cores will be used for parallel processing
                              n_cohort_calib = NULL, # Number of individuals to simulate for calibration simulations
                              seed_calib = NULL, # Random seed for calibration
                              showWarnings = FALSE) {
  
  # Update model parameters from function inputs
  updated_params <- list(v_strats = l_params_model$v_strats[1])
  if (!is.null(n_cohort_calib)) {
    updated_params <- c(updated_params, list(n_cohort = n_cohort_calib))
  }
  
  if (!is.null(seed_calib)) {
    updated_params <- c(updated_params, list(seed = seed_calib))
  }
  
  l_params_model <- update_param_list(
    l_params_model,
    updated_params
  )
  
  # Load targets
  l_true_targets <- load_calibration_targets(l_outcome_params)
  
  # Further process target data
  for (label in names(l_outcome_params)) {
    # Subset columns for consolidated target dataframe
    l_true_targets[[label]] <- l_true_targets[[label]] %>%
      dplyr::select(any_of(c("target_names", "target_groups", 
                             "target_index",  
                             "age_start", "age_end",
                             "targets", "se", 
                             "ci_lb", "ci_ub", 
                             "total_atrisk", "n_events",
                             "sex", "lesion_type")))
    
    # Subset targets to specified sex if applicable
    if ("sex" %in% colnames(l_true_targets[[label]])) {
      l_true_targets[[label]] <- l_true_targets[[label]] %>%
        filter(sex %in% c(l_params_model$sex, "both"))
    }
    
    # Augment parameters with v_ages for calculating outcomes
    if ("age_start" %in% colnames(l_true_targets[[label]])) {
      l_outcome_params[[label]]$v_ages <- c(l_true_targets[[label]]$age_start, 
                                            tail(l_true_targets[[label]]$age_end, 1))
    }
  }
  
  # Flatten targets in order of l_outcome_params
  df_targets_flattened <- rbindlist(l_true_targets[names(l_outcome_params)], fill = TRUE)
  
  # Create directory if it does not exist
  outpath_split <- unlist(strsplit(outpath, split = "/"))
  for (i in 1:length(outpath_split)) {
    dir.create(do.call(file.path, as.list(outpath_split[1:i])), showWarnings = showWarnings)
  }
  
  # Load parameter mapping
  prior_map <- read.csv(file_priors)
  prior_map$param_val = rowMeans(prior_map[, c("min", "max")])
  
  # Update default parameters with  mean of priors
  l_params_model <- update_param_from_map(l_params_model, 
                                          param_map = prior_map,
                                          update_distr = T)
  
  # Return list of calibration parameters
  l_params_calib <- list(
    l_params_model = l_params_model,
    df_targets = df_targets_flattened,
    l_outcome_params = l_outcome_params,
    l_censor_vars = l_censor_vars,
    prior_map = prior_map,
    outpath = outpath,
    n_cores_reserved_local = n_cores_reserved_local
  )
  
  return(l_params_calib)
}


# Run model and generate calibration target outputs
# l_censor_vars: list of lists of variables to use to create new variables for censoring outcomes
calc_calib_outputs <- function(m_cohort,
                               l_outcome_params,
                               l_censor_vars = NULL) {
  
  # Separate patient and lesion data as necessary
  if (is.data.table(m_cohort)) {
    m_patients <- m_cohort
  } else {
    m_patients <- m_cohort$patient_level
    m_lesions <- m_cohort$lesion_level
  }
  
  # Create censor variables
  if (!is.null(l_censor_vars)) {
    for (dt in names(l_censor_vars)) {
      for (varname in names(l_censor_vars[[dt]])) {
        get(dt)[, (varname) := do.call("pmin", c(.SD, na.rm = TRUE)),
                .SDcols = l_censor_vars[[dt]][[varname]]]
      }
    }
  }
  
  # Calculate outcomes
  l_results <- list()
  for (outcome in names(l_outcome_params)) {
    l_results[[outcome]] <- do.call(
      paste0("calc_", l_outcome_params[[outcome]][[2]]), 
      c(list(get(l_outcome_params[[outcome]][[3]])), tail(l_outcome_params[[outcome]], -3)))
  }
  
  # Return outputs
  return(l_results)
}


# Reshape calibration outputs to single vector
reshape_calib_outputs <- function(l_calib_outputs, var_to_keep = "value") {
  v_outputs <- c()
  for (df in l_calib_outputs) {
    if (var_to_keep %in% names(df)) {
      v_outputs <- c(v_outputs, df[[var_to_keep]])
    } else {
      if (is.null(dim(df))) {
        n_vals <- 1
      } else {
        n_vals <- nrow(df)
      }
      v_outputs <- c(v_outputs, rep(NA, n_vals))
    }
  }
  
  return(v_outputs)
}


# Wrapper for running calibration and outputting vector of outputs
params_to_calib_outputs <- function(l_params_model, v_params_update = NULL, 
                                    param_map = NULL,
                                    l_outcome_params, l_censor_vars,
                                    reshape_output = TRUE, 
                                    output_map = FALSE, 
                                    conf_level = 0.95) {
  # Update parameters
  if (!is.null(v_params_update)) {
    if (is.null(param_map)) {
      stop("Input parameter map")
    }
    l_params_update <- update_param_from_map(l_params_model, v_params_update, param_map)
  } else {
    l_params_update <- l_params_model
  }
  
  # Run decision model
  m_cohort <- run_base_model(l_params_update)
  
  # Get calibration targets and reshape to vector
  l_outputs <- calc_calib_outputs(m_cohort, 
                                  l_outcome_params = l_outcome_params,
                                  l_censor_vars = l_censor_vars)
  
  if (reshape_output) {
    v_calib_outputs <- reshape_calib_outputs(l_outputs)
    return(v_calib_outputs)
  } else {
    return(v_calib_outputs)
  }
}

