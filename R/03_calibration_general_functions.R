# Load general parameters for calibration
load_calib_params <- function(l_params_model, # Model parameters to update
                              target_files, # List of file paths with targets
                              l_outcome_params, # Path to file with outcome parameters
                              l_censor_vars, # List of variables to combine for censor variables
                              prior_file, # Path to parameter mapping as .rds file
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
  l_true_targets <- load_calibration_targets(target_files)
  
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
  prior_map <- readRDS(prior_file)
  
  # Update default parameters with  mean of priors
  l_params_model <- update_param_from_map(l_params_model, 
                                        v_params_update = rowMeans(prior_map[, c("min", "max")]), 
                                        param_map = prior_map)
  
  # Return list of calibration parameters
  l_params_calib <- list(
    l_params_model = l_params_model,
    df_true_targets = df_targets_flattened,
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
calc_calib_outputs <- function(m_patients,
                               l_outcome_params,
                               l_censor_vars = NULL,
                               m_lesions = NULL) {
  # Create censor variables
  if (!is.null(l_censor_vars)) {
    for (varname in names(l_censor_vars)) {
      m_patients[, (varname) := do.call("pmin", c(.SD, na.rm = TRUE)),
                 .SDcols = l_censor_vars[[varname]]]
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
  results <- run_base_model(l_params_update)
  
  # Get calibration targets and reshape to vector
  l_outputs <- calc_calib_outputs(results, 
                                  l_outcome_params = l_outcome_params,
                                  l_censor_vars = l_censor_vars)
  
  if (reshape_output) {
    v_calib_outputs <- reshape_calib_outputs(l_outputs)
    return(v_calib_outputs)
  } else {
    return(v_calib_outputs)
  }
}


calc_and_plot_calib_targets <- function(l_true_targets, 
                                        m_target_outputs,
                                        titles = list(prevalence = 'Prevalence',
                                                      incidence = 'Incidence',
                                                      stage_distr = 'Stage distribution')) {
  # Initialize data
  v_ages <- list()
  l_plot_df <- list()
  v_lengths <- c()
  for (targets in names(l_true_targets)) {
    # Get vector of ages
    if (targets %in% c('prevalence', 'incidence'))
      v_ages[[targets]] <- get_age_range(l_true_targets[[targets]])
    
    # Initialize plotting dataframe
    l_plot_df[[targets]] <- l_true_targets[[targets]] %>%
      mutate(label = 'True')
    
    # Get number of targets
    v_lengths[targets] <- nrow(l_true_targets[[targets]])
  }
  
  # Summarize target model outputs
  model_UB_95 = apply(m_target_outputs, 2, quantile, probs = c(1),  na.rm = TRUE)
  model_LB_95 = apply(m_target_outputs, 2, quantile, probs = c(0),  na.rm = TRUE)
  out_summary = data.frame(targets = rep(names(l_true_targets), v_lengths[names(l_true_targets)]),
                           model_LB_95, 
                           model_UB_95)
  
  # Plot true and comparison targets
  l_plots <- list()
  for (targets in names(l_true_targets)) {
    if (targets %in% c('prevalence', 'incidence')) {
      # Add median age
      l_plot_df[[targets]] <- l_plot_df[[targets]] %>%
        mutate(median_age = (age_end + age_start) / 2) %>%
        mutate(median_age = ifelse(label == 'True', floor(median_age), ceiling(median_age))) %>%
        cbind(out_summary[out_summary$targets == targets,])
      
      if (targets %in% c('incidence')) {
        # Add confidence intervals
        l_plot_df[[targets]] <- l_plot_df[[targets]] %>%
          mutate(ci_lb = incidence - 1.96 * se,
                 ci_ub = incidence + 1.96 * se)
        
        # Plot characteristics
        ylab <- paste('Incidence per', format(as.numeric(l_plot_df[[targets]]$unit[1]), big.mark=","))
      } else {
        # Plot characteristics
        ylab = 'Prevalence'
      }
      
      # Create plot
      l_plots[[targets]] <- ggplot(l_plot_df[[targets]], aes(x = median_age, 
                                                           y = !!as.name(targets))) + 
        geom_point() + 
        geom_errorbar(aes(x=median_age, ymin=ci_lb, ymax=ci_ub), alpha = 0.5) + 
        geom_ribbon(aes(ymin = model_LB_95,
                        ymax = model_UB_95),
                    fill = "black",
                    alpha = 0.3) +
        {if (targets %in% c('incidence')) scale_y_continuous(breaks = seq(0, max(l_plot_df[[targets]]$ci_ub), by = 100))} +
        labs(title = titles[[targets]],
             x = 'Age',
             y = ylab) +
        theme_bw()
    } else if (targets %in% 'stage_distr') {
      # Get model target outputs
      stage_data <- m_target_outputs[, out_summary$targets == 'stage_distr']
      stage_plot_df <- data.frame()
      for (i in 1:ncol(stage_data)) {
        temp_stage_plot_df <- data.frame(
          stage_dx = i,
          vals = stage_data[, i]
        )
        
        stage_plot_df <- rbind(stage_plot_df, temp_stage_plot_df)
      }
      
      # Create plot
      l_plots[[targets]] <- ggplot(l_plot_df[[targets]], aes(x = factor(stage_dx), 
                                                           y = pct)) + 
        geom_point() + 
        geom_errorbar(aes(x=stage_dx, ymin=ci_lb, ymax=ci_ub), alpha = 0.5) + 
        geom_violin(data = stage_plot_df, aes(x = factor(stage_dx), y = vals), fill = 'black', alpha = 0.3) +
        scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
        labs(title = titles[[targets]],
             x = 'Stage',
             y = 'Percentage')+
        theme_bw()
    } else print(paste('Warning: target', targets, 'not found'))
  }
  
  return(l_plots)
}


# Plot calibration targets
plot_calib_targets <- function(l_targets, 
                               titles = list(prevalence = 'Prevalence',
                                             incidence = 'Incidence',
                                             stage_distr = 'Stage distribution'),
                               colorlab = 'Label') {
  # Initialize data
  v_ages <- list()
  l_plot_df <- list()
  v_lengths <- c()
  first_df <- TRUE
  for (label in names(l_targets)) {
    for (targets in names(l_targets[[label]])) {
      if (first_df) {
        # Get vector of ages
        if (targets %in% c('prevalence', 'incidence'))
          v_ages[[targets]] <- get_age_range(l_targets[[label]][[targets]])
        
        # Initialize plotting dataframe
        l_plot_df[[targets]] <- l_targets[[label]][[targets]] %>%
          mutate(label = label)
        
        # Get number of targets
        v_lengths[targets] <- nrow(l_targets[[label]][[targets]])
        
        first_df <- FALSE
      } else {
        # Add to plotting dataframe
        l_plot_df[[targets]] <- rbind(l_plot_df[[targets]],
                                     l_targets[[label]][[targets]] %>%
                                       mutate(label = label))
      }
      
    }
  }
  
  # Plot true and comparison targets
  l_plots <- list()
  for (targets in names(l_plot_df)) {
    if (targets %in% c('prevalence', 'incidence')) {
      # Add median age
      l_plot_df[[targets]] <- l_plot_df[[targets]] %>%
        mutate(median_age = (age_end + age_start) / 2) 
      
      # Find difference in entries and add offset
      offset_scale <- Mode(diff(l_plot_df[[targets]]$median_age))
      # dodge <- position_dodge(width = offset_scale / (length(names(l_targets)) - 1))
      dodge <- position_dodge(width = offset_scale * 0.9)
      
      if (targets %in% c('incidence')) {
        # Add confidence intervals
        l_plot_df[[targets]] <- l_plot_df[[targets]] %>%
          mutate(ci_lb = incidence - 1.96 * se,
                 ci_ub = incidence + 1.96 * se)
        
        # Plot characteristics
        
        ylab <- paste('Incidence per', format(as.numeric(l_plot_df[[targets]]$unit[1]), big.mark=","))
      } else {
        # Plot characteristics
        ylab = 'Prevalence'
      }
      
      # Create plot
      l_plots[[targets]] <- ggplot(l_plot_df[[targets]], aes(x = median_age, 
                                                           y = !!as.name(targets),
                                                           color = label)) + 
        geom_point(position = dodge, alpha = 0.7) + 
        geom_errorbar(aes(x=median_age, ymin=ci_lb, ymax=ci_ub), alpha = 0.5, position = dodge) + 
        {if (targets %in% c('incidence')) scale_y_continuous(breaks = seq(0, max(l_plot_df[[targets]]$ci_ub), by = 100))} +
        labs(title = titles[[targets]],
             x = 'Age',
             y = ylab,
             color = colorlab) +
        theme_bw()
    } else if (targets %in% 'stage_distr') {
      # Find difference in entries and add offset
      offset_scale <- 1
      dodge <- position_dodge(width = offset_scale * 0.9)
      
      # Create plot
      l_plots[[targets]] <- ggplot(l_plot_df[[targets]], aes(x = factor(stage_dx), 
                                                           y = pct,
                                                           color = label)) + 
        geom_point(position = dodge, alpha = 0.7) + 
        geom_errorbar(aes(x=stage_dx, ymin=ci_lb, ymax=ci_ub), alpha = 0.5, position = dodge) + 
        coord_cartesian(ylim = c(0, 1)) +
        scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
        labs(title = titles[[targets]],
             x = 'Stage',
             y = 'Percentage',
             color = colorlab)+
        theme_bw()
    } else print(paste('Warning: target', targets, 'not found'))
  }
  
  return(l_plots)
}
