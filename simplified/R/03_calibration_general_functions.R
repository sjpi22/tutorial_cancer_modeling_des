# Dependencies:
# calc_pcl_prevalence from R/utils/epi_functions.R
# update_param_from_map from R/01_load_model_inputs.R
# run_model from R/02_decision_model_functions.R

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


# Run model and generate calibration target outputs
calc_calib_targets <- function(l_params_all, m_patients, 
                               v_ages,
                               n_screen_sample = NULL,
                               verbose = FALSE) {
  
  # Sample patient cohort for screening if 'n_screen_sample" is populated
  if (!is.null(n_screen_sample)) {
    if (n_screen_sample < nrow(m_patients)) {
      if(verbose) print(paste('Sampling', n_screen_sample, 'patients for lesion prevalence'))
      m_screen_sample <- m_patients[sample(.N, n_screen_sample)]
    } else {
      if(verbose) print(paste('Inputted screening sample size of', n_screen_sample, '< number of patients', nrow(m_patients)))
      m_screen_sample <- m_patients
    }
  } else {
    if(verbose) print('Using full sample for lesion prevalence')
    m_screen_sample <- m_patients
  }
  
  # Create censor variable for precancerous lesion prevalence screening
  m_screen_sample[, time_screen_censor := pmin(time_0_D, time_0_2, na.rm = TRUE)]
  
  # Calculate precancerous lesion prevalence with default model
  prevalence <- calc_prevalence(m_screen_sample, "time_0_1", "time_0_2", "time_screen_censor", 
                                v_ages[['prevalence']],
                                conf_level = l_params_all$conf_level) %>%
    dplyr::select(-c("person_years_cases", "person_years_total"))
  
  # Calculate cancer age-specific incidence
  if(verbose) print('Calculating age-specific incidence')
  incidence <- calc_incidence(m_patients, "time_0_2", "time_0_D", v_ages[['incidence']]) %>%
    dplyr::select(-c("total_atrisk", "age_diff"))
  
  # Calculate cancer stage distribution
  if(verbose) print('Calculating stage distribution')
  stage_distr <- calc_stage_distr(l_params_all, m_patients, "stage_dx", "time_0_2", "time_0_D")
  
  # Return outputs
  return(list(prevalence = prevalence,
              incidence = incidence,
              stage_distr = stage_distr))
}

# Reshape calibration targets to single vector
reshape_calib_targets <- function(l_calib_targets, output_se = FALSE, 
                                  output_map = FALSE,
                                  conf_level = 0.95) {
  # Initialize vector(s) of outputs
  v_targets <- c()
  if (output_se) {
    v_se <- c()
    z <- qnorm(conf_level + (1-conf_level)/2, 0, 1)
  }
  if (output_map) {
    target_map <- data.frame()
  }
  
  # Prevalence by lesion type
  prev_df <- l_calib_targets[['prevalence']] %>%
    mutate(varname = paste('prev', age_start, age_end, sep = '_'))
  
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
  
  # Get output mapping
  if (output_map) {
    target_map <- rbind(target_map,
                        data.frame(target_groups = paste("Prevalence of precancerous lesion"),
                                   target_index = (prev_df$age_start + prev_df$age_end)/2))
    
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
  
  # Get output mapping
  if (output_map) {
    target_map <- rbind(target_map,
                    data.frame(target_groups = "Cancer incidence",
                               target_index = (inc_df$age_start + inc_df$age_end)/2))
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
  
  # Get output mapping
  if (output_map) {
    target_map <- rbind(target_map,
                    data.frame(target_groups = "Stage at diagnosis",
                               target_index = stage_df$stage_dx))
  }
  
  # Return outputs
  if (output_se | output_map) {
    v_outputs <- list(v_targets = v_targets)
    if (output_se) {
      v_outputs$v_se <- v_se
    } 
    if (output_map) {
      v_outputs$target_map <- target_map
    }
    return(v_outputs)
  } else return(v_targets)
}

# Wrapper for running calibration and outputting vector of outputs
params_to_calib_targets <- function(l_params_all, v_params_update, param_map, 
                                    v_ages, reshape_output = TRUE,
                                    output_se = FALSE, output_map = FALSE, 
                                    conf_level = 0.95,
                                    verbose = FALSE) {
  # Update parameters
  l_params_update <- update_param_from_map(l_params_all, v_params_update, param_map)
  
  # Run decision model
  results <- run_model(l_params_update, verbose = verbose)
  results_noscreening <- results[['None']]
  
  # Get calibration targets and reshape to vector
  l_calib_targets <- calc_calib_targets(l_params_update, 
                                        results_noscreening,
                                        v_ages,
                                        verbose = verbose)
  if (reshape_output) {
    v_calib_targets <- reshape_calib_targets(l_calib_targets, output_se = output_se,
                                             output_map = output_map, conf_level = 0.95)
    return(v_calib_targets)
  } else {
    return(l_calib_targets)
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
  for (target in names(l_true_targets)) {
    # Get vector of ages
    if (target %in% c('prevalence', 'incidence'))
      v_ages[[target]] <- get_age_range(l_true_targets[[target]])
    
    # Initialize plotting dataframe
    l_plot_df[[target]] <- l_true_targets[[target]] %>%
      mutate(label = 'True')
    
    # Get number of targets
    v_lengths[target] <- nrow(l_true_targets[[target]])
  }
                
  # Summarize target model outputs
  model_UB_95 = apply(m_target_outputs, 2, quantile, probs = c(1),  na.rm = TRUE)
  model_LB_95 = apply(m_target_outputs, 2, quantile, probs = c(0),  na.rm = TRUE)
  out_summary = data.frame(target = rep(names(l_true_targets), v_lengths[names(l_true_targets)]),
                           model_LB_95, 
                           model_UB_95)
  
  # Plot true and comparison targets
  l_plots <- list()
  for (target in names(l_true_targets)) {
    if (target %in% c('prevalence', 'incidence')) {
      # Add median age
      l_plot_df[[target]] <- l_plot_df[[target]] %>%
        mutate(median_age = (age_end + age_start) / 2) %>%
        mutate(median_age = ifelse(label == 'True', floor(median_age), ceiling(median_age))) %>%
        cbind(out_summary[out_summary$target == target,])
      
      if (target %in% c('incidence')) {
        # Add confidence intervals
        l_plot_df[[target]] <- l_plot_df[[target]] %>%
          mutate(ci_lb = incidence - 1.96 * se,
                 ci_ub = incidence + 1.96 * se)
        
        # Plot characteristics
        ylab <- paste('Incidence per', format(as.numeric(l_plot_df[[target]]$unit[1]), big.mark=","))
      } else {
        # Plot characteristics
        ylab = 'Prevalence'
      }
      
      # Create plot
      l_plots[[target]] <- ggplot(l_plot_df[[target]], aes(x = median_age, 
                                                           y = !!as.name(target))) + 
        geom_point() + 
        geom_errorbar(aes(x=median_age, ymin=ci_lb, ymax=ci_ub), alpha = 0.5) + 
        geom_ribbon(aes(ymin = model_LB_95,
                        ymax = model_UB_95),
                    fill = "black",
                    alpha = 0.3) +
        {if (target %in% c('incidence')) scale_y_continuous(breaks = seq(0, max(l_plot_df[[target]]$ci_ub), by = 100))} +
        labs(title = titles[[target]],
             x = 'Age',
             y = ylab) +
        theme_bw()
    } else if (target %in% 'stage_distr') {
      # Get model target outputs
      stage_data <- m_target_outputs[, out_summary$target == 'stage_distr']
      stage_plot_df <- data.frame()
      for (i in 1:ncol(stage_data)) {
        temp_stage_plot_df <- data.frame(
          stage_dx = i,
          vals = stage_data[, i]
        )
        
        stage_plot_df <- rbind(stage_plot_df, temp_stage_plot_df)
      }
      
      # Create plot
      l_plots[[target]] <- ggplot(l_plot_df[[target]], aes(x = factor(stage_dx), 
                                                           y = pct)) + 
        geom_point() + 
        geom_errorbar(aes(x=stage_dx, ymin=ci_lb, ymax=ci_ub), alpha = 0.5) + 
        geom_violin(data = stage_plot_df, aes(x = factor(stage_dx), y = vals), fill = 'black', alpha = 0.3) +
        scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
        labs(title = titles[[target]],
             x = 'Stage',
             y = 'Percentage')+
        theme_bw()
    } else print(paste('Warning: target', target, 'not found'))
  }
  
  return(l_plots)
}

# Mode from https://stackoverflow.com/questions/2547402/how-to-find-the-statistical-mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
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
    for (target in names(l_targets[[label]])) {
      if (first_df) {
        # Get vector of ages
        if (target %in% c('prevalence', 'incidence'))
          v_ages[[target]] <- get_age_range(l_targets[[label]][[target]])
        
        # Initialize plotting dataframe
        l_plot_df[[target]] <- l_targets[[label]][[target]] %>%
          mutate(label = label)
        
        # Get number of targets
        v_lengths[target] <- nrow(l_targets[[label]][[target]])
        
        first_df <- FALSE
      } else {
        # Add to plotting dataframe
        l_plot_df[[target]] <- rbind(l_plot_df[[target]],
                                     l_targets[[label]][[target]] %>%
                                       mutate(label = label))
      }
      
    }
  }
  
  # Plot true and comparison targets
  l_plots <- list()
  for (target in names(l_plot_df)) {
    if (target %in% c('prevalence', 'incidence')) {
      # Add median age
      l_plot_df[[target]] <- l_plot_df[[target]] %>%
        mutate(median_age = (age_end + age_start) / 2) 
      
      # Find difference in entries and add offset
      offset_scale <- Mode(diff(l_plot_df[[target]]$median_age))
      # dodge <- position_dodge(width = offset_scale / (length(names(l_targets)) - 1))
      dodge <- position_dodge(width = offset_scale * 0.9)
      
      if (target %in% c('incidence')) {
        # Add confidence intervals
        l_plot_df[[target]] <- l_plot_df[[target]] %>%
          mutate(ci_lb = incidence - 1.96 * se,
                 ci_ub = incidence + 1.96 * se)
        
        # Plot characteristics
        
        ylab <- paste('Incidence per', format(as.numeric(l_plot_df[[target]]$unit[1]), big.mark=","))
      } else {
        # Plot characteristics
        ylab = 'Prevalence'
      }
      
      # Create plot
      l_plots[[target]] <- ggplot(l_plot_df[[target]], aes(x = median_age, 
                                                           y = !!as.name(target),
                                                           color = label)) + 
        geom_point(position = dodge, alpha = 0.7) + 
        geom_errorbar(aes(x=median_age, ymin=ci_lb, ymax=ci_ub), alpha = 0.5, position = dodge) + 
        {if (target %in% c('incidence')) scale_y_continuous(breaks = seq(0, max(l_plot_df[[target]]$ci_ub), by = 100))} +
        labs(title = titles[[target]],
             x = 'Age',
             y = ylab,
             color = colorlab) +
        theme_bw()
    } else if (target %in% 'stage_distr') {
      # Find difference in entries and add offset
      offset_scale <- 1
      dodge <- position_dodge(width = offset_scale * 0.9)
      
      # Create plot
      l_plots[[target]] <- ggplot(l_plot_df[[target]], aes(x = factor(stage_dx), 
                                                           y = pct,
                                                           color = label)) + 
        geom_point(position = dodge, alpha = 0.7) + 
        geom_errorbar(aes(x=stage_dx, ymin=ci_lb, ymax=ci_ub), alpha = 0.5, position = dodge) + 
        coord_cartesian(ylim = c(0, 1)) +
        scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
        labs(title = titles[[target]],
             x = 'Stage',
             y = 'Percentage',
             color = colorlab)+
        theme_bw()
    } else print(paste('Warning: target', target, 'not found'))
  }
  
  return(l_plots)
}
