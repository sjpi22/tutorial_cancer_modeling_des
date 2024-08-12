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

# Load general parameters for calibration
load_calib_params <- function(l_params_all, # Model parameters to update
                              target_files, # List of file paths with targets
                              prior_file = 'data/priors.rds', # Path to parameter mapping as .rds file
                              outpath = 'output', # Path for output
                              n_cores_reserved_local = 2, # Number of cores to remove when running in parallel locally; all other detected cores will be used for parallel processing
                              n_cohort_calib = NULL, # Number of individuals to simulate for calibration simulations
                              seed_calib = NULL, # Random seed for calibration
                              showWarnings = FALSE) {
  
  # Update model parameters
  l_params_all <- update_param_list(l_params_all,
                                    list(v_strats = l_params_all$v_strats[1]))
  
  if (!is.null(n_cohort_calib)) {
    l_params_all <- update_param_list(l_params_all,
                                      list(n_cohort = n_cohort_calib))
  }
  
  if (!is.null(seed_calib)) {
    l_params_all <- update_param_list(l_params_all,
                                      list(seed = seed_calib))
  }
  
  # Load targets
  l_true_targets <- recursive_read_csv(target_files)
  
  # Get vector of ages for prevalence and incidence
  v_ages_prevalence <- get_age_range(l_true_targets$prevalence)
  v_ages_incidence <- get_age_range(l_true_targets$incidence)
  v_ages <- list(prevalence = v_ages_prevalence,
                 incidence = v_ages_incidence)
  
  # Parameters for calculating targets
  # list(prevalence = list(start_var,
  #                     end_var,
  #                     censor_var,
  #                     v_ages,
  #                     conf_level),
  #      incidence = list(time_var, censor_var, v_ages, strat_var, rate_unit),
  #      stage_distr = list(groups_expected, grouping_var, event_var, censor_var, conf_level))
  
  # Create directory if it does not exist
  outpath_split <- unlist(strsplit(outpath, split = "/"))
  for (i in 1:length(outpath_split)) {
    dir.create(do.call(file.path, as.list(outpath_split[1:i])), showWarnings = showWarnings)
  }
  
  # Load parameter mapping
  prior_map <- readRDS(prior_file)
  
  # Return list of calibration parameters
  l_params_calib <- list(l_params_all = l_params_all,
                         l_true_targets = l_true_targets,
                         v_ages = v_ages,
                         prior_map = prior_map,
                         outpath = outpath,
                         n_cores_reserved_local = n_cores_reserved_local)
  
}

add_IMABC_params <- function(l_params_calib,
                             alpha_current = 0.0001,
                             alpha_stop = 0.05,
                             use_seed = FALSE,
                             N_start_multiplier = 1000,
                             optional_args = list()) {
  
  # Define Model Parameters/Priors
  param_df <- with(l_params_calib, {
    data.frame(
      name_var = prior_map$var_id,
      dist_var = prior_map$distr,
      min = prior_map$min,
      max = prior_map$max
    )
  })
  
  priors <- as.priors(
    param_df,
    parameter_name = "name_var", dist_base_name = "dist_var"
  )
  
  # Reshape true targets, get mapping, and calculate confidence intervals for true targets given alpha values
  l_true_reshaped <- reshape_calib_targets(l_params_calib$l_true_targets,
                                           output_map = TRUE)
  
  # Calculate confidence intervals
  ci <- calc_target_ci(l_params_calib$l_true_targets,
                       alpha = list(current = alpha_current, stopping = alpha_stop),
                       reshape_outputs = TRUE)
  
  # Define Target Values
  target_df <- with(l_true_reshaped, {
    
    target_df <- data.frame(
      target_groups = make.names(target_map$target_groups),
      target_names = names(v_targets),
      targets = unname(v_targets),
      current_lower_bounds = ci$current$lb,
      current_upper_bounds = ci$current$ub,
      stopping_lower_bounds = ci$stopping$lb,
      stopping_upper_bounds = ci$stopping$ub
    )
    
    return(target_df)
  })
  
  targets <- as.targets(target_df)
  
  # Define Target Function
  fn <- function(v_params_update) {
    v_targets <- with(l_params_calib, {
      params_to_calib_targets(l_params_all, v_params_update, prior_map, 
                              v_ages)
    })
    return(v_targets)
  }
  
  target_fun <- define_target_function(
    targets, priors, FUN = fn, use_seed = use_seed
  )
  
  # Assign N_start
  N_start = N_start_multiplier * length(priors)
  
  # IMABC core parameters
  imabc_inputs <- list(
    target_fun = target_fun,
    priors = priors,
    targets = targets,
    N_start = N_start,
    seed = l_params_calib$l_params_all$seed)
  
  imabc_inputs <- c(imabc_inputs,
                    optional_args)
  
  # Return augmented calibration parameter list
  l_params_calib <- c(l_params_calib, 
                      fn = fn,
                      use_seed = use_seed,
                      imabc_inputs = list(imabc_inputs))
  
  return(l_params_calib)
}

# Calculate 1-alpha confidence intervals given alpha as a single value, vector, or list
calc_target_ci <- function(l_calib_targets,
                           alpha = 0.05,
                           reshape_outputs = TRUE) {
  # Check if alpha is a list
  if (is.list(alpha)) {
    ci <- list()
    for (label in names(alpha)) {
      ci[[label]] <- calc_target_ci_internal(l_calib_targets, alpha[[label]], reshape_outputs = TRUE)
    }
  } else {
    ci <- calc_target_ci_internal(l_calib_targets, alpha, reshape_outputs = TRUE)
  }
  
  return(ci)
}

# Calculate confidence intervals of calibration targets given alpha as a single value or vector
calc_target_ci_internal <- function(l_calib_targets,
                                    alpha = 0.05, 
                                    reshape_outputs = TRUE) {
  
  # If alpha is a single value, turn into vector
  if (is.numeric(alpha) & (length(alpha) == 1)) {
    alpha <- rep(alpha, length(l_calib_targets))
  } else {
    # Check that alpha is a vector of length equal to number of targets
    assert_that(is.numeric(alpha) & (length(alpha) == length(l_calib_targets)),
                msg = "alpha should be a single value or vector with length equal to number of targets.")
  }
  
  # Name vector by targets
  names(alpha) <- names(l_calib_targets)
  
  # Loop over targets
  ci <- list()
  for (target in names(l_calib_targets)) {
    # Define parameters
    if (target == 'prevalence') {
      # Get binomial exact confidence interval with 'binom' package
      df_confint <- with(l_calib_targets[[target]], {
        binom.confint(x = prevalence * n_total, n = n_total, conf.level = 1-alpha[target], methods = "exact")
      })
      
      # Extract CI as list
      ci[[target]] <- list(lb = df_confint$lower,
                           ub = df_confint$upper)
    } else if (target == 'stage_distr') {
      # Get binomial exact confidence interval with 'binom' package
      df_confint <- with(l_calib_targets[[target]], {
        binom.confint(x = count, n = sum(count), conf.level = 1-alpha, methods = "exact")
      })
      
      # Extract CI as list
      ci[[target]] <- list(lb = df_confint$lower,
                           ub = df_confint$upper)
    } else if (target == 'incidence') {
      # Calculate CI with normal approximation
      ci[[target]] <- with(l_calib_targets[[target]], {
        # Width of confidence interval
        ci_width <- qnorm((1-alpha) + alpha/2) * se
        
        # Create CI
        ci <- list()
        ci$lb <- pmax(0, incidence - ci_width) # Lower bound
        ci$ub <- incidence + ci_width # Upper bound
        
        return(ci)
      })
      
    }
  }
  
  # If reshape_outputs is TRUE, flatten target outputs
  if (reshape_outputs) {
    lb <- ub <- c()
    for (target in ci) {
      lb <- c(lb, target$lb)
      ub <- c(ub, target$ub)
    }
    return(list(lb = lb,
                ub = ub))
  } else return (ci)
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
  m_screen_sample[, time_screen_censor := pmin(time_H_D, time_H_C, na.rm = TRUE)]
  
  # Calculate precancerous lesion prevalence with default model
  prevalence <- calc_prevalence(m_screen_sample, "time_H_P", "time_H_C", "time_screen_censor", 
                                v_ages[['prevalence']],
                                conf_level = l_params_all$conf_level) %>%
    dplyr::select(-c("person_years_cases", "person_years_total"))
  
  # Calculate cancer age-specific incidence
  if(verbose) print('Calculating age-specific incidence')
  incidence <- calc_incidence(m_patients, "time_H_C", "time_H_D", v_ages[['incidence']]) %>%
    dplyr::select(-c("total_atrisk", "age_diff"))
  
  # Calculate cancer stage distribution
  if(verbose) print('Calculating stage distribution')
  stage_distr <- calc_stage_distr(l_params_all, m_patients, "stage_dx", "time_H_C", "time_H_D")
  
  # Return outputs
  return(list(prevalence = prevalence,
              incidence = incidence,
              stage_distr = stage_distr))
}

# Reshape calibration targets to single vector
reshape_calib_targets <- function(l_calib_targets, 
                                  output_map = FALSE) {
  # Initialize vector(s) of outputs
  v_targets <- c()
  if (output_map) {
    target_map <- data.frame()
  }
  
  # Loop over target list
  for (target in names(l_calib_targets)) {
    # Age-specific targets
    if (grepl("^(prevalence|incidence)", target)) {
      target_df <- l_calib_targets[[target]] %>%
        mutate(varname = paste(target, age_start, age_end, sep = '_'))
      
      # Select vector of targets
      add_v_target <- target_df[[target]]
      
      # Get output mapping
      if (output_map) {
        target_map <- rbind(target_map,
                            data.frame(target_groups = target,
                                       target_index = (target_df[['age_start']] + target_df[['age_end']])/2))
      } 
    } else { # Non age-specific targets
      # Stage distribution
      target_df <- l_calib_targets[[target]] %>%
        mutate(varname = paste(target, stage_dx, sep = '_'))
      
      # Select vector of targets
      add_v_target <- target_df[['pct']]
      
      # Get output mapping
      if (output_map) {
        target_map <- rbind(target_map,
                            data.frame(target_groups = target,
                                       target_index = target_df$stage_dx))
      }
    }
    
    # Append to target vector
    names(add_v_target) <- target_df$varname
    v_targets = c(v_targets, add_v_target)
  }
  
  # Return outputs
  if (output_map) {
    v_outputs <- list(v_targets = v_targets)
    if (output_map) {
      v_outputs$target_map <- target_map
    }
    return(v_outputs)
  } else return(v_targets)
}


# Wrapper for running calibration and outputting vector of outputs
params_to_calib_targets <- function(l_params_all, v_params_update, param_map, 
                                    v_ages, reshape_output = TRUE, 
                                    output_map = FALSE, 
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
    v_calib_targets <- reshape_calib_targets(l_calib_targets,
                                             output_map = output_map)
    return(v_calib_targets)
  } else {
    return(l_calib_targets)
  }
  
}

# Simulate model outputs for parameter samples
# Note: If running in parallel (run_parallel = TRUE), requires packages 
# doParallel and foreach, and requires registering cores using
# registerDoParallel()
param_sample_to_outputs <- function(m_param_samp, 
                                    fn = 'params_to_calib_targets', 
                                    param_arg_name = 'v_params_update',
                                    fn_other_args = list(),
                                    run_parallel = TRUE) {
  
  # Define function arguments
  fn_args <- list('.')
  names(fn_args) <- param_arg_name
  fn_args <- c(fn_args, fn_other_args)
  
  # Simulate with parallel processing
  if(run_parallel) {
    stime <- system.time({
      m_outputs <- foreach(
        i=1:nrow(m_param_samp), 
        .combine=rbind, 
        .inorder=FALSE, 
        .packages=c("data.table","tidyverse")) %dopar% {
          # Get row of parameters and calculate targets
          v_params_update <- m_param_samp[i,]
          fn_args[[param_arg_name]] <- v_params_update
          v_calib_targets <- do.call(what = fn,
                                     args = fn_args)
          # Call item to save
          t(v_calib_targets)
        }
    })
    print(stime)
    closeAllConnections()
  } else { # Simulate with apply()
    
    fn <- function(v_params_update) {
      v_outputs <- do.call(what = fn,
                           args = fn_args,
                           quote = TRUE)
      return(v_outputs)
    }
    
    start_time <- Sys.time()
    m_outputs <- apply(X, MARGIN = 1, FUN = fn)
    end_time <- Sys.time()
    print(paste('Simulation time:', end_time - start_time))
  }
  return(m_outputs)
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
