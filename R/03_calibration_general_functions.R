# Load general parameters for calibration
load_calib_params <- function(l_params_model, # Model parameters to update
                              l_outcome_params, # List of outcome parameters
                              l_censor_vars, # List of variables to combine for censor variables
                              file_priors, # Path to parameter mapping as .rds file
                              outpath = 'output', # Path for output
                              n_cores_reserved_local = 2, # Number of cores to remove when running in parallel locally; all other detected cores will be used for parallel processing
                              n_cohort_calib = NULL, # Number of individuals to simulate for calibration simulations
                              seed_calib = NULL, # Random seed for calibration
                              showWarnings = FALSE
) {
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
  l_targets <- load_calibration_targets(l_outcome_params)
  
  # Further process target data
  v_outcomes_categorical <- c()
  for (label in names(l_outcome_params)) {
    # Subset columns for consolidated target dataframe
    l_targets[[label]] <- l_targets[[label]] %>%
      dplyr::select(any_of(c("target_names", "target_groups", 
                             "target_index",  
                             "age_start", "age_end",
                             "targets", "se", 
                             "ci_lb", "ci_ub", 
                             "total_atrisk", "n_events",
                             "sex", "lesion_type")))
    
    # Subset targets to specified sex if applicable
    if ("sex" %in% colnames(l_targets[[label]])) {
      l_targets[[label]] <- l_targets[[label]] %>%
        filter(sex %in% c(l_params_model$sex, "both"))
    }
    
    # Augment parameters with v_ages for calculating outcomes
    if ("age_start" %in% colnames(l_targets[[label]])) {
      l_outcome_params[[label]]$lit_params$v_ages <- c(l_targets[[label]]$age_start, 
                                                       tail(l_targets[[label]]$age_end, 1))
    }
    
    # Keep list of categorical variables
    if (l_outcome_params[[label]]$categorical == T) v_outcomes_categorical <- c(v_outcomes_categorical, label)
  }
  
  # Flatten targets in order of l_outcome_params
  df_targets_flattened <- rbindlist(l_targets[names(l_outcome_params)], fill = TRUE)
  
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
    v_outcomes_categorical = v_outcomes_categorical,
    l_censor_vars = l_censor_vars,
    prior_map = prior_map,
    outpath = outpath,
    n_cores_reserved_local = n_cores_reserved_local
  )
  
  return(l_params_calib)
}


# Generate parameter sample using Latin hypercube sampling (requires LHS package)
lhs_param_samp <- function(prior_map, # Map of uniform parameter priors with min and max column,
                           n_samp = NULL, # Number of samples to generate
                           n_samp_per_param = NULL, # Number of samples per parameter
                           seed = NULL # Random seed for reproducibility
) { 
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed, kind = "L'Ecuyer-CMRG")
  }
  
  # Get number of parameters to calibrate
  n_param <- nrow(prior_map)
  
  # Set number of samples if needed
  if (is.null(n_samp)) {
    if (!is.null(n_samp_per_param)) {
      n_samp <- n_samp_per_param * n_param
    } else {
      stop("Either n_samp or n_samp_per_param must be provided")
    }
  }
  
  # Sample unit Latin Hypercube
  m_lhs_unit <- randomLHS(n_samp, n_param)
  
  # Rescale to min/max of each parameter
  m_param_samp <- matrix(nrow = n_samp, ncol = n_param)
  for (i in 1:n_param) {
    m_param_samp[, i] <- qunif(
      m_lhs_unit[, i],
      min = prior_map$min[i],
      max = prior_map$max[i])
  }
  colnames(m_param_samp) <- prior_map$var_id
  return(m_param_samp)
}


# Plot coverage (i.e., simulated outputs against targets) and save (requires ggplot2 package)
plot_coverage <- function(
    df_targets, # Dataframe of targets with columns target_index (factor) for x axis, plot_grps for groups of targets to plot in the same graph, targets (target values), and se (standard error)
    m_outputs, # Matrix of outputs with columns corresponding to target_names
    file_fig_coverage = NULL, # File path to save coverage figure
    v_quantiles = c(50, 95), # Inner box and whisker quantiles for boxplot
    plt_size_text = 18 # Size of text in plot
) { 
  # Calculate simulated output mean and box plot quantiles
  df_targets[["model_mean"]] <- colMeans(m_outputs)
  for (i in v_quantiles) {
    df_targets[[paste0("model_LB_", i)]] <- apply(m_outputs, 2, FUN = quantile, probs = (1 - i/100)/2, simplify = TRUE)
    df_targets[[paste0("model_UB_", i)]] <- apply(m_outputs, 2, FUN = quantile, probs = (1 + i/100)/2, simplify = TRUE)
  }
  
  # Convert outputs from wide to long
  out_full_bc_cat <- data.frame(m_outputs) %>%
    pivot_longer(
      cols = everything(), 
      names_to = "target_names",
      values_to = "value"
    ) %>%
    mutate(target_index = rep(df_targets$target_index, nrow(m_outputs)),
           plot_grps = rep(df_targets$plot_grps, nrow(m_outputs)))
  
  # Plot distribution of outputs against targets
  plt <- ggplot(data = df_targets, aes(x = target_index)) + 
    geom_errorbar(
      aes(y    = targets, 
          ymin = targets - se, 
          ymax = targets + se),
      width = 0.3, linewidth = 0.9, color = "red", 
      position = position_nudge(x = -0.2)) +
    theme(legend.position = "none") +
    geom_boxplot(aes(ymin = model_LB_95, 
                     lower = model_LB_50, 
                     middle = model_mean, 
                     upper = model_UB_50, 
                     ymax = model_UB_95),
                 stat = "identity",
                 alpha = 0.5,
                 width = 0.3,
                 position = position_nudge(x = 0.2)) +
    facet_wrap(~plot_grps, scales="free",
               labeller = labeller(plot_grps = label_wrap_gen(plt_size_text*1.5))) +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(), legend.position="none") +
    scale_fill_manual(values = c("grey10", "grey30"))+
    scale_y_continuous(breaks = number_ticks(5))+
    theme_bw(base_size = plt_size_text + 5) +
    theme(plot.title = element_text(size = plt_size_text, face = "bold"),
          axis.text.x = element_text(size = plt_size_text),
          axis.text.y = element_text(size = plt_size_text),
          axis.title = element_text(size = plt_size_text),
          panel.grid.major = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0)) +
    labs(x     = "", y     = "")
  
  # Save plot and adjust size based on number of rows and columns
  if (!is.null(file_fig_coverage)) {
    n_plot_grps <- length(unique(df_targets$plot_grps))
    n_plot_rows <- ceiling(n_plot_grps/3)
    n_plot_cols <- ceiling(n_plot_grps/n_plot_rows)
    ggsave(file_fig_coverage, plot = plt,
           width = n_plot_cols*4, height = 4*n_plot_rows)
  }
  return(plt)
}