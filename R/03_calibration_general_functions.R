# Load general parameters for calibration
load_calib_params <- function(l_params_model, # Model parameters to update
                              l_params_outcome, # List of outcome parameters
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
  l_targets <- load_calibration_targets(l_params_outcome)
  
  # Further process target data
  v_outcomes_categorical <- c()
  for (label in names(l_params_outcome)) {
    # Subset columns for consolidated target dataframe
    l_targets[[label]] <- l_targets[[label]] %>%
      dplyr::select(any_of(c("target_names", "target_groups", 
                             "target_index", "target_index_cat",
                             "age_start", "age_end",
                             "targets", "se", 
                             "ci_lb", "ci_ub", 
                             "n_cases", "n_total",
                             "n_events", "person_years_total",
                             "sex", "lesion_type")))
    
    # Subset targets to specified sex if applicable
    if ("sex" %in% colnames(l_targets[[label]])) {
      l_targets[[label]] <- l_targets[[label]] %>%
        filter(sex %in% c(l_params_model$sex, "both"))
    }
    
    # Augment parameters with v_ages for calculating outcomes
    if ("age_start" %in% colnames(l_targets[[label]])) {
      l_params_outcome[[label]]$lit_params$v_ages <- c(l_targets[[label]]$age_start, 
                                                       tail(l_targets[[label]]$age_end, 1))
    }
    
    # Keep list of categorical variables
    if (l_params_outcome[[label]]$categorical == T) v_outcomes_categorical <- c(v_outcomes_categorical, label)
  }
  
  # Flatten targets in order of l_params_outcome
  df_targets_flattened <- rbindlist(l_targets[names(l_params_outcome)], fill = TRUE)
  
  # Create directory if it does not exist
  outpath_split <- unlist(strsplit(outpath, split = "/"))
  for (i in 1:length(outpath_split)) {
    dir.create(do.call(file.path, as.list(outpath_split[1:i])), showWarnings = showWarnings)
  }
  
  # Load parameter mapping
  prior_map <- read.csv(file_priors)
  prior_map$param_val = rowMeans(prior_map[, c("min", "max")])
  
  # Update stage distribution parameters
  l_params_model[["p_cancer"]] <- l_targets[["stage_distr"]]$targets
  
  # Update default parameters with mean of priors
  l_params_model <- update_param_from_map(l_params_model, 
                                          param_map = prior_map,
                                          update_distr = T)
  
  # Return list of calibration parameters
  l_params_calib <- list(
    l_params_model = l_params_model,
    df_targets = df_targets_flattened,
    l_params_outcome = l_params_outcome,
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


# Plot coverage (i.e., simulated outputs against targets) and save
# Note: requires ggplot2
plot_coverage <- function(
    df_targets, # Dataframe of targets with columns target_index (factor) for x axis, plot_grps for groups of targets to plot in the same graph, targets (target values), se (standard error), and optional categorical indicator variable
    m_outputs = NULL, # Matrix of outputs with columns corresponding to target_names; if NULL, assumes that df_targets already contains summary statistics
    file_fig_coverage = NULL, # File path to save coverage figure
    target_range = "se", # "se" for standard error and "ci" for confidence interval
    v_quantiles = c(50, 95), # Inner box and whisker quantiles for boxplot
    plt_size_text = 18, # Size of text in plot
    labeller_multiplier = 4, # Multiplier for title label text size
    n_cols_max = 3 # Maximum number of columns in plot
) { 
  # Get range to plot for targets
  if (target_range == "se") {
    df_targets <- df_targets %>%
      mutate(diff_lb = se,
             diff_ub = se)
  } else {
    df_targets <- df_targets %>%
      mutate(diff_lb = targets - ci_lb,
             diff_ub = ci_ub - targets)
  }
  
  # Calculate simulated output mean and box plot quantiles
  if (!is.null(m_outputs)) {
    df_targets[["model_mean"]] <- colMeans(m_outputs)
    for (i in v_quantiles) {
      df_targets[[paste0("model_LB_", i)]] <- apply(m_outputs, 2, FUN = quantile, probs = (1 - i/100)/2, simplify = TRUE)
      df_targets[[paste0("model_UB_", i)]] <- apply(m_outputs, 2, FUN = quantile, probs = (1 + i/100)/2, simplify = TRUE)
    }
  }
  
  # Determine whether to plot categorical and/or continuous outputs
  fl_plt <- c()
  l_plts <- list()
  if (!"categorical" %in% colnames(df_targets) | sum(df_targets$categorical, na.rm = T) > 0) {
    # No categorical column or at least one categorical variables - plot categorical variables
    fl_plt["cat"] <- T
    if (!"categorical" %in% colnames(df_targets) | sum(df_targets$categorical, na.rm = T) == nrow(df_targets)) {
      # No categorical column or no continuous variables - default to plotting only categorical variables
      fl_plt["cont"] <- F
      df_targets_cat <- df_targets
    } else {
      # Categorical column and at least one continuous variables - plot both categorical and continuous variables
      fl_plt["cont"] <- T
      df_targets_cat <- df_targets %>%
        filter(categorical == T)
      df_targets_cont <- df_targets %>%
        filter(categorical == F)
    }
    
    # Factorize target index
    if ("target_index_cat" %in% colnames(df_targets_cat)) {
      # Replace missing values of categorical target index values
      v_target_index <- ifelse(is.na(df_targets_cat$target_index_cat), 
                               df_targets_cat$target_index, 
                               df_targets_cat$target_index_cat)
      
      # Replace column
      df_targets_cat$target_index <- v_target_index
    }
    df_targets_cat$target_index <- factor(df_targets_cat$target_index)
    
  } else {
    # No categorical variables in df_targets - plot continuous only
    fl_plt["cat"] <- F
    fl_plt["cont"] <- T
    df_targets_cont <- df_targets %>%
      filter(categorical == F)
  }
  
  # Plot distribution of categorical outputs against targets as box plots
  if (fl_plt["cat"]) {
    l_plts[["cat"]] <- ggplot(data = df_targets_cat, 
                              aes(x = target_index)) + 
      geom_errorbar(aes(y    = targets, 
                        ymin = targets - diff_lb, 
                        ymax = targets + diff_ub),
                    width = 0.3, linewidth = 0.9, color = "red", 
                    position = position_nudge(x = -0.2)) +
      geom_boxplot(aes(ymin = model_LB_95, 
                       lower = model_LB_50, 
                       middle = model_mean, 
                       upper = model_UB_50, 
                       ymax = model_UB_95),
                   stat = "identity",
                   alpha = 0.5,
                   width = 0.3,
                   position = position_nudge(x = 0.2)) +
      facet_wrap(~plot_grps, scales = "free", ncol = n_cols_max,
                 labeller = labeller(plot_grps = label_wrap_gen(plt_size_text*labeller_multiplier/length(unique(df_targets_cat$plot_grps))))) +
      scale_y_continuous(breaks = number_ticks(5)) +
      theme_bw(base_size = plt_size_text + 5) +
      theme(plot.title = element_text(size = plt_size_text, face = "bold"),
            axis.text.x = element_text(size = plt_size_text),
            axis.text.y = element_text(size = plt_size_text),
            axis.title = element_text(size = plt_size_text),
            panel.grid.major = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA),
            strip.background = element_blank(),
            strip.text = element_text(hjust = 0),
            legend.position = "none") +
      labs(x     = "", y     = "")
  }
  
  # Plot distribution of continuous outputs against targets as ribbons
  if (fl_plt["cont"]) {
    l_plts[["cont"]] <- ggplot(data = df_targets_cont, 
                               aes(x = target_index)) + 
      geom_errorbar(aes(y    = targets, 
                        ymin = targets - diff_lb, 
                        ymax = targets + diff_ub),
                    width = 0.3, linewidth = 0.9, color = "red") +
      geom_ribbon(aes(y    = model_mean,
                      ymin = model_LB_95,
                      ymax = model_UB_95),
                  fill = "black",
                  alpha = 0.3) +
      geom_ribbon(aes(y    = model_mean,
                      ymin = model_LB_50,
                      ymax = model_UB_50),
                  fill = "black",
                  alpha = 0.5) +
      facet_wrap(~plot_grps, scales = "free", ncol = n_cols_max,
                 labeller = labeller(plot_grps = label_wrap_gen(plt_size_text*labeller_multiplier/length(unique(df_targets_cont$plot_grps))))) +
      scale_y_continuous(breaks = number_ticks(5)) +
      theme_bw(base_size = plt_size_text + 5) +
      theme(plot.title = element_text(size = plt_size_text, face = "bold"),
            axis.text.x = element_text(size = plt_size_text),
            axis.text.y = element_text(size = plt_size_text),
            axis.title = element_text(size = plt_size_text),
            panel.grid.major = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA),
            strip.background = element_blank(),
            strip.text = element_text(hjust = 0),
            legend.position = "none") +
      labs(x     = "Age", y     = "")
  }
  
  # Depending on number of plots, create final plot layout
  if (length(l_plts) == 1) {
    plt <- l_plts[[1]]
  } else {
    plt <- l_plts[["cont"]] / l_plts[["cat"]]
  }
  
  # Save plot and adjust size based on number of rows and columns
  if (!is.null(file_fig_coverage)) {
    # Get number of rows and columns of graphs
    n_plot_grps <- length(unique(df_targets$plot_grps))
    n_plot_rows <- ceiling(n_plot_grps/n_cols_max)
    n_plot_cols <- ceiling(n_plot_grps/n_plot_rows)
    
    # Save plot
    ggsave(file_fig_coverage, plot = plt,
           width = n_plot_cols*4, height = 4*n_plot_rows)
  }
  return(plt)
}


# Plot correlation of posteriors
plot_posterior_corr <- function(df_post, file_fig_corr = NULL) {
  gg_calib_post_pair_corr <- GGally::ggpairs(
    df_post,
    upper = list(continuous = wrap("cor",
                                   color = "black",
                                   size = 5)),
    diag = list(continuous = wrap("barDiag",
                                  alpha = 0.8)),
    lower = list(continuous = wrap("points",
                                   alpha = 0.3,
                                   size = 0.5)),
    labeller = "label_parsed") +
    theme_bw(base_size = 18) +
    theme(axis.title.x = element_blank(),
          axis.text.x  = element_text(size = 6),
          axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          strip.background = element_rect(fill = "white",
                                          color = "white"),
          strip.text = element_text(hjust = 0))
  
  if (!is.null(file_fig_corr)) {
    ggsave(filename = file_fig_corr,
           gg_calib_post_pair_corr,
           width = 36, height = 24)
  }
  
  return(gg_calib_post_pair_corr)
}
