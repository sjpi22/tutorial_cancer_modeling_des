###########################  Load Calibration Parameters  ##########################
#
#  Objective: Program to load general calibration parameters, perform Monte 
#  Carlo error analysis for sample size, and check coverage of targets
########################### <<<<<>>>>> #########################################

rm(list = ls()) # Clean environment
options(scipen = 999) # View data without scientific notation

#### 1.Libraries and functions  ==================================================

###### 1.1 Load packages
library(readxl)
library(data.table)
library(tidyverse)
library(lhs)
library(doParallel)
library(foreach)
library(assertthat)

###### 1.2 Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)


#### 2. General parameters ========================================================

###### 2.1 Configurations
# Run file to process configurations
source("configs/process_configs.R")

# Extract relevant parameters from configs
params_model <- configs$params_model
params_calib <- configs$params_calib
file_params_calib <- configs$paths$file_params_calib
file_plot_labels <- configs$paths$file_plot_labels

# Get list of relevant output file paths and load to global environment
l_filepaths <- update_config_paths("files_coverage", configs$paths)
list2env(l_filepaths, envir = .GlobalEnv)

# Load Monte Carlo and coverage analysis parameters from configs file
list2env(configs$params_montecarlo, envir = .GlobalEnv)
list2env(configs$params_coverage, envir = .GlobalEnv)


#### 3. Pre-processing  ===========================================

# Load model parameters
l_params_init <- do.call(load_model_params, c(
  params_model
))

# Load calibration parameters
l_params_calib <- do.call(load_calib_params, c(
  l_params_model = list(l_params_init),
  params_calib
))

# Load plot labels
df_plot_labels <- read.csv(file_plot_labels)


#### 4. Monte Carlo error analysis  ===========================================

# Set Monte Carlo simulation parameters
l_params_mc <- l_params_calib$l_params_model
l_params_mc$n_cohort <- n_init
l_params_mc$seed <- NULL # Remove seed that is reset every time model is run
set.seed(l_params_calib$l_params_model$seed) # Set seed externally

df_res_mc <- matrix(nrow = n_mc_reps, ncol = nrow(l_params_calib$df_targets))
for (i in 1:n_mc_reps) {
  df_res_mc[i, ] <- params_to_calib_outputs(l_params_mc, 
                                            l_outcome_params = l_params_calib$l_outcome_params,
                                            l_censor_vars = l_params_calib$l_censor_vars)
}

# Get column-wise mean and SD
mean_mc <- colMeans(df_res_mc)
sd_mc <- apply(df_res_mc, 2, sd)

# Calculate required sample size
n_target <- n_init * (sd_mc / l_params_calib$df_targets$se)^2

# Set final N as multiple of maximum required sample size, rounded up to second highest digit
n_final <- mc_multiplier * max(n_target)
n_digits <- floor(log(n_final, base = 10))
n_final <- ceiling(n_final / 10^(n_digits - 1)) * 10^(n_digits - 1)

# Update parameters
l_params_calib$l_params_model$n_cohort <- n_final


#### 5. Save parameters  ===========================================

saveRDS(l_params_calib, file = file_params_calib)


#### 6. Coverage analysis  ===========================================

# Run coverage analysis if necessary
if (check_coverage == T) {
  # Set number of cores to use
  registerDoParallel(cores = detectCores(logical = TRUE) - l_params_calib$n_cores_reserved_local)
  
  # Generate random set of inputs
  m_param_samp <- with(l_params_calib, {
    # Get number of params to calibrate and number of samples
    n_param <- nrow(prior_map)
    n_samp <- n_samp_coverage
    
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
  )
  
  # Run model for each input parameter sample and get corresponding outputs with parallel processing
  stime <- system.time({
    m_outputs <- foreach(
      i=1:nrow(m_param_samp), 
      .combine=rbind, 
      .inorder=TRUE, 
      .packages=c("data.table","tidyverse")) %dopar% {
        # Get row of parameters and calculate outputs
        v_params_update <- m_param_samp[i,]
        v_calib_outputs <- with(l_params_calib, {
          params_to_calib_outputs(
            l_params_model = l_params_model,
            v_params_update = v_params_update,
            param_map = prior_map,
            l_outcome_params = l_outcome_params,
            l_censor_vars = l_censor_vars
          )
        })
        # Call item to save
        t(v_calib_outputs)
      }
  })
  print(stime)
  closeAllConnections()
  
  # Set column names
  colnames(m_outputs) <- l_params_calib$df_target$target_names
  
  # Save parameter sample and corresponding outputs
  saveRDS(list(m_param_samp = m_param_samp, 
               m_calib_outputs = m_outputs,
               runtime = stime), 
          file = file_coverage)
  
  # Process targets
  df_targets <- l_params_calib$df_target %>%
    mutate(target_index = factor(target_index)) %>% # Create plot labels
    left_join(df_plot_labels, by = "target_groups")
  df_targets$plot_grps <- factor(df_targets$plot_grps, levels = df_plot_labels$plot_grps)
  
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
  plot_coverage <- ggplot(data = df_targets) + 
    geom_errorbar(
      aes(x    = target_index, 
          y    = targets, 
          ymin = targets - se, 
          ymax = targets + se),
      width = 0.4, linewidth = 0.9, color="red") +
    theme(legend.position = "none") +
    geom_violin(data = out_full_bc_cat,
                aes(x    = target_index,
                    y    = value),
                alpha = 0.4) +
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
  plot_coverage
  
  # Save plot and adjust size based on number of rows and columns
  n_plot_grps <- length(unique(df_targets$plot_grps))
  n_plot_rows <- ceiling(n_plot_grps/3)
  n_plot_cols <- ceiling(n_plot_grps/n_plot_rows)
  ggsave(file_fig_coverage, plot = plot_coverage,
         width = n_plot_cols*4, height = 4*n_plot_rows)
}
