###########################  Internal Validation  #########################################
#
#  Objective: Validate BayCANN posteriors by plotting fit of calibration outputs
########################### <<<<<>>>>> ##############################################

rm(list = ls()) # Clean environment
options(scipen = 999) # View data without scientific notation

#### 1.Libraries and functions  ==================================================

###### 1.1 Load packages
library(tidyverse)
library(patchwork)

###### 1.2 Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)


#### 2. General parameters ========================================================

###### 2.1 Configurations
# Load configs
file_configs <- file.path("configs", "configs_colorectal.yaml")
configs <- load_configs(file_configs)

# Extract relevant parameters from configs
file_params_calib <- configs$paths$file_params_calib
file_plot_labels <- configs$paths$file_plot_labels
n_chains <- configs$params_baycann$params_stan$n_chains
plot_by_chain <- configs$params_baycann$params_validation$plot_by_chain
if (is.null(plot_by_chain)) plot_by_chain <- FALSE

# Get list of BayCANN output file paths and load to global environment
l_filepaths <- update_config_paths("files_baycann", configs$paths)
list2env(l_filepaths, envir = .GlobalEnv)

# Load coverage analysis parameters from configs file
list2env(configs$params_coverage, envir = .GlobalEnv)


#### 3. Pre-processing actions  ===========================================

# Load model parameters
l_params_calib <- readRDS(file_params_calib)

# Load calibration outputs
l_outputs <- readRDS(file_outputs)$l_outputs

# Load plot labels
df_plot_labels <- read.csv(file_plot_labels)

# Process target data
df_targets <- l_params_calib$df_target %>%
  left_join(df_plot_labels, by = "target_groups") %>%
  group_by(target_groups) %>%
  mutate(categorical = (target_groups %in% c(l_params_calib$v_outcomes_categorical) | n()==1))
df_targets$plot_grps <- factor(df_targets$plot_grps, levels = df_plot_labels$plot_grps)

# Calculate quantiles and column labels from inner quantile vector
v_quantiles_lb <- (1 - v_quantiles/100)/2
names(v_quantiles_lb) <- paste0("model_LB_", v_quantiles)
v_quantiles_ub <- (1 + v_quantiles/100)/2
names(v_quantiles_ub) <- paste0("model_UB_", v_quantiles)
v_quantiles_calc <- sort(c(v_quantiles_lb, v_quantiles_ub))


#### 4. Plots ===========================================

# Extract calibration outputs and convert to data frame
if(is.list(l_outputs)) {
  m_outputs <- do.call(rbind, lapply(l_outputs, function(u) {
    reshape_outputs(u[["outputs_base"]][names(l_params_calib$l_params_outcome)])
  }))
} else {
  m_outputs <- l_outputs
}

# Calculate mean of outputs
df_targets$model_mean <- colMeans(m_outputs)

# Get quantiles of calibration outputs
m_output_quantiles <- t(apply(m_outputs, 2, function(u) {
  quantile(u, probs = v_quantiles_calc)
}))
colnames(m_output_quantiles) <- names(v_quantiles_calc)

# Append quantiles to df_targets
df_targets_outputs <- cbind(df_targets, m_output_quantiles)

# Get and save target plot
plot_targets(df_targets = df_targets_outputs,
             target_range = "ci",
             outfile = file_fig_targets,
             plt_size_text = plt_size_text)

# Get and save coverage plot
plt_coverage <- plot_coverage(df_targets = df_targets_outputs,
                              target_range = "ci",
                              outfile = file_fig_validation,
                              plt_size_text = plt_size_text)
plt_coverage

# Plot by chain if applicable
if (plot_by_chain) {
  n_rows_per_chain <- nrow(m_outputs) / n_chains
  for (i in 1:n_chains) {
    # Get quantiles of calibration outputs
    m_output_quantiles <- t(apply(m_outputs[((i-1)*n_rows_per_chain + 1):(i*n_rows_per_chain),], 2, function(u) {
      quantile(u, probs = v_quantiles_calc)
    }))
    colnames(m_output_quantiles) <- names(v_quantiles_calc)
    
    # Append quantiles to df_targets
    df_targets_outputs <- cbind(df_targets, m_output_quantiles)
    
    # Modify file path
    file_fig_validation_chain <- paste0(tools::file_path_sans_ext(file_fig_validation),
                                        "_chain", i, ".",
                                        tools::file_ext(file_fig_validation))
    
    # Make coverage plot
    plt_coverage_chain <- plot_coverage(df_targets = df_targets_outputs,
                                        target_range = "ci",
                                        outfile = file_fig_validation_chain)
  }
}
