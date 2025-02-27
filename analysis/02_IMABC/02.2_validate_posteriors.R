###########################  Internal Validation  #########################################
#
#  Objective: Validate IMABC posteriors by plotting fit of calibration outputs
########################### <<<<<>>>>> ##############################################

rm(list = ls()) # Clean environment
options(scipen = 999) # View data without scientific notation

#### 1.Libraries and functions  ==================================================

###### 1.1 Load packages
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggdist)
library(GGally)

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
file_params_calib <- configs$paths$file_params_calib
file_plot_labels <- configs$paths$file_plot_labels

# Get list of relevant output file paths and load to global environment
l_filepaths <- update_config_paths("files_imabc", configs$paths)
list2env(l_filepaths, envir = .GlobalEnv)

# Load coverage analysis parameters from configs file
list2env(configs$params_coverage, envir = .GlobalEnv)


#### 3. Pre-processing actions  ===========================================

# Load model parameters
l_params_calib <- readRDS(file_params_calib)

# Load IMABC calibration outputs
calibration_outputs <- readRDS(file_posterior)

# Extract posteriors
m_params <- calibration_outputs$good_parm_draws %>%
  dplyr::select(l_params_calib$prior_map$var_id)

# Extract unweighted calibration outputs
imabc_targets_unweighted <- calibration_outputs$good_sim_target %>%
  dplyr::select(-c("iter", "draw", "step"))

# Load plot labels
df_plot_labels <- read.csv(file_plot_labels)

# Process target data
df_targets <- l_params_calib$df_target %>%
  left_join(df_plot_labels, by = "target_groups") %>%
  group_by(target_groups) %>%
  mutate(categorical = (target_groups %in% c(l_params_calib$v_outcomes_categorical) | n()==1))
df_targets$plot_grps <- factor(df_targets$plot_grps, levels = df_plot_labels$plot_grps)


#### 4. Internal validation  ===========================================

# Calculate weighted mean of outputs
df_targets$model_mean <- apply(imabc_targets_unweighted, 2, weighted.mean, w = calibration_outputs$good_parm_draws$sample_wt)

# Calculate quantiles and column labels from inner quantile vector
v_quantiles_lb <- (1 - v_quantiles/100)/2
names(v_quantiles_lb) <- paste0("model_LB_", v_quantiles)
v_quantiles_ub <- (1 + v_quantiles/100)/2
names(v_quantiles_ub) <- paste0("model_UB_", v_quantiles)
v_quantiles_calc <- sort(c(v_quantiles_lb, v_quantiles_ub))

# Get weighted quantiles of calibration outputs
m_output_quantiles <- t(apply(imabc_targets_unweighted, 2, function(u) {
  weighted_quantile(
    x = u,
    probs = v_quantiles_calc,
    weights = calibration_outputs$good_parm_draws$sample_wt
  )
}))
colnames(m_output_quantiles) <- names(v_quantiles_calc)

# Append quantiles to df_targets
df_targets <- cbind(df_targets, m_output_quantiles)

# Get and save coverage plot
plt_coverage <- plot_coverage(df_targets = df_targets,
                              file_fig_coverage = file_fig_validation)
plt_coverage

# Plot correlation graph
df_post <- m_params
gg_calib_post_pair_corr <- plot_posterior_corr(df_post,
                                               file_fig_corr = file_fig_corr)
gg_calib_post_pair_corr
