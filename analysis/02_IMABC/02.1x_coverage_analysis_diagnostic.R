###########################  Coverage Analysis Diagnostics  ##########################
#
#  Objective: Program to check which parameter sets were in bounds
# 
########################### <<<<<>>>>> #########################################

rm(list = ls()) # Clean environment
options(scipen = 999) # View data without scientific notation

#### 1.Libraries and functions  ==================================================

###### 1.1 Load packages
library(tools)
library(data.table)
library(tidyverse)
library(cobs)
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
# Load configs
file_configs <- file.path("configs", "configs_bladder.yaml")
configs <- load_configs(file_configs)

# Extract relevant parameters from configs
params_model <- configs$params_model
params_calib <- configs$params_calib
file_params_calib <- configs$paths$file_params_calib
file_plot_labels <- configs$paths$file_plot_labels
alpha_ci <- 1 - configs$params_model$conf_level

# Get list of relevant output file paths and load to global environment
l_filepaths <- update_config_paths("files_coverage", configs$paths)
list2env(l_filepaths, envir = .GlobalEnv)

# Load coverage analysis parameters from configs file
list2env(configs$params_coverage, envir = .GlobalEnv)


#### 3. Pre-processing  ===========================================

# Load model and calibration parameters
l_params_calib <- readRDS(file_params_calib)

# Load plot labels
df_plot_labels <- read.csv(file_plot_labels)

# Process targets
df_targets <- l_params_calib$df_target %>%
  # Set CIs if not already set
  mutate(ci_lb = ifelse(is.na(ci_lb), targets - se*qnorm(1 - alpha_ci/2), ci_lb),
         ci_ub = ifelse(is.na(ci_ub), targets + se*qnorm(1 - alpha_ci/2), ci_ub)) %>% 
  # Merge plot labels
  left_join(df_plot_labels, by = "target_groups")
df_targets$plot_grps <- factor(df_targets$plot_grps, levels = df_plot_labels$plot_grps)

# Read coverage data
df_coverage <- readRDS(file_coverage)

# Were there any parameter sets within the bounds?
df_above_lb <- sweep(df_coverage$m_calib_outputs, 2, df_targets$current_lower_bounds, FUN = ">=") 
df_below_ub <- sweep(df_coverage$m_calib_outputs, 2, df_targets$current_upper_bounds, FUN = "<=") 
df_within_bounds <- df_above_lb & df_below_ub
test <- apply(df_within_bounds, 1, min) # If all elements in row were within bounds, vector entry for that index will be 1
print(sum(test)) # Number within bounds

# Find parameter sets with the most in bounds
n_in_bounds <- rowSums(df_within_bounds)

# Examine top 3
idx_top <- order(n_in_bounds, decreasing = TRUE)[1:3]

plt_coverage <- plot_coverage(df_targets = df_targets, 
                              m_outputs = df_coverage$m_calib_outputs[idx_top, ], 
                              target_range = "ci",
                              plt_size_text = plt_size_text,
                              labeller_multiplier = 6)
plt_coverage

df_coverage$m_param_samp[idx_top[1],]

