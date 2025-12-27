###########################  Unit test: IMABC and BayCANN functions  ##########################
#
#  Objective: Visual checks for IMABC and BayCANN functions
########################### <<<<<>>>>> #########################################

#### 1.Libraries and functions  ==================================================
# Clear workspace
rm(list = ls())

# Options
options(scipen=999)

# Load packages
library(readxl)
library(data.table)
library(tidyverse)
library(VGAM)
library(doParallel)
library(foreach)
library(assertthat)

# Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)

# IMABC function
fn <- function(v_params_update) {
  v_targets <- with(l_params_calib_orig, {
    params_to_outputs(
      l_params_model = l_params_true,
      param_map = prior_map,
      l_params_outcome = l_params_outcome,
      l_censor_vars = l_censor_vars
    )
  })
  return(v_targets)
}


#### 2. General parameters ========================================================

###### 2.1 file paths
file_truth <- "_ground_truth/params_model.rds"

###### 2.2 other parameters
n_samp <- 100 # Number of samples to simulate
n_cohort_orig <- 100000
n_screen_sample_orig <- 20000
seed_orig <- 2024


#### 3. Pre-processing =====================================

# Load configs
file_configs <- file.path("configs", "configs_simulated.yaml")
configs <- load_configs(file_configs)

# Extract relevant parameters from configs
file_plot_labels <- configs$paths$file_plot_labels
file_params <- configs$paths$file_params_calib

# Load model and ground truth params
l_params_calib <- readRDS(file_params)
l_params_true <- readRDS(file_truth)

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

# Make copy of calibration parameters based on ground truth
l_params_calib_orig <- copy(l_params_calib)
l_params_calib_orig$l_params_all$n_cohort <- n_cohort_orig
l_params_calib_orig$l_params_all$seed <- seed_orig
l_params_calib_orig$l_outcome_params$prevalence$n_screen_sample <- n_screen_sample_orig

# Set seed
set.seed(seed_orig, kind = "L'Ecuyer-CMRG")

# If running locally, use all available cores except for reserved ones
registerDoParallel(cores = detectCores(logical = TRUE) - l_params_calib$n_cores_reserved_local)


#### 4. Validate IMABC and BayCANN functions =====================================

##### 4.1 IMABC
# Run ground truth values with IMABC function and join to targets
v_calib_outputs_rep <- fn(v_params_update)
df_targets$rep_targets <- v_calib_outputs_rep

# Plot distribution of IMABC outputs against targets
plt_coverage_imabc <- ggplot(data = df_targets) + 
  geom_point( # IMABC outputs
    aes(x = target_index, 
        y = rep_targets),
    color = "black") +
  geom_errorbar( # True targets
    aes(x    = target_index, 
        y    = targets, 
        ymin = targets - se, 
        ymax = targets + se),
    width = 0.4, linewidth = 0.9, color="red") +
  facet_wrap(~ plot_grps, scales="free") +
  scale_fill_manual(values = c("grey10", "grey30"))+
  scale_y_continuous(breaks = number_ticks(5))+
  theme_bw(base_size = 23) +
  theme(plot.title = element_text(size = 22, face = "bold"),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0), 
        legend.position="none") +
  labs(x     = "", y     = "", title = "IMABC")
plt_coverage_imabc


##### 4.2 BayCANN
# Set null seed
l_params_calib$l_params_all$seed <- NULL

# Simulate outputs with true parameters
stime <- system.time({
  m_outputs <- foreach(
    i=1:n_samp, 
    .combine=rbind, 
    .inorder=TRUE, 
    .packages=c("data.table","tidyverse")) %dopar% {
      # Get row of parameters and calculate outputs
      v_calib_outputs <- with(l_params_calib_orig, {
        params_to_outputs(
          l_params_model = l_params_true,
          param_map = prior_map,
          l_params_outcome = l_params_outcome,
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
colnames(m_outputs) <- df_targets$target_names

# Plot distribution of BayCANN outputs against targets
plt_coverage_baycann <- plot_coverage(df_targets = df_targets, 
                              m_outputs = m_outputs,
                              target_range = "ci")
plt_coverage_baycann
