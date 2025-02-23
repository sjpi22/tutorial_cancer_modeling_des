###########################  Generate IMABC outputs   ##########################
#
#  Objective: Script to generate calibration target outputs for parameter 
#  posterior distributions calibrated with IMABC
########################### <<<<<>>>>> #########################################

rm(list = ls()) # Clean environment
options(scipen = 999) # View data without scientific notation

#### 1.Libraries and functions  ==================================================

###### 1.1 Load packages
library(tidyverse)
library(dplyr)
library(doBy)
library(patchwork)

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

# Load IMABC parameters from configs file
list2env(configs$params_imabc, envir = .GlobalEnv)


#### 3. Pre-processing actions  ===========================================

# Load model parameters
l_params_calib <- readRDS(file_params_calib)

# Load IMABC outputs and extract parameter samples
calibration_outputs <- readRDS(file_posterior)
m_param_samp <- calibration_outputs$good_parm_draws %>%
  dplyr::select(l_params_calib$prior_map$var_id)

# Load plot labels
df_plot_labels <- read.csv(file_plot_labels)

# Process target data
df_targets <- l_params_calib$df_target %>%
  mutate(target_index = factor(target_index)) %>% # Create plot labels
  left_join(df_plot_labels, by = "target_groups")
df_targets$plot_grps <- factor(df_targets$plot_grps, levels = df_plot_labels$plot_grps)


#### 4. Internal validation  ===========================================

# Sample with replacement
n_imabc_sample <- nrow(calibration_outputs$good_parm_draws) * 0.3
indices_imabc_sample <- sample(1:nrow(calibration_outputs$good_parm_draws), 
                               n_imabc_sample, 
                               replace = TRUE, 
                               prob = calibration_outputs$good_parm_draws$sample_wt)

# IMABC posterior parameters
imabc_params <- calibration_outputs$good_parm_draws[indices_imabc_sample, ] %>%
  dplyr::select(all_of(l_params_calib$prior_map$var_id))

df_params_long <- imabc_params %>%
  pivot_longer(everything())

# Get posterior outputs
imabc_targets_unweighted <- calibration_outputs$good_sim_target %>%
  dplyr::select(-c("iter", "draw", "step"))

m_outputs <- as.matrix(imabc_targets_unweighted[indices_imabc_sample, ])

plt_coverage <- plot_coverage(df_targets = df_targets,
                              m_outputs = m_outputs,
                              file_fig_coverage = NULL)
plt_coverage


# Categorize categorical targets
df_targets <- df_targets %>%
  group_by(target_groups) %>%
  mutate(categorical = (target_groups %in% c(l_params_calib$v_outcomes_categorical) | n()==1))
  
# Plot continuous items
out_summary_cont <- out_summary %>%
  filter(categorical == 0)

plot_targets_cont <- ggplot(data = out_summary_cont, 
                            aes(x    = target_index, 
                                y    = true_val, 
                                ymin = true_val - true_se, 
                                ymax = true_val + true_se))+ 
  geom_errorbar(width=.4, linewidth=0.9, color="red") +
  theme(legend.position="none") +
  geom_ribbon(data = out_summary_cont,
              aes(x    = target_index,
                  y    = true_val,
                  ymin = model_LB_95,
                  ymax = model_UB_95),
              fill = "black",
              alpha = 0.3) +
  geom_ribbon(data = out_summary_cont,
              aes(x    = target_index,
                  y    = true_val,
                  ymin = model_LB_50,
                  ymax = model_UB_50),
              fill = "black",
              alpha = 0.5) +
  facet_wrap(~ plot_grp, scales="free", ncol = 3,
             labeller = labeller(plot_grp = label_wrap_gen(12))) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(), legend.position="none") +
  scale_fill_manual(values = c("grey10", "grey30"))+
  scale_y_continuous(breaks = number_ticks(5))+
  theme_bw(base_size = 23) +
  theme(plot.title = element_text(size = 22, face = "bold"),
        axis.text.x = element_text(size = 18, angle = 90),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0)) +
  labs(x     = "Age", y     = "")

# Categorical items
out_summary_cat <- out_summary %>%
  filter(categorical == 1) %>%
  mutate(target_index = factor(target_index))

# Get full data, filter to categorical, and convert wide to long
out_full_cat <- data.frame(m_outputs[, out_summary$categorical == TRUE]) %>%
  pivot_longer(
    cols = everything(), 
    names_to = "target_groups",
    values_to = "value"
  ) %>%
  mutate(target_index = rep(out_summary_cat$target_index, nrow(m_outputs)),
         plot_grp = rep(out_summary_cat$plot_grp, nrow(m_outputs)))

# Single items as violin plots
plot_targets_cat <- ggplot(data = out_summary_cat) + 
  geom_errorbar(
    aes(x    = target_index, 
        y    = true_val, 
        ymin = true_val - true_se, 
        ymax = true_val + true_se),
    width=.4, linewidth=0.9, color="red") +
  theme(legend.position="none") +
  geom_violin(data = out_full_cat,
              aes(x    = target_index,
                  y    = value),
              alpha = 0.4) +
  facet_wrap(~ plot_grp, scales="free", 
             ncol = length(unique(out_summary_cat$plot_grp)),
             labeller = labeller(plot_grp = label_wrap_gen(12))) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(), legend.position="none") +
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
        strip.text = element_text(hjust = 0)) +
  labs(x     = "", y     = "")

plot_all <- plot_targets_cont / plot_targets_cat
plot_all
ggsave(file_fig_validation, plot_all, width = 10, height = 8)
