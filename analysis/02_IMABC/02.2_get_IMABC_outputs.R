###########################  IMABC outputs   ###################################
#
#  Objective: Script to regenerate calibration target outputs for IMABC 
# calibrated parameters with cohort size consistent with targets
########################### <<<<<>>>>> #########################################


#### 1.Libraries and functions  ==================================================
#* Clean environment
rm(list = ls())

library(tidyverse)
library(dplyr)
library(doBy)
library(patchwork)

###### 1.1 Load functions =================================================

# Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)


#### 2. General parameters ========================================================

###### 2.1 file paths
file_truth <- "ground_truth/true_param_map.rds" 
file_params <- "data/calibration_params.rds"
outpath <- "output/calibration/IMABC"
file_imabc_params <- file.path(outpath, "params_IMABC.rds")
file_posterior <- file.path(outpath, "calibrated_posteriors_IMABC.rds")
file_fig_validation <- file.path(outpath, "plots/fig_internal_validation.png")
file_fig_prior <- file.path(outpath, "plots/fig_param_calibration.png")


#### 3. Pre-processing actions  ===========================================

# Load model parameters
l_params_calib <- readRDS(file_params)

# Load IMABC outputs and extract parameter samples
calibration_outputs <- readRDS(file_posterior)
m_param_samp <- calibration_outputs$good_parm_draws %>%
  dplyr::select(l_params_calib$prior_map$var_id)

# Load true parameters
df_true_params <- readRDS(file_truth)


#### 4. Internal validation  ===========================================

# Sample with replacement
n_imabc_sample <- nrow(calibration_outputs$good_parm_draws) * nrow(calibration_outputs$good_parm_draws)
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

# Summarize targets
# Note: weight for 50th and 95th with resampling
collapse_mean  <- t(summaryBy( . ~ index , FUN=c(weighted.mean), w = calibration_outputs$good_parm_draws$sample_wt, data=imabc_targets_unweighted, keep.names=TRUE))
collapse_UB_95 <- apply(m_outputs, 2, FUN=quantile, probs = 0.975, simplify = TRUE)
collapse_LB_95 <- apply(m_outputs, 2, FUN=quantile, probs = 0.025, simplify = TRUE)
collapse_UB_50 <- apply(m_outputs, 2, FUN=quantile, probs = 0.75, simplify = TRUE)
collapse_LB_50 <- apply(m_outputs, 2, FUN=quantile, probs = 0.25, simplify = TRUE)

out_summary <-  data.frame(
  l_params_calib$df_true_targets,
  model_mean = collapse_mean,
  model_UB_95 = collapse_UB_95,
  model_LB_95 = collapse_LB_95,
  model_UB_50 = collapse_UB_50,
  model_LB_50 = collapse_LB_50) %>%
  group_by(target_groups) %>%
  mutate(categorical = (target_groups %in% c('stage_distr') | n()==1),
         plot_label = as.character(target_index)) %>%
  rename(true_val = targets, true_se = se) %>%
  mutate(plot_grp = case_when(target_groups == 'prevalence' ~ 'Preclinical cancer prevalence',
                              target_groups == 'incidence' ~ 'CRC incidence per 100k',
                              target_groups == 'stage_distr' ~ 'CRC stage distribution'))

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
ggsave(file_fig_validation, plot_all, width = 10, height = 8)

# Plot IMABC parameters against true parameters and priors
plot_params <- ggplot(df_params_long, aes(value)) + 
  geom_histogram() + 
  geom_vline(data = df_true_params %>%
               rename(name = var_id), aes(xintercept=param_val), color = 'red') +
  geom_vline(data = l_params_calib$prior_map %>%
               rename(name = var_id), aes(xintercept=min), color = 'blue') +
  geom_vline(data = l_params_calib$prior_map %>%
               rename(name = var_id), aes(xintercept=max), color = 'blue') +
  facet_wrap(~name, scales = "free")
ggsave(file_fig_prior, plot_params, width = 10, height = 8)
