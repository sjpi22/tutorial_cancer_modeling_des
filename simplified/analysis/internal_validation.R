###########################  Calibrated model   #########################################
#
#  Objective: Script to calculate decision model outputs with BayCANN calibrated parameters
########################### <<<<<>>>>> ##############################################


#### 1.Libraries and functions  ==================================================
library(tidyverse)
library(doBy)
library(dplyr)
library(GGally)
library(patchwork)

###### 1.1 Load functions =================================================
#* Clean environment
# rm(list = ls())

# Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)


#### 2. General parameters ========================================================

# Control variables for running on cluster and/or parallelized
run_parallel <- TRUE
verbose <- FALSE

###### 2.1 file paths 
prior_file <- "data/priors.RData"

targets_files <- list(prevalence = "data/prevalence_asymptomatic_cancer.csv",
                      incidence = "data/incidence_symptomatic_cancer.csv",
                      stage_distr = "data/stage_distr.csv")

path_baycann_post <- "output/calibrated_posteriors_BayCANN.csv"
path_baycann_target <- "output/calibration_targets_BayCANN.csv"
path_imabc <- "output/imabc/IMABC_calibration.RData"

###### 2.1 model parameters 

# Load default data
l_params_all <- load_default_params()

#### 3. Pre-processing actions  ===========================================

# Set random seed
set.seed(l_params_all$seed)

# Load parameter mapping
load(prior_file) 

# Load targets
l_true_targets <- recursive_read_csv(targets_files)

# Reshape true targets and get mapping
l_true_reshaped <- reshape_calib_targets(l_true_targets, output_se = TRUE, output_map = TRUE)

# Get vector of ages for prevalence and incidence
v_ages_prevalence <- get_age_range(l_true_targets$prevalence)
v_ages_incidence <- get_age_range(l_true_targets$incidence)
v_ages <- list(prevalence = v_ages_prevalence,
               incidence = v_ages_incidence)

# Load true parameters
load("_solutions/true_param_map.RData") # Output: param_map

# Load priors
load("data/priors.RData")

################################################################################
###  Generate IMABC outputs                                          ###
################################################################################

# Note: check convergence of IMABC
# Note: Try mutliple centers

# Load IMABC results
load(path_imabc)

# Sample with replacement
n_imabc_sample <- nrow(calibration_results$good_parm_draws) * 10
indices_imabc_sample <- sample(1:nrow(calibration_results$good_parm_draws), 
                               n_imabc_sample, 
                               replace = TRUE, 
                               prob = calibration_results$good_parm_draws$sample_wt)

# IMABC posterior
imabc_params_unweighted <- calibration_results$good_parm_draws %>%
  dplyr::select(all_of(prior_map$var_id))

imabc_params <- imabc_params_unweighted[indices_imabc_sample, ]

imabc_params_long <- imabc_params %>%
  pivot_longer(everything())

# Plot IMABC parameters against true parameters and priors
plot_params_imabc <- ggplot(imabc_params_long, aes(value)) + 
  geom_histogram() + 
  geom_vline(data = param_map %>%
               rename(name = var_id), aes(xintercept=param_val), color = 'red') +
  geom_vline(data = prior_map %>%
               rename(name = var_id), aes(xintercept=prior_min), color = 'blue') +
  geom_vline(data = prior_map %>%
               rename(name = var_id), aes(xintercept=prior_max), color = 'blue') +
  facet_wrap(~name, scales = "free") +
  ggtitle('IMABC')
  

# Correlation plot
gg_calib_post_pair_corr <- GGally::ggpairs(imabc_params,
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
        axis.text.x  = element_text(size=6),
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_rect(fill = "white",
                                        color = "white"),
        strip.text = element_text(hjust = 0))
gg_calib_post_pair_corr

# Get posterior mean of outputs
imabc_targets_unweighted <- calibration_results$good_sim_target %>%
  dplyr::select(-c("iter", "draw", "step"))

imabc_targets <- imabc_targets_unweighted[indices_imabc_sample, ]

# Summarize targets
# Note: weight for 50th and 95th with resampling
collapse_mean  <- summaryBy( . ~ index , FUN=c(weighted.mean), w = calibration_results$good_parm_draws$sample_wt, data=imabc_targets_unweighted, keep.names=TRUE)
collapse_UB_95 <- summaryBy( . ~ index , FUN=quantile, probs = 0.975, data=imabc_targets, keep.names=TRUE)
collapse_LB_95 <- summaryBy( . ~ index , FUN=quantile, probs = 0.025, data=imabc_targets, keep.names=TRUE)
collapse_UB_50 <- summaryBy( . ~ index , FUN=quantile, probs = 0.75, data=imabc_targets, keep.names=TRUE)
collapse_LB_50 <- summaryBy( . ~ index , FUN=quantile, probs = 0.25, data=imabc_targets, keep.names=TRUE)

out_summary <-  data.frame(l_true_reshaped$target_map,
                           true_val = l_true_reshaped$v_targets,
                           true_se = l_true_reshaped$v_se,
                           model_mean = unname(t(collapse_mean)),
                           model_UB_95 = unname(t(collapse_UB_95)),
                           model_LB_95 = unname(t(collapse_LB_95)),
                           model_UB_50 = unname(t(collapse_UB_50)),
                           model_LB_50 = unname(t(collapse_LB_50))) %>%
  mutate(categorical = (target_groups %in% c('Stage at diagnosis', 'Prevalence of lesion type b')))

# Plot 
plot_targets_imabc <- ggplot(data = out_summary, 
               aes(x    = target_index, 
                   y    = true_val, 
                   ymin = true_val - true_se, 
                   ymax = true_val + true_se))+ 
  geom_errorbar(width=.4, size=0.9, color="red") +
  theme(legend.position="none") +
  geom_errorbar(data = out_summary[out_summary$categorical==1,],
                aes(x    = target_index-0.07,
                    y    = model_mean,
                    ymin = model_LB_95,
                    ymax = model_UB_95),width=.4, size=0.7, color="black", alpha = 0.7, position = "dodge2") +
  geom_ribbon(data = out_summary[out_summary$categorical==0,],
              aes(x    = target_index,
                  y    = true_val,
                  ymin = model_LB_95,
                  ymax = model_UB_95),
              fill = "black",
              alpha = 0.3) +
  geom_ribbon(data = out_summary[out_summary$categorical==0,],
              aes(x    = target_index,
                  y    = true_val,
                  ymin = model_LB_50,
                  ymax = model_UB_50),
              fill = "black",
              alpha = 0.5) +
  facet_wrap(~ target_groups,scales="free", ncol = 3) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(), legend.position="none") +
  scale_fill_manual(values = c("grey10", "grey30"))+
  scale_y_continuous(breaks = number_ticks(5))+
  theme_bw(base_size = 23) +
  theme(plot.title = element_text(size = 22, face = "bold"),
        axis.text.x = element_text(size = 12, angle = 90),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0)) +
  labs(title = "IMABC calibration", 
       x     = "", y     = "")


################################################################################
###  Generate BayCANN plots                                          ###
################################################################################

# Load BayCANN calibration outputs
baycann_posteriors <- read_csv(path_baycann_post) %>%
  dplyr::select(-lp)
baycann_targets <- read_csv(path_baycann_target)

# Rename posterior columns and create long df
names(baycann_posteriors) <- prior_map$var_id

baycann_posteriors_long <- baycann_posteriors %>%
  pivot_longer(everything())

# Plot BayCANN parameters against true parameters and priors
plot_params_baycann <- ggplot(baycann_posteriors_long, aes(value)) + 
  geom_histogram() + 
  geom_vline(data = param_map %>%
               rename(name = var_id), aes(xintercept=param_val), color = 'red') +
  geom_vline(data = prior_map %>%
               rename(name = var_id), aes(xintercept=prior_min), color = 'blue') +
  geom_vline(data = prior_map %>%
               rename(name = var_id), aes(xintercept=prior_max), color = 'blue') +
  facet_wrap(~name, scales = "free") +
  ggtitle('BayCANN')


# Summarize targets
collapse_mean_bc  <- summaryBy( . ~ index , FUN=c(mean), data=baycann_targets,keep.names=TRUE)
collapse_UB_95_bc <- summaryBy( . ~ index , FUN=quantile, probs = 0.975, data=baycann_targets, keep.names=TRUE)
collapse_LB_95_bc <- summaryBy( . ~ index , FUN=quantile, probs = 0.025, data=baycann_targets, keep.names=TRUE)
collapse_UB_50_bc <- summaryBy( . ~ index , FUN=quantile, probs = 0.75, data=baycann_targets, keep.names=TRUE)
collapse_LB_50_bc <- summaryBy( . ~ index , FUN=quantile, probs = 0.25, data=baycann_targets, keep.names=TRUE)

out_summary_bc <-  data.frame(l_true_reshaped$target_map,
                           true_val = l_true_reshaped$v_targets,
                           true_se = l_true_reshaped$v_se,
                           model_mean = unname(t(collapse_mean_bc)),
                           model_UB_95 = unname(t(collapse_UB_95_bc)),
                           model_LB_95 = unname(t(collapse_LB_95_bc)),
                           model_UB_50 = unname(t(collapse_UB_50_bc)),
                           model_LB_50 = unname(t(collapse_LB_50_bc))) %>%
  mutate(categorical = (target_groups %in% c('Stage at diagnosis', 'Prevalence of lesion type b')))

# Plot 
plot_targets_baycann <- ggplot(data = out_summary_bc, 
                aes(x    = target_index, 
                    y    = true_val, 
                    ymin = true_val - true_se, 
                    ymax = true_val + true_se))+ 
  geom_errorbar(width=.4, size=0.9, color="red") +
  theme(legend.position="none") +
  geom_errorbar(data = out_summary_bc[out_summary_bc$categorical==1,],
                aes(x    = target_index-0.07,
                    y    = model_mean,
                    ymin = model_LB_95,
                    ymax = model_UB_95),width=.4, size=0.7, color="black", alpha = 0.7, position = "dodge2") +
  geom_ribbon(data = out_summary_bc[out_summary_bc$categorical==0,],
              aes(x    = target_index,
                  y    = true_val,
                  ymin = model_LB_95,
                  ymax = model_UB_95),
              fill = "black",
              alpha = 0.3) +
  geom_ribbon(data = out_summary_bc[out_summary_bc$categorical==0,],
              aes(x    = target_index,
                  y    = true_val,
                  ymin = model_LB_50,
                  ymax = model_UB_50),
              fill = "black",
              alpha = 0.5) +
  facet_wrap(~ target_groups,scales="free", ncol = 3) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(), legend.position="none") +
  scale_fill_manual(values = c("grey10", "grey30"))+
  scale_y_continuous(breaks = number_ticks(5))+
  theme_bw(base_size = 23) +
  theme(plot.title = element_text(size = 22, face = "bold"),
        axis.text.x = element_text(size = 12, angle = 90),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0)) +
  labs(title = "BayCANN calibration", 
       x     = "", y     = "")

# Plot both together

(plot_params_imabc / plot_params_baycann) +
  plot_layout(
    guides = "collect"
  ) + plot_annotation(title = 'Prior vs. posterior distributions') 

(plot_targets_imabc / plot_targets_baycann) +
  plot_layout(
    guides = "collect"
  ) + plot_annotation(title = 'Calibration targets') 


