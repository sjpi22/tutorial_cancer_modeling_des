###########################  Unit tests  ##########################
#
#  Objective: Tests and sanity checks to validate functions
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
library(survival)
library(survminer)
library(doParallel)
library(foreach)
library(assertthat)

# Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)


#### 2. General parameters ========================================================

###### 2.1 file paths
file_truth <- "ground_truth/true_param_map.rds"
file_params <- "data/calibration_params.rds"
outpath <- "output/calibration/BayCANN"

###### 2.2 other parameters
n_samp <- 100 # Number of samples to simulate
n_cohort_orig <- 100000
n_screen_sample_orig <- 20000
seed_orig <- 2024


#### 3. Validate BayCANN sample generation =====================================

# Load model and ground truth params
l_params_calib <- readRDS(file_params)
df_true_params <- readRDS(file_truth)
v_params_update <- df_true_params$param_val

# Make copy of calibration parameters based on ground truth
l_params_calib_orig <- copy(l_params_calib)
l_params_calib_orig$l_params_all$n_cohort <- n_cohort_orig
l_params_calib_orig$l_params_all$seed <- seed_orig
l_params_calib_orig$l_outcome_params$prevalence$n_screen_sample <- n_screen_sample_orig

# Set seed
set.seed(seed_orig)

# Replicate ground truth values with IMABC function
fn <- function(v_params_update) {
  v_targets <- with(l_params_calib_orig, {
    params_to_calib_outputs(
      l_params_all = l_params_all,
      v_params_update = v_params_update,
      param_map = prior_map,
      l_outcome_params = l_outcome_params,
      l_censor_vars = l_censor_vars
    )
  })
  return(v_targets)
}
v_calib_outputs_rep <- fn(v_params_update)

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
      v_calib_outputs <- with(l_params_calib, {
        params_to_calib_outputs(
          l_params_all = l_params_all,
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

# Process targets
df_targets <- l_params_calib$df_true_targets %>%
  mutate(target_index = factor(target_index)) %>% # Create plot labels
  mutate(plot_grp = case_when(target_groups == "stage_distr" ~ "Stage distribution",
                              target_groups == "incidence" ~ "Incidence per 100k by age",
                              target_groups == "prevalence" ~ "Prevalence by age"))
df_targets$rep_targets <- v_calib_outputs_rep

# Convert outputs from wide to long
out_full_bc_cat <- data.frame(m_outputs) %>%
  pivot_longer(
    cols = everything(), 
    names_to = "target_names",
    values_to = "value"
  ) %>%
  mutate(target_index = rep(df_targets$target_index, nrow(m_outputs)),
         plot_grp = rep(df_targets$plot_grp, nrow(m_outputs)))

# Plot distribution of outputs against targets
plot_coverage <- ggplot(data = df_targets) + 
  geom_point(
    aes(x = target_index, 
        y = rep_targets),
    color = "black") +
  geom_errorbar(
    aes(x    = target_index, 
        y    = targets, 
        ymin = targets - se, 
        ymax = targets + se),
    width = 0.4, linewidth = 0.9, color="red") +
  theme(legend.position="none") +
  geom_violin(data = out_full_bc_cat,
              aes(x    = target_index,
                  y    = value),
              alpha = 0.4) +
  facet_wrap(~ plot_grp, scales="free") +
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
plot_coverage


################################################################################
# Parameters
################################################################################

#### Modifiable parameters ####
n_cohort <- 100000 # Number to simulate in cohort
n_screen_sample <- 10000

# Randomization 
seed <- 1 # Random seed for generating data

# Outcome reporting
v_ages_prevalence <- seq(30, 80, 10) # Age for lesion prevalence
v_ages <- list(prevalence = v_ages_prevalence,
               incidence = c(v_ages_prevalence, 100)) # Age for cancer incidence 
v_time_surv <- seq(0, 10) # Times from event to calculate relative survival

################################################################################
# Test character cancer stages
################################################################################

# Load default data
l_params_all <- load_default_params(v_cancer = c('a', 'b'),
                                    file.surv = NULL)

# Add true survival distribution of exponential from diagnosis
l_params_all$time_Ca_Dc$distr <- "exp"
l_params_all$time_Ca_Dc$params <- list(rate = 0.1)
l_params_all$time_Cb_Dc$distr <- "exp"
l_params_all$time_Cb_Dc$params <- list(rate = 0.3)

# Map variables to parameters for tuning - make dataframe of all parameters with "src = unknown"
param_map <- make_param_map(l_params_all)

# Set "true" parameters
v_param_update <- c(0.42, 0.25, 0.5, 3, 300)

# Update params
l_params_all <- update_param_from_map(l_params_all, v_param_update, param_map)
l_params_all <- update_param_list(l_params_all,
                                  list(seed = seed,
                                       n_cohort = n_cohort,
                                       v_strats = l_params_all$v_strats[1]))

# Update parameter map
param_map$param_val <- v_param_update

#### Initialize population and disease natural history ####
results <- run_model(l_params_all)
results_noscreening <- results[['None']]

################################################################################
# Test more than 2 cancer stages
################################################################################

# Update cancer stages
l_params_all$v_cancer = c('a', 'b', 'c', 'd')

# Add progression within preclinical stage variables
l_params_all$time_Pb_Pc$distr <- "exp"
l_params_all$time_Pb_Pc$params <- list(rate = 0.2)
l_params_all$time_Pc_Pd$distr <- "exp"
l_params_all$time_Pc_Pd$params <- list(rate = 0.3)

# Add detection variables
l_params_all$time_Pc_C$distr <- "exp"
l_params_all$time_Pc_C$params <- list(rate = 0.4)
l_params_all$time_Pd_C$distr <- "exp"
l_params_all$time_Pd_C$params <- list(rate = 0.6)

# Add true survival distribution of exponential from diagnosis
l_params_all$time_Cc_Dc$distr <- "exp"
l_params_all$time_Cc_Dc$params <- list(rate = 0.4)
l_params_all$time_Cd_Dc$distr <- "exp"
l_params_all$time_Cd_Dc$params <- list(rate = 0.5)

#### Initialize population and disease natural history ####
results <- run_model(l_params_all)
results_noscreening <- results[['None']]

# Check sums
are_equal(results_noscreening[stage_dx == 'a'], results_noscreening[time_Pa_Pb > time_Pa_C])
are_equal(results_noscreening[stage_dx == 'b'], results_noscreening[time_Pb_Pc > time_Pb_C])
are_equal(results_noscreening[stage_dx == 'c'], results_noscreening[time_Pc_Pd > time_Pc_C])

are_equal(abs(sum(results_noscreening[stage_dx == 'a', time_H_C] - rowSums(results_noscreening[stage_dx == 'a', c("time_H_P", "time_Pa_C")]))) < 1e-6, TRUE)
are_equal(abs(sum(results_noscreening[stage_dx == 'a', time_H_Dc] - rowSums(results_noscreening[stage_dx == 'a', c("time_H_P", "time_Pa_C", "time_C_Dc")]))) < 1e-6, TRUE)
are_equal(abs(sum(results_noscreening[stage_dx == 'b', time_H_C] - rowSums(results_noscreening[stage_dx == 'b', c("time_H_P", "time_Pa_Pb", "time_Pb_C")]))) < 1e-6, TRUE)
are_equal(abs(sum(results_noscreening[stage_dx == 'b', time_H_Dc] - rowSums(results_noscreening[stage_dx == 'b', c("time_H_P", "time_Pa_Pb", "time_Pb_C", "time_C_Dc")]))) < 1e-6, TRUE)
are_equal(abs(sum(results_noscreening[stage_dx == 'c', time_H_Dc] - rowSums(results_noscreening[stage_dx == 'c', c("time_H_P", "time_Pa_Pb", "time_Pb_Pc", "time_Pc_C", "time_C_Dc")]))) < 1e-6, TRUE)
are_equal(abs(sum(results_noscreening[stage_dx == 'd', time_H_Dc] - rowSums(results_noscreening[stage_dx == 'd', c("time_H_P", "time_Pa_Pb", "time_Pb_Pc", "time_Pc_Pd", "time_Pd_C", "time_C_Dc")]))) < 1e-6, TRUE)

are_equal(abs(sum(results_noscreening[, time_H_D] - pmin(results_noscreening[, time_H_Dc], results_noscreening[, time_H_Do], na.rm = TRUE))) < 1e-6, TRUE)


################################################################################
# Test having precancerous lesion stage
################################################################################
l_params_all <- load_default_params(v_states = c('H', 'L', 'P', 'C', 'D'),
                                    file.surv = NULL)

# Add true survival distribution of exponential from diagnosis
l_params_all$time_C1_Dc$distr <- "exp"
l_params_all$time_C1_Dc$params <- list(rate = 0.1)
l_params_all$time_C2_Dc$distr <- "exp"
l_params_all$time_C2_Dc$params <- list(rate = 0.3)

# Map variables to parameters for tuning - make dataframe of all parameters with "src = unknown"
param_map <- make_param_map(l_params_all)

# Set "true" parameters
v_param_update <- c(0.42, 0.25, 0.5, 0.1, 3, 300)

# Update params
l_params_all <- update_param_from_map(l_params_all, v_param_update, param_map)
l_params_all <- update_param_list(l_params_all,
                                  list(seed = seed,
                                       n_cohort = n_cohort,
                                       v_strats = l_params_all$v_strats[1]))

# Update parameter map
param_map$param_val <- v_param_update

#### Initialize population and disease natural history ####
results <- run_model(l_params_all)
results_noscreening <- results[['None']]

# Check sums
are_equal(results_noscreening[stage_dx == 2], results_noscreening[time_P1_P2 < time_P1_C])
are_equal(abs(sum(results_noscreening[stage_dx == 2, time_H_Dc] - rowSums(results_noscreening[stage_dx == 2, c("time_H_L", "time_L_P", "time_P1_P2", "time_P2_C", "time_C_Dc")]))) < 1e-6, TRUE)
