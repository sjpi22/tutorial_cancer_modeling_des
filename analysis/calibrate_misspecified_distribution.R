# Program to compare differences in outputs with misspecified distributions

################################################################################
# Setup
################################################################################
# Clear workspace
rm(list = ls())

# Options
options(scipen=999,
        digits=2)

# Load packages
library(readxl)
library(data.table)
library(tidyverse)
library(patchwork)
library(lhs)

# Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)

# Common parameters
v_ages_prevalence <- c(seq(30, 80, 10), 100) # Age for lesion prevalence
v_ages <- list(prevalence = v_ages_prevalence,
               incidence = v_ages_prevalence) # Age for cancer incidence 

# Load default parameters
l_params_model <- load_default_params()

################################################################################
# True distributions: After onset, all gamma distributions
################################################################################

l_params_true <- copy(l_params_model)

l_params_true$time_0_1$params <- list(shape=2, scale=75)

l_params_true$time_1ii_2$distr <- "gamma"
l_params_true$time_1ii_2$params <- list(shape=5/2, scale=2)

l_params_true$time_1i_2$distr <- "gamma"
l_params_true$time_1i_2$params <- list(shape=7/2, scale=2)

l_params_true$time_1i_1ii$distr <- "gamma"
l_params_true$time_1i_1ii$params <- list(shape=4/2, scale=2)


# Map variables to parameters for tuning - make dataframe of all parameters with "src = unknown"
param_map_true <- make_param_map(l_params_true)

#### Initialize population and disease natural history ####
results_true <- run_model(l_params_true)
results_noscreening_true <- results_true[['None']]

# Get prevalence, incidence, and stage outputs
l_outputs_true <- calc_calib_targets(l_params_true, 
                                      results_noscreening_true, 
                                      v_ages)


################################################################################
# Experiments with incidence
################################################################################

# Incidence is similar to density
l_outputs_true$incidence

################################################################################
# Get coverage for exponential distributions
################################################################################

# Map variables to parameters for tuning - make dataframe of all parameters with "src = unknown"
param_map_model <- make_param_map(l_params_model)

# Weibull distribution for time_0_1 should be at least the prevalence in state 1
-log(1 - l_outputs_true$prevalence)
1 - exp(-age_end / beta)^alpha >= prevalence


#### Initialize population and disease natural history ####
results_model <- run_model(l_params_model)
results_noscreening_model <- results_model[['None']]

# Get prevalence, incidence, and stage outputs
l_outputs_model <- calc_calib_targets(l_params_model, 
                                     results_noscreening_model, 
                                     v_ages)

lhs::randomLHS()

l_targets <- list(True = l_outputs_true,
                  Model = l_outputs_model)
plot_targets <- plot_calib_targets(l_targets)

(plot_targets$prevalence + plot_targets$incidence + plot_targets$stage_distr) + 
  plot_layout(
    ncol = 3,
    guides = "collect"
  ) + plot_annotation(title = 'Calibration targets') &
  theme(legend.position = 'bottom')
