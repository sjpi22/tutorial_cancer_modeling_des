# Disease simulator for heterogeneous disease with two precursor lesions: a and b
# Default version with data generating function consistent with assumed model

################################################################################
# Setup
################################################################################
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

# Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)

################################################################################
# Parameters
################################################################################

#### Modifiable parameters ####
debug <- FALSE
data_outpath <- 'data/'
true_param_outpath <- '_solutions/'

# Check if directory exists, make if not
dir.create(file.path(data_outpath), showWarnings = FALSE)

# Model structure
if(debug) {
  n_cohort <- 1000 # Number to simulate in cohort
  n_screen_sample <- 100
} else {
  n_cohort <- 100000 # Number to simulate in cohort
  n_screen_sample <- 500
}

# Randomization 
seed <- 1 # Random seed for generating data
jitter_multiplier_low <- 0.5 # Jitter for default parameters - lower bound in uniform distribution
jitter_multiplier_high <- 2 # Jitter for default parameters - upper bound in uniform distribution

# Outcome reporting
v_age_prevalence <- c(seq(20, 85, 5), 100) # Age for lesion prevalence
v_age_cancer_incidence <- v_age_prevalence # Age for cancer incidence 
v_time_surv <- seq(0, 10) # Times from event to calculate relative survival


#### Load and update parameters ####

# Load default data
l_params_all <- load_default_params(file.surv = NULL)

# Jitter disease distribution parameters uniformly from 0.5*default to 2*default
assertthat::are_equal((seed == l_params_all$seed), FALSE)
set.seed(seed)

# Map variables to parameters for tuning - make dataframe of all parameters with "src = unknown"
param_map <- make_param_map(l_params_all)

# Apply jitter to unknown parameters
v_param_update <- sapply(param_map$param_val, function(x) runif(1, jitter_multiplier_low*x, jitter_multiplier_high*x))

# Update params
l_params_all <- update_param_from_map(l_params_all, v_param_update, param_map)
l_params_all <- update_param_list(l_params_all,
                                  list(seed = seed,
                                       n_cohort = n_cohort,
                                       v_strats = l_params_all$v_strats[1]))

# Update parameter map
param_map$param_val <- v_param_update

################################################################################
# Generate population data 
################################################################################

#### Initialize population and disease natural history ####

results <- run_model(l_params_all)
results_noscreening <- results[['None']]

################################################################################
# Generate outputs
################################################################################

# Get prevalence, incidence, and stage outputs
l_outputs <- calc_calib_targets(l_params_all, 
                                results_noscreening, 
                                v_age_prevalence, 
                                v_age_cancer_incidence,
                                n_screen_sample)

#### Relative survival by stage and years from diagnosis ####

# Filter to individuals diagnosed with cancer in lifetime
m_cohort_cancer_dx <- results_noscreening[time_0_3 <= time_0_D] 

# Create survival object for death due to cancer
cancer_surv_obj <- with(m_cohort_cancer_dx, {
  Surv(time_3_D, fl_death_cancer)
})

# Get Kaplan-Meier fit
cancer_surv_fit = survfit(cancer_surv_obj ~ stage_dx, data = m_cohort_cancer_dx)
output_surv <- with(summary(cancer_surv_fit, times = 0:10),
                    data.frame(
                      stage = strata,
                      years_from_dx = time,
                      surv = surv
                    )) %>%
  mutate(stage = sapply(stage, as.character)) %>%
  mutate(stage = substring(stage, nchar(stage), nchar(stage)))


#### Save outputs to data folder ####
if(!debug) {
  # Lesion prevalence
  write.csv(l_outputs$prevalence, paste0(data_outpath, "prevalence_lesion.csv"), row.names = FALSE)
  
  # Cancer incidence
  write.csv(l_outputs$incidence, paste0(data_outpath, "incidence_cancer.csv"), row.names = FALSE)
  
  # Stage distribution
  write.csv(l_outputs$stage_distr, paste0(data_outpath, "stage_distr.csv"), row.names = FALSE)
  
  # Survival
  write.csv(output_surv, paste0(data_outpath, "relative_survival_cancer.csv"), row.names = FALSE)
  
  # Save true parameters in current folder
  save(param_map, file = paste0(true_param_outpath, "true_param_map.RData"))
}

