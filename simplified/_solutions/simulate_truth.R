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
  n_screen_sample <- 10000
}

# Randomization 
seed <- 1 # Random seed for generating data

# Outcome reporting
v_ages_prevalence <- c(seq(30, 80, 10), 100) # Age for lesion prevalence
v_ages <- list(prevalence = v_ages_prevalence,
               incidence = v_ages_prevalence) # Age for cancer incidence 
v_time_surv <- seq(0, 10) # Times from event to calculate relative survival


#### Load and update parameters ####

# Load default data
l_params_all <- load_default_params(file.surv = NULL)

# Add true survival distribution of exponential from diagnosis
l_params_all$time_2i_Dc$distr <- "exp"
l_params_all$time_2i_Dc$params <- list(rate = 0.1)
l_params_all$time_2ii_Dc$distr <- "exp"
l_params_all$time_2ii_Dc$params <- list(rate = 0.3)

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

# Establish priors
prior_map <- param_map %>%
  mutate(shift = runif(nrow(param_map))) %>%
  mutate(prior_distr = "unif",
         prior_min = param_val * (0.7 + (shift - 0.5) * 0.6),
         prior_max = param_val * (1.3 + (shift - 0.5) * 0.6)) %>%
  dplyr::select(-c("param_val", "shift"))

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
                                v_ages,
                                n_screen_sample)

#### Relative survival by stage and years from diagnosis ####

# Filter to individuals diagnosed with cancer in lifetime
m_cohort_cancer_dx <- results_noscreening[time_0_2 <= time_0_D] 

# Create survival object for death due to cancer
cancer_surv_obj <- with(m_cohort_cancer_dx, {
  Surv(time_2_D, fl_death_cancer)
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
  write.csv(l_outputs$prevalence, paste0(data_outpath, "prevalence_asymptomatic_cancer.csv"), row.names = FALSE)
  
  # Cancer incidence
  write.csv(l_outputs$incidence, paste0(data_outpath, "incidence_symptomatic_cancer.csv"), row.names = FALSE)
  
  # Stage distribution
  write.csv(l_outputs$stage_distr, paste0(data_outpath, "stage_distr.csv"), row.names = FALSE)
  
  # Survival
  write.csv(output_surv, paste0(data_outpath, "relative_survival_cancer.csv"), row.names = FALSE)
  
  # Save priors
  save(prior_map, file = paste0(data_outpath, "priors.RData"))
  
  # Save true parameters in current folder
  save(param_map, file = paste0(true_param_outpath, "true_param_map.RData"))
}

