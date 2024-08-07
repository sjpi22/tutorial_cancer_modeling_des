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

# Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)

################################################################################
# Parameters
################################################################################

# File paths
data_outpath <- 'data'
true_param_outpath <- 'ground_truth'

# Check if directory exists, make if not
dir.create(file.path(data_outpath), showWarnings = FALSE)

# Model parameters
n_cohort <- 1000000                           # Number to simulate in cohort
seed <- 1                                     # Random seed for generating data
v_param_update <- c(3, 300, 0.5, 0.42, 0.25)  # True values for unknown parameters (in order of param_map)
rate_C1_Dc <- 0.1                             # Rate for survival from early-stage diagnosis exponential distribution
rate_C2_Dc <- 0.3                             # Rate for survival from late-stage diagnosis exponential distribution

# Outcome reporting
v_ages_prevalence <- seq(30, 80, 10) # Ages for lesion prevalence
v_ages <- list(prevalence = v_ages_prevalence,
               incidence = c(v_ages_prevalence, 100)) # Age for cancer incidence 
v_time_surv <- seq(0, 10)            # Times from event to calculate relative survival


#### Load and update parameters ####

# Load default data
l_params_all <- load_default_params(file.surv = NULL)

# Add true survival distribution of exponential from diagnosis
l_params_all$time_C1_Dc$distr <- "exp"
l_params_all$time_C1_Dc$params <- list(rate = rate_C1_Dc)
l_params_all$time_C2_Dc$distr <- "exp"
l_params_all$time_C2_Dc$params <- list(rate = rate_C2_Dc)

# Map variables to parameters for tuning - make dataframe of all parameters with "src = unknown"
param_map <- make_param_map(l_params_all)

# Update params
l_params_all <- update_param_from_map(l_params_all, v_param_update, param_map)
l_params_all <- update_param_list(l_params_all,
                                  list(seed = seed,
                                       n_cohort = n_cohort,
                                       v_strats = l_params_all$v_strats[1]))
param_map$param_val <- v_param_update

################################################################################
# Generate population data 
################################################################################

results <- run_model(l_params_all)
results_noscreening <- results[['None']]

################################################################################
# Generate outputs
################################################################################

# Establish priors
prior_map <- param_map %>%
  mutate(shift = runif(nrow(param_map))) %>%
  mutate(distr = "unif",
         min = param_val * (0.7 + (shift - 0.5) * 0.6),
         max = param_val * (1.3 + (shift - 0.5) * 0.6)) %>%
  dplyr::select(-c("param_val", "shift"))

# Get prevalence, incidence, and stage outputs
l_outputs <- calc_calib_targets(l_params_all, 
                                results_noscreening, 
                                v_ages)

#### Relative survival by stage and years from diagnosis ####

# Filter to individuals diagnosed with cancer in lifetime
m_cohort_cancer_dx <- results_noscreening[time_H_C < time_H_D] 

# Create survival object for death due to cancer
cancer_surv_obj <- with(m_cohort_cancer_dx, {
  Surv(time_C_D, fl_Dc)
})

# Get Kaplan-Meier fit
cancer_surv_fit = survfit(cancer_surv_obj ~ stage_dx, data = m_cohort_cancer_dx)
output_surv <- with(summary(cancer_surv_fit, times = v_time_surv),
                    data.frame(
                      stage = strata,
                      years_from_dx = time,
                      surv = surv
                    )) %>%
  mutate(stage = sapply(stage, as.character)) %>%
  mutate(stage = substring(stage, nchar(stage), nchar(stage)))


#### Save outputs to data folder ####
# Lesion prevalence
write.csv(l_outputs$prevalence, file.path(data_outpath, "prevalence_asymptomatic_cancer.csv"), row.names = FALSE)

# Cancer incidence
write.csv(l_outputs$incidence, file.path(data_outpath, "incidence_symptomatic_cancer.csv"), row.names = FALSE)

# Stage distribution
write.csv(l_outputs$stage_distr, file.path(data_outpath, "stage_distr.csv"), row.names = FALSE)

# Survival
write.csv(output_surv, file.path(data_outpath, "relative_survival_cancer.csv"), row.names = FALSE)

# Save priors
saveRDS(prior_map, file = file.path(data_outpath, "priors.rds"))

# Save true parameters in current folder
saveRDS(param_map, file = file.path(true_param_outpath, "true_param_map.rds"))
