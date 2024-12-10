# Ground truth disease simulator

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
library(assertthat)

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
n_cohort <- 100000                            # Number to simulate in cohort
n_screen_sample <- 20000                      # Number to sample for screening
seed <- 2024                                  # Random seed for generating ground truth data
v_cancer <- seq(4)                            # Cancer stages in order
rate_Cx_Dc <- c(0.05, 0.1, 0.2, 0.3)          # Exponential rate for survival from stage at diagnosis
v_param_update <- c(3, 200, 0.2, 0.06, 0.3, 0.1, 0.3, 0.2, 0.4)  # True values for unknown parameters (in order of param_map)

# Outcome reporting
v_ages_prevalence <- seq(30, 80, 10) # Ages for lesion prevalence
v_ages_incidence = c(v_ages_prevalence, 90) # Age for cancer incidence 
v_time_surv <- seq(0, 10)            # Times from event to calculate relative survival

# Prior generation
prior_pct_width_init <- 0.3 # Percentage width of initial randomly generated prior bounds
prior_pct_multiplier <- 0.2 # Final multiplier adjustment to increase bounds of priors

# Define parameters for calculating outcomes
l_outcome_params <- list(
  prevalence = list(
    outcome_name = "prevalence",
    m_patients = "m_patients",
    start_var = "time_H_P", 
    end_var = "time_H_C", 
    censor_var = "time_screen_censor", 
    v_ages = v_ages_prevalence,
    n_screen_sample = n_screen_sample),
  incidence = list(
    outcome_name = "incidence",
    m_patients = "m_patients",
    time_var = "time_H_C", 
    censor_var = "time_H_D", 
    v_ages = v_ages_incidence),
  stage_distr = list(
    outcome_name = "distr",
    m_patients = "m_patients",
    grouping_var = "stage_dx", 
    event_var = "time_H_C", 
    censor_var = "time_H_D",
    groups_expected = 1:4)
)

l_censor_vars <- list(
  time_screen_censor = c("time_H_C", "time_H_D")
)

#### Load and update parameters ####

# Load default data
l_params_all <- load_model_params(
  seed = seed,
  n_cohort = n_cohort,
  v_cancer = v_cancer,
  file.surv = NULL
)

# Add ground truth distribution for survival from diagnosis
for (i in v_cancer) {
  l_params_all$d_time_C_Dc[[i]]$distr <- "exp"
  l_params_all$d_time_C_Dc[[i]]$params$rate = rate_Cx_Dc[i]
}

# Map variables to parameters for tuning - make dataframe of all parameters with "src = unknown"
param_map <- make_param_map(l_params_all)

# Update params
l_params_all <- update_param_from_map(l_params_all, v_param_update, param_map)
param_map$param_val <- v_param_update

# Set seed
set.seed(l_params_all$seed)

################################################################################
# Generate population data and outputs
################################################################################

# Simulate cohort
m_cohort <- run_base_model(l_params_all)

# Establish priors by adding random noise around the true parameter values
prior_map <- param_map %>%
  mutate(shift = runif(nrow(param_map))) %>%
  mutate(distr = "unif",
         min = param_val * (1 - prior_pct_width_init/2 + (shift - 0.5) * prior_pct_width_init),
         max = param_val * (1 + prior_pct_width_init/2 + (shift - 0.5) * prior_pct_width_init)) %>%
  mutate(min = min * (1 - prior_pct_multiplier),
         max = max * (1 + prior_pct_multiplier)) %>%
  dplyr::select(-c("param_val", "shift"))

# Get prevalence, incidence, and stage outputs
l_outputs <- calc_calib_outputs(m_cohort, 
                                l_outcome_params = l_outcome_params,
                                l_censor_vars = l_censor_vars)


#### Relative survival by stage and years from diagnosis ####

# Filter to individuals diagnosed with cancer in lifetime
m_cohort_cancer_dx <- m_cohort[time_H_C < time_H_D] 

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
# Preclinical cancer prevalence
write.csv(l_outputs$prevalence %>%
            rename(targets = value) %>%
            dplyr::select(-c("person_years_cases", "person_years_total")), 
          file.path(data_outpath, "prevalence_asymptomatic_cancer.csv"), row.names = FALSE)

# Cancer incidence
write.csv(l_outputs$incidence %>%
            rename(targets = value) %>%
            dplyr::select(-c("total_atrisk", "n_events")), 
          file.path(data_outpath, "incidence_symptomatic_cancer.csv"), row.names = FALSE)

# Stage distribution
write.csv(l_outputs$stage_distr %>%
          rename(targets = value), 
          file.path(data_outpath, "stage_distr.csv"), row.names = FALSE)

# Survival
write.csv(output_surv, file.path(data_outpath, "relative_survival_cancer.csv"), row.names = FALSE)

# Save priors
saveRDS(prior_map, file = file.path(data_outpath, "priors.rds"))

# Save true parameters in current folder
saveRDS(param_map, file = file.path(true_param_outpath, "true_param_map.rds"))
