###########################  Ground truth disease simulator  ################
#
#  Objective: Simulate targets for preclinical cancer prevalence and cancer 
#  incidence with ground truth model
########################### <<<<<>>>>> ##############################################

rm(list = ls()) # Clean environment
options(scipen = 999) # View data without scientific notation

#### 1.Libraries and functions  ==================================================

###### 1.1 Load packages
library(readxl)
library(data.table)
library(tidyverse)
library(survival)
library(assertthat)

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
params_model <- configs$params_model
params_calib <- configs$params_calib

###### 2.2 File paths
file_true_params <- file.path("_ground_truth", "true_params.xlsx")
file_constant_priors <- file.path("_ground_truth", "constant_priors.xlsx")

###### 2.3 Other parameters
# Simulation parameters and outcome reporting
seed <- 2025 # Random seed for generating ground truth data
conf_level <- 0.95 # Confidence level for calculating outcomes
v_ages <- list( # Age ranges for calculating outcomes
  prevalence = seq(30, 80, 10),
  incidence = seq(30, 90, 10),
  prevalence_lesion = seq(30, 80, 10)
)
v_time_surv <- seq(0, 10) # Times from event to calculate relative survival
n_cohort <- c(screen = 10000, pop = 100000) # Number to simulate for screen vs. population samples
v_outcomes_cs <- c("prevalence", "nlesions") # Outcome types to calculate cross-sectionally
l_outcome_grps <- list( # Outcomes to calculate together for screen vs. population samples
  screen = c("prevalence_lesion", "n_lesions", "prevalence"),
  pop = c("incidence", "stage_distr")
)

# Prior generation
prior_pct_width_init <- 0.8 # Percentage width of initial randomly generated prior bounds
prior_pct_multiplier <- 0.2 # Final multiplier adjustment to increase bounds of priors


#### 3. Pre-processing ========================================================

# Set seed
set.seed(seed)

# Load constant priors
df_constant_priors <- read_xlsx(file_constant_priors)

# Load ground truth model parameters (with file_surv set to NULL as survival data is generated in this script)
l_params_model <- do.call(load_model_params, c(
  modifyList(params_model,
             list(file.surv = NULL),
             keep.null = T),
  list(seed = NULL,
       file.distr = file_true_params)
))

# Map variables to parameters for tuning - make dataframe of all parameters with "src = unknown"
param_map <- make_param_map(l_params_model)

# Establish priors by adding random noise around the true parameter values
prior_map <- param_map %>%
  mutate(shift = runif(nrow(param_map))) %>%
  mutate(distr = "unif",
         min = param_val * (1 - prior_pct_width_init/2 + (shift - 0.5) * prior_pct_width_init),
         max = param_val * (1 + prior_pct_width_init/2 + (shift - 0.5) * prior_pct_width_init)) %>%
  mutate(min = round(min * (1 - prior_pct_multiplier), 2),
         max = round(max * (1 + prior_pct_multiplier), 2)) %>%
  dplyr::select(-c("param_val", "shift")) %>%
  # Bind constant priors
  bind_rows(df_constant_priors) %>%
  relocate(idx, .after = var_name)

# Extract outcome parameters
l_params_outcome <- params_calib$l_params_outcome
l_censor_vars <- params_calib$l_censor_vars

# Process outcome parameters
for (target in names(l_params_outcome)) {
  # Modify outcome list to output uncertainty with assigned confidence level
  l_params_outcome[[target]][["lit_params"]][["output_uncertainty"]] <- TRUE
  l_params_outcome[[target]][["lit_params"]][["conf_level"]] <- params_model$conf_level
  
  # Modify outcome list to change to cross-sectional functions if not already
  if (l_params_outcome[[target]][["outcome_type"]] %in% v_outcomes_cs) {
    l_params_outcome[[target]][["lit_params"]][["method"]] <- "cs"
  }
  
  # If applicable, add age ranges for calculating outcomes
  l_params_outcome[[target]]$lit_params$v_ages <- v_ages[[target]]
}

# Check if data directory exists, make if not
dir.create(dirname(params_calib$l_params_outcome[[1]]$file_path), showWarnings = FALSE)


#### 4. Generate population data and outputs ========================================================

###### 4.1 Simulate data and outputs
# Calculate outcomes from different simulated cohorts
l_cohorts <- list()
l_results <- list()
for (grp in names(l_outcome_grps)) {
  # Update cohort size
  l_params_model$n_cohort <- n_cohort[grp]
  
  # Calculate outputs in group
  l_results_grp <- params_to_outputs(l_params_model = l_params_model, 
                                     l_params_outcome = l_params_outcome[l_outcome_grps[[grp]]], 
                                     l_censor_vars = l_censor_vars,
                                     reshape_output = FALSE, 
                                     individual_data = TRUE,
                                     conf_level = conf_level)
  
  # Append results
  l_cohorts[[grp]] <- l_results_grp$m_cohort
  l_results <- c(l_results, l_results_grp$outputs)
}


###### 4.2 Calculate relative survival by stage and years from diagnosis
# Set patient-level data
m_cohort <- l_cohorts$pop
if (is.data.table(m_cohort)) {
  m_patients <- m_cohort
} else {
  m_patients <- m_cohort$patient_level
}

# Filter to individuals diagnosed with cancer in lifetime
m_cohort_cancer_dx <- m_patients[time_H_C < time_H_D] 

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


#### 5. Save outputs ========================================================

# Save calibration targets
for (target in names(l_results)) {
  # Process data based on target type
  df_target <- l_results[[target]] %>%
    rename(targets = value) 
  
  # Save data to file path
  write.csv(df_target,
            l_params_outcome[[target]][["file_path"]],
            row.names = FALSE)
}

# Save disease-specific relative survival from diagnosis
write.csv(output_surv, params_model$file.surv, row.names = FALSE)

# Save priors
write.csv(prior_map, file = params_calib$file_priors, row.names = FALSE)
