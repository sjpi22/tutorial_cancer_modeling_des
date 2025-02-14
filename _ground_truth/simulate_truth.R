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

###### 2.2 File paths
path_truth <- "_ground_truth"
path_data <- paths$data
file_distr <- file.path(path_truth, "true_params.xlsx")

# Check if data directory exists, make if not
dir.create(file.path(path_data), showWarnings = FALSE)

###### 2.3 Other parameters
# Simulation parameters and outcome reporting
seed <- 2025 # Random seed for generating ground truth data
n_screen_sample <- 20000 # Subset of cohort for prevalence outputs
v_time_surv <- seq(0, 10) # Times from event to calculate relative survival

# Prior generation
prior_pct_width_init <- 0.8 # Percentage width of initial randomly generated prior bounds
prior_pct_multiplier <- 0.2 # Final multiplier adjustment to increase bounds of priors


#### 3. Pre-processing ========================================================

# Set seed
set.seed(seed)

# Load ground truth model parameters (with file_surv set to NULL as survival data is generated in this script)
l_params_model <- do.call(load_model_params, c(
  modifyList(params_model,
             list(file.surv = NULL),
             keep.null = T),
  list(seed = seed,
       file.distr = file_distr)
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
  dplyr::select(-c("param_val", "shift"))

# Update screening sample for prevalence
for (target in names(params_calib$l_outcome_params)) {
  if (params_calib$l_outcome_params[[target]]["outcome_type"] == "prevalence") {
    params_calib$l_outcome_params[[target]] <- c(params_calib$l_outcome_params[[target]], 
                                                 list(n_screen_sample = n_screen_sample))
  }
}


#### 4. Generate population data and outputs ========================================================

###### 4.1 Simulate data and outputs
# Simulate cohort
m_cohort <- run_base_model(l_params_model)

# Separate patient and lesion data as necessary
if (params_model$lesion_state == T) {
  m_patients <- m_cohort$patient_level
  m_lesions <- m_cohort$lesion_level
} else {
  m_patients <- m_cohort$patient_level
  m_lesions <- NULL
}

# Get calibration outputs
l_outputs <- calc_calib_outputs(m_patients, 
                                l_outcome_params = params_calib$l_outcome_params,
                                l_censor_vars = params_calib$l_censor_vars,
                                m_lesions = m_lesions)

###### 4.2 Calculate relative survival by stage and years from diagnosis
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
for (target in names(l_outputs)) {
  # Process data based on target type
  df_target <- l_outputs[[target]] %>%
    rename(targets = value) 
  
  if (params_calib$l_outcome_params[[target]][["outcome_type"]] == "prevalence") {
    df_target <- df_target %>%
      dplyr::select(-c("person_years_cases", "person_years_total"))
  } else if (params_calib$l_outcome_params[[target]][["outcome_type"]] == "incidence") {
    df_target <- df_target %>%
      dplyr::select(-c("total_atrisk", "n_events"))
  }
  
  # Save data to file path
  write.csv(df_target,
            params_calib$l_outcome_params[[target]][["file_path"]],
            row.names = FALSE)
}

# Save disease-specific relative survival from diagnosis
write.csv(output_surv, params_model$file.surv, row.names = FALSE)

# Save priors
write.csv(prior_map, file = params_calib$file_priors, row.names = FALSE)
