###########################  Ground truth disease simulator  ################
#
#  Objective: Simulate cancer targets, background mortality data, and relative
#  survival data for ground truth model
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
library(VGAM) # For Gompertz

###### 1.2 Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)


#### 2. General parameters ========================================================

###### 2.1 File paths
file_configs <- file.path("configs", "configs_simulated.yaml")
file_true_params <- file.path("_ground_truth", "true_params.xlsx")
file_constant_priors <- file.path("_ground_truth", "constant_priors.xlsx")
file_model_params <- file.path("_ground_truth", "params_model.rds")

###### 2.2 Configurations
# Load configs
configs <- load_configs(file_configs)

# Extract relevant parameters from configs
params_model <- configs$params_model
params_calib <- configs$params_calib
file_targets <- configs$params_calib$file_targets
file_surv <- params_model$file.surv
l_files_mort <- params_model$file.mort

###### 2.3 Other parameters
# Simulation parameters and outcome reporting
seed <- 2025 # Random seed for generating ground truth data
conf_level <- 0.95 # Confidence level for calculating outcomes
v_ages <- list( # Age ranges for calculating outcomes
  prevalence_preclin = seq(30, 80, 10),
  incidence_clin = seq(30, 90, 10),
  prevalence_lesion = seq(30, 80, 10),
  n_lesion = c(50, 80)
)
v_time_surv <- seq(0, 10) # Times from event to calculate relative survival
n_cohort <- c(screen = 10000, pop = 100000) # Number to simulate for screen vs. population samples
v_outcomes_cs <- c("prevalence", "nlesions") # Outcome types to calculate cross-sectionally
l_outcome_grps <- list( # Outcomes to calculate together for screen vs. population samples
  screen = c("prevalence_lesion", "n_lesion", "prevalence_preclin"),
  pop = c("incidence_clin", "stage_distr")
)

# Prior generation
prior_pct_width_init <- 0.8 # Percentage width of initial randomly generated prior bounds
prior_pct_multiplier <- 0.2 # Final multiplier adjustment to increase bounds of priors

# Distribution parameters for background mortality
shape_male <- 0.000078
scale_male <- 0.08
multiplier_female <- 1/1.02
max_age <- 110
n_total <- 100000

# Plot to visualize background mortality distributions
curve(pgompertz(x, shape = shape_male, scale = scale_male), 0, 110, ylim = c(0, 1)) # male
curve(pgompertz(x, shape = shape_male*multiplier_female, scale = scale_male*multiplier_female), 0, 110, ylim = c(0, 1)) # female


#### 3. Pre-processing ========================================================

# Set seed
set.seed(seed)

# Load constant priors
df_constant_priors <- read_xlsx(file_constant_priors)

# Load ground truth model parameters (with file_surv set to NULL as survival data is generated in this script)
l_params_model <- do.call(load_model_params, c(
  modifyList(params_model,
             list(file.mort = NULL,
                  file.surv = NULL),
             keep.null = T),
  list(seed = NULL,
       file.distr = file_true_params)
))

# Set distributions for background mortality
l_params_model$d_time_H_Do <- list()
l_params_model$d_time_H_Do$male <- list(
  distr = "gompertz",
  params = list(
    shape = shape_male,
    scale = scale_male
  ),
  src = "known"
)

l_params_model$d_time_H_Do$female <- list(
  distr = "gompertz",
  params = list(
    shape = shape_male*multiplier_female,
    scale = scale_male*multiplier_female
  ),
  src = "known"
)

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


#### 4. Generate population data and outputs ========================================================

###### 4.1 Simulate targets
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

# Combine calibration targets
df_target_full <- data.table()
for (target in names(l_results)) {
  # Process data based on target type
  df_target <- l_results[[target]] %>%
    rename(targets = value) %>%
    mutate(target_groups = target)
  
  if (l_params_outcome[[target]]$outcome_type %in% c("prevalence", "incidence")) {
    # Set target index and name
    df_target <- df_target %>%
      mutate(target_index = (age_start + age_end)/2,
             target_names = paste(target_groups, age_start, age_end, sep="_"))
    
    # Rename column for total n
    if (l_params_outcome[[target]]$outcome_type %in% c("incidence")) {
      df_target <- df_target %>%
        rename(n_total = person_years_total,
               n_cases = n_events)
    }
  } else if (l_params_outcome[[target]]$outcome_type %in% c("distr")) {
    df_target <- df_target %>%
      rename(target_index = stage_dx) %>%
      mutate(target_names = paste(target_groups, target_index, sep="_"))
  } else if (l_params_outcome[[target]]$outcome_type %in% c("nlesions")) {
    df_target <- df_target %>%
      mutate(age_start = v_ages[[target]][1],
             age_end = v_ages[[target]][2]) %>%
      rename(target_index = n_lesions,
             target_index_cat = n_lesions_cat) %>% 
      mutate(target_groups = target_groups,
             target_names = paste(target_groups, target_index, sep="_"))
  }
  
  # Bind to full dataframe
  df_target_full <- bind_rows(df_target_full, df_target)
}

# Subset columns for targets
df_target_full <- df_target_full %>%
  dplyr::select(any_of(c("target_names", "target_groups", 
                         "target_index", "target_index_cat",
                         "age_start", "age_end",
                         "targets", "se",
                         "ci_lb", "ci_ub", 
                         "n_cases", "n_total",
                         "sex", "lesion_type")))


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


###### 4.3 Calculate background mortality
# Create necessary variables
v_ages <- seq(0, max_age) # vector of ages
shape <- c(male = shape_male, female = shape_male*multiplier_female)
scale <- c(male = scale_male, female = scale_male*multiplier_female)

l_df_mort <- list()
for (sex in c("male", "female")) {
  # Create mortality data table
  df_mort <- data.table(
    Age = v_ages,
    pct_dead = query_distr("p", v_ages, 
                           l_params_model$d_time_H_Do[[sex]]$distr,
                           l_params_model$d_time_H_Do[[sex]]$params)
  )
  
  # Calculate probability of dying (qx) and number of survivors (lx)
  df_mort[, `:=` (qx = diff(c(pct_dead, 1)),
                  lx = n_total*(1-pct_dead))]
  
  # Remove pct_dead variable
  df_mort[, pct_dead := NULL]
  
  # Save data table
  l_df_mort[[sex]] <- df_mort
}


#### 5. Save outputs ========================================================

# Save true model parameters
saveRDS(l_params_model, file = file_model_params)

# Check if data directory exists, make if not
dir.create(dirname(file_targets), showWarnings = FALSE)

# Save targets as excel
write.csv(df_target_full, file_targets, row.names = FALSE)

# Save disease-specific relative survival from diagnosis
write.csv(output_surv, file_surv, row.names = FALSE)

# Save priors
write.csv(prior_map, file = params_calib$file_priors, row.names = FALSE)

# Save background mortality
for (sex in names(l_df_mort)) {
  # Write to csv file
  write.table(l_df_mort[[sex]], 
              file = l_files_mort[[sex]],
              row.names=FALSE)
}
