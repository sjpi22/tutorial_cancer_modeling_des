# Disease simulator for heterogeneous disease with two precursor lesions: a and b
# Alternative version with data generating function different from assumed model

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
data_outpath <- 'data/alt/'
true_param_outpath <- 'disease_simulator/'
label <- 'alt'

# Check if directory exists, make if not
dir.create(file.path(data_outpath), showWarnings = FALSE)

# Model structure
if(debug) {
  n_cohort <- 1000 # Number to simulate in cohort
  n_screen_sample <- 100
} else {
  n_cohort <- 1000000 # Number to simulate in cohort
  n_screen_sample <- 500
}

# Randomization 
seed <- 8 # Random seed for generating data
jitter_multiplier_low <- 0.5 # Jitter for default parameters - lower bound in uniform distribution
jitter_multiplier_high <- 2 # Jitter for default parameters - upper bound in uniform distribution

# Probability distributions
shape_time_2_Du <- 5.2 # Shape for gamma distribution for time from cancer onset to death without treatment
scale_time_2_Du <- 1.3 # Scale for gamma distribution for time from cancer onset to death without treatment
param_haz_detect <- exp(1) # Parameter for distribution of progression at symptomatic detection
prog_stage_mapping <- c("1" = 0, "2" = 0.24, "3" = 0.53, "4" = 0.74) # Progression at transition to cancer stage
p_surgery <- function(prog) {pmax(pmin(1, 1.2/(1+exp(12 * (prog - 0.5))) - 0.1), 0)} # Probability of receiving surgery as function of cancer progression at diagnosis (sigmoid function)
p_die_surgery <- function(age, prog) {pmax(pmin(1, 0.04 * (1 + 2*prog)/(1+exp(-0.07 * (age - 60)))), 0)} # Probability of dying during surgery as function of age and stage at diagnosis
p_cure <- function(prog) {pmax(pmin(1, 1 - exp(-0.003 * (exp(-8*(prog - 1.2)) - 1))), 0)} # Probability of cure as function of progression
cdf_time_3_Dc <- function(time, prog) {1/(1+exp(-0.5 * (time - 10 + prog * 12)))} # CDF of time from diagnosis to death as function of progression at diagnosis (if not cured)
time_3_Dc <- function(cdf, prog) {pmax(0, log(1/cdf - 1) / (-0.5) + 10 - prog * 12)} # Time from diagnosis to death from cancer (if not cured); solve for time by setting above formula equal to cdf

# Outcome reporting
v_age_lesion_a <- c(seq(20, 85, 5), 100) # Age for lesion A prevalence
v_age_lesion_b <- c(50, 80) # Ages for lesion B prevalence
v_ages_prevalence <- list(a = v_age_lesion_a,
                          b = v_age_lesion_b) # List of ages for lesion-specific prevalence
v_age_cancer_incidence <- v_age_lesion_a # Age for cancer incidence 
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

# Add variables specific to disease generating process
time_2_Du <- list(distr = "gamma", 
                  params = list(shape = shape_time_2_Du, scale = scale_time_2_Du))

prog_3s <- list(distr = "hazard", 
                params = list(haz_rate = function(x) hazard_detect(x, a = param_haz_detect), 
                              haz_cumulative = function(x) cumulative_hazard_detect(x, a = param_haz_detect)),
                description = "Degree of disease progression at symptom-based diagnosis")


################################################################################
# Generate population data 
################################################################################

#### Initialize population and disease natural history ####

# Initialize matrix of patient data
m_cohort_init <- initialize_cohort(l_params_all)

# Generate baseline characteristics
simulate_baseline_data(m_cohort_init, l_params_all)

# Generate precancerous lesions
m_lesions <- simulate_lesion_data(m_cohort_init, l_params_all)

# Consolidate lesion-level data to patient-level data for cancer onset and mortality
m_cohort_cancer <- calc_cancer_onset(m_cohort_init, m_lesions, l_params_all)


#### Generate disease-related mortality ####

# Time to death from other causes
m_cohort_cancer[, time_2_Du := query_distr("r", .N, time_2_Du$distr, time_2_Du$params)]

# Progression of disease at diagnosis
m_cohort_cancer[, prog_3s := query_distr("r", .N, prog_3s$distr, prog_3s$params)]

# Map to stage at diagnosis
m_cohort_cancer[, stage_dx := findInterval(prog_3s, prog_stage_mapping)]

# Time to and age at diagnosis
m_cohort_cancer[, time_2_3 := prog_3s * time_2_Du]
m_cohort_cancer[, time_0_3 := time_0_2 + time_2_3]

# Remaining time to death if not treated
m_cohort_cancer[, time_3_Du := time_2_Du - time_2_3]

# Flag surgery and immediate death from surgery
m_cohort_cancer[, p_surgery := p_surgery(prog_3s)]
m_cohort_cancer[, fl_tx_surgery := rbinom(.N, size = 1, prob = p_surgery)]
m_cohort_cancer[, p_die_surgery := p_die_surgery(time_0_3, prog_3s)]
m_cohort_cancer[, fl_die_surgery := fl_tx_surgery * rbinom(.N, size = 1, prob = p_die_surgery)]

# Flag cure
m_cohort_cancer[, p_cure := p_cure(prog_3s)]
m_cohort_cancer[, fl_cure := rbinom(.N, size = 1, prob = p_cure)]

# Flag life extension from treatment if not cured
m_cohort_cancer[fl_cure == 0, cdf_nocure := runif(.N)]
m_cohort_cancer[fl_cure == 0, time_3_Dc_nocure := time_3_Dc(cdf_nocure, prog_3s)]

# Calculate death from cancer
m_cohort_cancer[, time_3_Dc := pmin(ifelse(fl_die_surgery, 1/52, Inf), ifelse(fl_cure, Inf, time_3_Du + time_3_Dc_nocure))]
m_cohort_cancer[, time_0_Dc := time_0_3 + time_3_Dc]


#### Summarize mortality ####

# Join cancer data to patient-level data
dupe_names <- intersect(names(m_cohort_init), names(m_cohort_cancer))[-1]
m_cohort_final <- merge(m_cohort_init, m_cohort_cancer[, (dupe_names) := NULL], by = "pt_id", all.x = TRUE)

# Calculate death all-cause death and cause of death
m_cohort_final[, time_0_D := pmin(time_0_Do, time_0_Dc, na.rm = TRUE)]
m_cohort_final[, fl_death_cancer := (time_0_Do > pmin(time_0_Dc, Inf, na.rm = TRUE))]

# Calculate survival from cancer diagnosis
m_cohort_final[time_0_3 <= time_0_D, time_3_D := time_0_D - time_0_3]
assertthat::are_equal(nrow(m_cohort_final[(fl_death_cancer == 1) & is.na(time_3_D), ]), 0) # Check that all patients with death due to cancer have time from diagnosis to death populated


################################################################################
# Generate outputs
################################################################################

# Get prevalence, incidence, and stage outputs
l_outputs <- calc_calib_targets(l_params_all, 
                                m_cohort_final, 
                                m_lesions, 
                                v_ages_prevalence, 
                                v_age_cancer_incidence,
                                n_screen_sample)

#### Relative survival by stage and years from diagnosis ####

# Filter to individuals diagnosed with cancer in lifetime
m_cohort_cancer_dx <- m_cohort_final[time_0_3 <= time_0_D] 

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
  for (lesiontype in l_params_all$v_lesions) {
    write.csv(l_outputs$prevalence[[lesiontype]], paste0(data_outpath, "prevalence_lesion_", lesiontype, ".csv"), row.names = FALSE)
  }
  
  # Cancer incidence
  write.csv(l_outputs$incidence, paste0(data_outpath, "incidence_cancer.csv"), row.names = FALSE)
  
  # Stage distribution
  write.csv(l_outputs$stage_distr, paste0(data_outpath, "stage_distr.csv"), row.names = FALSE)
  
  # Survival
  write.csv(output_surv, paste0(data_outpath, "relative_survival_cancer.csv"), row.names = FALSE)
  
  # Save true parameters in current folder
  save(param_map, file = paste0(true_param_outpath, "true_param_map_", label, ".RData"))
}
