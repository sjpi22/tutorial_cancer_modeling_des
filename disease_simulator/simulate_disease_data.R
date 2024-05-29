# Disease simulator for heterogeneous disease with two precursor lesions: a and b, each with indolent (1) and aggressive (2)

################################################################################
# Setup
################################################################################
# Load packages
library(readxl)
library(data.table)
library(tidyverse)

# Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)

################################################################################
# Parameters
################################################################################
# Parameters to change
seed <- 100 # Random seed
n_cohort <- 100000 # Number to simulate in cohort
jitter_multiplier_low <- 0.5 # Jitter for default parameters - lower bound in uniform distribution
jitter_multiplier_high <- 2 # Jitter for default parameters - upper bound in uniform distribution

shape_time_2_Du <- 5.2 # Shape for gamma distribution for time from cancer onset to death without treatment
scale_time_2_Du <- 1.3 # Scale for gamma distribution for time from cancer onset to death without treatment
param_haz_detect <- exp(1) # Parameter for distribution of progression at symptomatic detection
prog_stage_mapping <- c("1" = 0, "2" = 0.24, "3" = 0.53, "4" = 0.74) # Progression at transition to cancer stage
p_surgery <- function(prog) {pmax(pmin(1, 1.2/(1+exp(12 * (prog - 0.5))) - 0.1), 0)} # Probability of receiving surgery as function of cancer progression at diagnosis (sigmoid function)
p_die_surgery <- function(age, prog) {pmax(pmin(1, 0.04 * (1 + 2*prog)/(1+exp(-0.07 * (age - 60)))), 0)} # Probability of dying during surgery as function of age and stage at diagnosis
p_cure <- function(prog) {pmax(pmin(1, 1 - exp(-0.003 * (exp(-8*(prog - 1.2)) - 1))), 0)} # Probability of cure as function of progression
cdf_time_3_Dc <- function(time, prog) {1/(1+exp(-0.5 * (time - 10 + prog * 12)))} # CDF of time from diagnosis to death as function of progression at diagnosis (if not cured)
time_3_Dc <- function(cdf, prog) {pmax(0, log(1/cdf - 1) / (-0.5) + 10 - prog * 12)} # Time from diagnosis to death from cancer (if not cured); solve for time by setting above formula equal to cdf

# Load default data
l_params_all <- load_default_params(file.incidence = NULL,
                                    file.prevalence = NULL,
                                    file.surv = NULL)

# Jitter disease distribution parameters uniformly from 0.5*default to 2*default
set.seed(seed)

# Apply jitter to parameters of cancer-related variables
vars_cancer <- grep("time_1|time_2", names(l_params_all))
for(p in names(l_params_all)[vars_cancer]) {
  params <- l_params_all[[p]]$params
  for(i in 1:length(params)) {
    if(is.numeric(params[[i]])) {
      print(params[[i]])
      l_params_all[[p]]$params[[i]] <- runif(1, jitter_multiplier_low*params[[i]], jitter_multiplier_high*params[[i]])
    }
  }
}


# Update params
l_params_all <- update_param_list(l_params_all,
                                  list(seed = seed,
                                       n_cohort = n_cohort,
                                       v_strats = l_params_all$v_strats[1]))
 

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

# Initialize matrix of patient data
m_cohort_init <- initialize_cohort(l_params_all)

# Generate baseline characteristics
simulate_baseline_data(m_cohort_init, l_params_all)

# Generate precancerous lesions
m_lesions <- simulate_lesion_data(m_cohort_init, l_params_all)

# Consolidate lesion-level data to patient-level data for cancer onset and mortality
m_cohort_cancer <- calc_cancer_onset(m_cohort_init, m_lesions, l_params_all)

# Time to death from other causes
m_cohort_cancer[, time_2_Du := query_distr("r", .N, time_2_Du$distr, time_2_Du$params)]

# Progression of disease at diagnosis
m_cohort_cancer[, prog_3s := query_distr("r", .N, prog_3s$distr, prog_3s$params)]

# Map to stage at diagnosis
m_cohort_cancer[, stage_dx := findInterval(prog_3s, prog_stage_mapping)]

# Time to and age at diagnosis
m_cohort_cancer[, time_2_3 := prog_3s * time_2_Du]
m_cohort_cancer[, time_0_3 := time_0_2 + time_2_3]

# Flag surgery and immediate death from surgery
m_cohort_cancer[, p_surgery := p_surgery(prog_3s)]
m_cohort_cancer[, fl_tx_surgery := rbinom(.N, size = 1, prob = p_surgery)]
m_cohort_cancer[, p_die_surgery := p_die_surgery(time_0_3, prog_3s)]
m_cohort_cancer[, fl_die_surgery := fl_tx_surgery * rbinom(.N, size = 1, prob = p_die_surgery)]

# Flag cure
m_cohort_cancer[, p_cure := p_cure(prog_3s)]
m_cohort_cancer[, fl_cure := rbinom(.N, size = 1, prob = p_cure)]

# Flag life extension from treatment if not cured
m_cohort_cancer[, cdf_nocure := runif(.N)]
m_cohort_cancer[, time_3_Dc_nocure := time_3_Dc(cdf_nocure, prog_3s)]

# Calculate death from cancer
m_cohort_cancer[, time_3_Dc := ifelse(fl_die_surgery, 0, ifelse(fl_cure, Inf, time_3_Dc_nocure))]
m_cohort_cancer[, time_0_Dc := time_0_3 + time_3_Dc]

# Join cancer data to patient-level data
dupe_names <- intersect(names(m_cohort_init), names(m_cohort_cancer))[-1]
m_cohort_final <- merge(m_cohort_init, m_cohort_cancer[, (dupe_names) := NULL], by = "pt_id", all.x = TRUE)

# Calculate death all-cause death and cause of death
m_cohort_final[, time_0_D := pmin(time_0_Do, time_0_Dc, na.rm = TRUE)]
m_cohort_final[, fl_death_cancer := (time_0_Do > pmin(time_0_Dc, Inf, na.rm = TRUE))]


################################################################################
# Generate outputs
################################################################################

# Lesion prevalence in screening
calc_prevalence(m_cohort_final, "time_0_1", )

# Cancer incidence by age and stage

# Relative survival by stage

# Save outputs to data folder
