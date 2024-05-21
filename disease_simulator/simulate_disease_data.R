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
seed <- 100
n_cohort <- 100
jitter_multiplier_low <- 0.5
jitter_multiplier_high <- 2

# Load default data
l_params_all <- load_default_params(file.incidence = NULL,
                                    file.prevalence = NULL,
                                    file.surv = NULL)

# Jitter disease distribution parameters uniformly from 0.5*default to 2*default
set.seed(seed)
df_vars <- l_params_all$df_vars %>%
  filter(vargroup != "cancer") %>% # Remove cancer, since we will calculate cancer parameters differently here
  filter(distr != "direct") # Remove directly calculated variables

# Add cancer parameters
df_vars <- add_row(
  df_vars,
  varname = "time_2_Du", 
  description = "Time from cancer onset to death if untreated",
  vargroup = "cancer", 
  distr = "gamma", 
  params = list(list(shape = 4, scale = 4))
) 

# Apply jitter to current parameters
updated_params <- lapply(df_vars$params, function(p) {
  for(i in 1:length(p)) {
    if(is.numeric(p[[i]])) {
      p[[i]] = runif(1, jitter_multiplier_low*p[[i]], jitter_multiplier_high*p[[i]])
    }
  }
  return(p)
})

df_vars$params[grep("lesion|cancer", df_vars$vargroup)] <- updated_params[grep("lesion|cancer", df_vars$vargroup)]

# Update params
l_params_all <- update_param_list(l_params_all,
                                  list(seed = seed,
                                       n_cohort = n_cohort,
                                       v_strats = l_params_all$v_strats[1]))
l_params_all$df_vars <- df_vars
  
# Generate data 
data <- run_model(l_params_all)

m_cohort <- data$None$m_cohort
m_lesion <- data$None$m_lesion

################################################################################
# Simulate incidence of precursor lesions by age
################################################################################

# Add
add_row(varname = "prog_3s", 
        description = "Degree of disease progression at symptom-based diagnosis",
        vargroup = "cancer", 
        distr = "hazard", 
        params = list(list(haz_rate = function(x) hazard_detect(x, a = exp(1)), 
                           haz_cumulative = function(x) cumulative_hazard_detect(x, a = exp(1))))) %>%
  add_row(varname = "time_2_3", 
          description = "Time from cancer onset to diagnosis",
          vargroup = "cancer", 
          distr = "direct", 
          params = list("time_2_Du * prog_3s")) %>%
  add_row(varname = "time_0_3", 
          description = "Time from birth to cancer diagnosis",
          vargroup = "cancer", 
          distr = "direct", 
          params = list("time_0_2 * time_2_3")) %>%
  add_row(varname = "p_death_surgery",
          description = "Probability of death during surgery",
          vargroup = "cancer",
          distr = "direct",
          params = list()) %>%
  add_row(varname = "fl_death_surgery",
          description = "Death at surgery",
          vargroup = "cancer",
          distr = "binom",
          params = list(list(size = 1, prob = "p_death_surgery")))

################################################################################
# Simulate cancer onset by age
################################################################################

################################################################################
# Simulate time to death without treatment
################################################################################

################################################################################
# Simulate hazard of symptomatic detection
################################################################################


################################################################################
# Simulate survival after treatment
################################################################################

# Outputs: lesion prevalence, cancer incidence, survival