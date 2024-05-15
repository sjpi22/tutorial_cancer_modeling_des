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
# n_cohort <- 1e7
n_cohort <- 100
jitter_multiplier_low <- 0.5
jitter_multiplier_high <- 2

# Load default data
l_params_all <- load_default_params(file.incidence = NULL,
                                    file.prevalence = NULL,
                                    file.surv = NULL)

# Jitter disease distribution parameters uniformly from 0.5*default to 2*default
set.seed(seed)
df_vars <- l_params_all$df_vars

updated_params <- lapply(df_vars$params, function(p) {
  for(i in 1:length(p)) {
    print(p)
    if(is.numeric(p[[i]])) {
      p[[i]] = runif(1, jitter_multiplier_low*p[[i]], jitter_multiplier_high*p[[i]])
    }
  }
  return(p)
})

df_vars$params[df_vars$vargroup %in% c("lesion", "cancer")] <- updated_params[df_vars$vargroup %in% c("lesion", "cancer")]

# Update params
l_params_all <- update_param_list(l_params_all,
                                  list(seed = seed,
                                       n_cohort = n_cohort,
                                       v_strats = l_params_all$v_strats[1],
                                       df_vars = df_vars))
  
# Generate data 
data <- run_model(l_params_all)

################################################################################
# Simulate incidence of precursor lesions by age
################################################################################

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