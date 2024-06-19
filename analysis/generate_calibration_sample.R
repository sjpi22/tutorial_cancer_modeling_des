# Program to generate model outputs with various inputs for model calibration

# Load packages
library(readxl)
library(data.table)
library(tidyverse)
library(lhs)

# Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)

################################################################################
###  Prepare parameters for calibration                                      ###
################################################################################

# Variables
n_samp <- 1000

# Load default data
l_params_all <- load_default_params()

# Make cohort small for testing
l_params_all <- update_param_list(l_params_all,
                                  list(n_cohort = 10000,
                                       v_strats = l_params_all$v_strats[1]))

# Set seed
set.seed(l_params_all$seed)

# Get all variables to tune
vars <- names(l_params_all)[as.logical(pmax(startsWith(names(l_params_all), "time_"),
                                            startsWith(names(l_params_all), "n_add"),
                                            startsWith(names(l_params_all), "size_"),
                                            startsWith(names(l_params_all), "p_time_"),
                                            startsWith(names(l_params_all), "stage_"),
                                            startsWith(names(l_params_all), "time_")))]

# Remove death from other causes
vars <- vars[!startsWith(vars, "time_0_Do")]

# Map variables to parameters for tuning - make dataframe of all parameters with "src = unknown"
param_map <- lapply(names(l_params_all), function(x) {
  # Find distribution variables, which have 'params' item
  var <- l_params_all[[x]]
  if ('params' %in% names(var)) {
    # Only keep distributions that are unknown (not from data or assumption)
    if (var$src == "unknown") {
      # Get params
      par <- var$params
      
      # Consolidate params and names
      if (var$distr == 'empirical' & 'probs' %in% names(par)) {
        # For empirical probabilities, save probabilities as calibrate-able value
        df <- data.frame(var_name = x,
                         param_name = 'probs',
                         param_index = par$xs,
                         param_val = par$probs)
        
        # Remove last row (whose value is 1 - sum of other values)
        return(df[-nrow(df),])
      } else {
        df <- data.frame(var_name = x,
                         param_name = names(par),
                         param_index = 1,
                         param_val = unname(unlist(par)))
        
        return(df)
      }
    }
  }
})
param_map <- rbindlist(param_map) %>%
  # Set variable identifier
  mutate(var_id = paste(var_name, param_name, param_index, sep = '.'))

# Set prior distributions (set uniform for everything)
param_map$prior_distr <- "unif"
param_map$prior_min <- 1/3 * param_map$param_val
param_map$prior_max <- 3 * param_map$param_val


# Get number of params to calibrate
n_param <- nrow(param_map)
  
################################################################################
###  Generate a random sample of input values                                ###
################################################################################

# Sample unit Latin Hypercube
m_lhs_unit <- randomLHS(n_samp, n_param)

# Rescale to min/max of each parameter
m_param_samp <- matrix(nrow = n_samp, ncol = n_param)
for (i in 1:n_param) {
  m_param_samp[, i] <- qunif(m_lhs_unit[, i],
                             min = param_map$prior_min[i],
                             max = param_map$prior_max[i])
}
colnames(m_param_samp) <- param_map$var_id

################################################################################
###  Generate corresponding outputs                                          ###
################################################################################
# Run model for each input parameter sample and get corresponding targets
for (i in 1:n_samp) {
  v_params_update <- m_param_samp[i,]
  v_calib_targets <- params_to_calib_targets(l_params_all, v_params_update, param_map)
}

# Save data files
save(param_map, m_param_samp, file = 'data/calibration_sample.RData')
