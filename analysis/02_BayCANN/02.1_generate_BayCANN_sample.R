###########################  Generate BayCANN Sample  ##########################
#
#  Objective: Program to simulate parameter inputs and model outputs for 
#  BayCANN model calibration
########################### <<<<<>>>>> #########################################

rm(list = ls()) # Clean environment
options(scipen = 999) # View data without scientific notation

#### 1.Libraries and functions  ==================================================

###### 1.1 Load packages
library(data.table)
library(tidyverse)
library(lhs)
library(doParallel)
library(foreach)
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
file_params_calib <- configs$paths$file_params_calib
file_plot_labels <- configs$paths$file_plot_labels
alpha_ci <- 1 - configs$params_model$conf_level

# Get list of relevant output file paths and load to global environment
l_filepaths <- update_config_paths("files_baycann", configs$paths)
list2env(l_filepaths, envir = .GlobalEnv)

# Load BayCANN and coverage analysis parameters from configs file
list2env(configs$params_baycann$params_sampling, envir = .GlobalEnv)
list2env(configs$params_coverage, envir = .GlobalEnv)


#### 3. Pre-processing  ===========================================

# Load model and calibration parameters
l_params_calib <- readRDS(file_params_calib)

# Load plot labels
df_plot_labels <- read.csv(file_plot_labels)

# Set seed
set.seed(l_params_calib$l_params_model$seed, kind = "L'Ecuyer-CMRG")

# Set number of cores to use
if(!is.na(as.integer(Sys.getenv("SLURM_NTASKS_PER_NODE")))) {
  # If using Sherlock, use environment variable to set the number of cores to use
  registerDoParallel(cores = (Sys.getenv("SLURM_NTASKS_PER_NODE")))
  print("Running on Sherlock")
} else {
  # If running locally, use all available cores except for reserved ones
  registerDoParallel(cores = detectCores(logical = TRUE) - l_params_calib$n_cores_reserved_local)
  print("Running locally")
}

# Show the number of parallel workers to be used
print(paste("# parallel workers:", getDoParWorkers())) 


#### 4. Generate inputs and outputs  ===========================================

# Generate parameter sample with LHS
m_param_samp <- lhs_param_samp(prior_map = l_params_calib$prior_map, 
                               n_samp_per_param = n_samp_per_param)

# Run model for each input parameter sample and get corresponding outputs with parallel processing
stime <- system.time({
  m_outputs <- foreach(
    i=1:nrow(m_param_samp), 
    .combine=rbind, 
    .inorder=TRUE, 
    .packages=c("data.table","tidyverse")) %dopar% {
      # Get row of parameters and calculate outputs
      v_params_update <- m_param_samp[i,]
      v_calib_outputs <- with(l_params_calib, {
        params_to_outputs(
          l_params_model = l_params_model,
          v_params_update = v_params_update,
          param_map = prior_map,
          l_params_outcome = l_params_outcome,
          l_censor_vars = l_censor_vars
        )
      })
      # Call item to save
      t(v_calib_outputs)
    }
})
print(stime)
closeAllConnections()
  
# Set column names
colnames(m_outputs) <- l_params_calib$df_target$target_names

# Check for any NaN
validate_that(
  sum(sapply(m_param_samp, function(x) any(is.nan(x)))) == 0, 
  msg = "Parameters include NaN")
validate_that(
  sum(sapply(m_outputs, function(x) any(is.nan(x)))) == 0, 
  msg = "Outputs include NaN")

# Save parameter sample and corresponding outputs
saveRDS(list(m_param_samp = m_param_samp, 
             m_calib_outputs = m_outputs,
             runtime = stime), 
        file = file_sample)


#### 5. Generate inputs and outputs  ===========================================

# Load BayCANN sample
BayCANN_sample <- readRDS(file_sample)

# Process targets
df_targets <- l_params_calib$df_target %>%
  mutate(ci_lb = ifelse(is.na(ci_lb), targets - se*qnorm(1 - alpha_ci/2), ci_lb),
         ci_ub = ifelse(is.na(ci_ub), targets + se*qnorm(1 - alpha_ci/2), ci_ub)) %>% # Create plot labels
  left_join(df_plot_labels, by = "target_groups")
df_targets$plot_grps <- factor(df_targets$plot_grps, levels = df_plot_labels$plot_grps)

# Plot and save coverage
plt_coverage <- plot_coverage(df_targets = df_targets, 
                              m_outputs = BayCANN_sample$m_calib_outputs, 
                              file_fig_coverage = file_fig_coverage,
                              target_range = "ci",
                              plt_size_text = plt_size_text,
                              labeller_multiplier = 6)
plt_coverage
