###########################  Load Calibration Parameters  ##########################
#
#  Objective: Program to load general calibration parameters, perform Monte 
#  Carlo error analysis for sample size, and check coverage of targets
########################### <<<<<>>>>> #########################################

rm(list = ls()) # Clean environment
options(scipen = 999) # View data without scientific notation

#### 1.Libraries and functions  ==================================================

###### 1.1 Load packages
library(readxl)
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
params_model <- configs$params_model
params_calib <- configs$params_calib
file_params_calib <- configs$paths$file_params_calib
file_plot_labels <- configs$paths$file_plot_labels

# Get list of relevant output file paths and load to global environment
l_filepaths <- update_config_paths("files_coverage", configs$paths)
list2env(l_filepaths, envir = .GlobalEnv)

# Load Monte Carlo and coverage analysis parameters from configs file
list2env(configs$params_montecarlo, envir = .GlobalEnv)
list2env(configs$params_coverage, envir = .GlobalEnv)

###### 2.2 Other parameters
v_outcomes_cs <- c("prevalence", "nlesions") # Outcome types that can be calculated cross-sectionally


#### 3. Pre-processing  ===========================================

# Load model parameters
l_params_init <- do.call(load_model_params, c(
  params_model
))

# Load calibration parameters
l_params_calib <- do.call(load_calib_params, c(
  l_params_model = list(l_params_init),
  params_calib
))

# Load plot labels
df_plot_labels <- read.csv(file_plot_labels)

# Create cross-sectional version of outcome parameters
l_params_outcome_cs <- l_params_calib$l_params_outcome
v_labels_cs <- c()
for (target in names(l_params_outcome_cs)) {
  if (l_params_outcome_cs[[target]][["outcome_type"]] %in% v_outcomes_cs) {
    # Change to cross-sectional function
    l_params_outcome_cs[[target]][["outcome_type"]] <- paste0(l_params_outcome_cs[[target]][["outcome_type"]], "_cs")
    
    # Keep track of outcome labels that can be cross-sectional
    v_labels_cs <- c(v_labels_cs, target)
  }
}

# Create outcome list to compare longitudinal and cross-sectional calculations
l_params_outcomes_cs_renamed <- l_params_outcome_cs[v_labels_cs]
names(l_params_outcomes_cs_renamed) <- paste0(names(l_params_outcomes_cs_renamed), "_cs")

# Set number of cores to use
registerDoParallel(cores = detectCores(logical = TRUE) - l_params_calib$n_cores_reserved_local)
print(getDoParWorkers())


#### 4. Monte Carlo error analysis  ===========================================

# Set Monte Carlo simulation parameters
l_params_mc <- l_params_calib$l_params_model
l_params_mc$n_cohort <- n_init
l_params_mc$seed <- NULL # Remove seed that is reset every time model is run
set.seed(l_params_calib$l_params_model$seed, kind = "L'Ecuyer-CMRG") # Set seed externally; kind set for parallel package

# Test sample size with longitudinally vs. cross-sectionally calculated vs. other outcomes
stime_mc <- system.time({
  df_res_mc <- foreach(
    i=1:n_mc_reps, 
    .combine=rbind, 
    .inorder=FALSE, 
    .packages=c("data.table","tidyverse")) %dopar% {
      # Calculate outputs in any order
      v_outputs_mc <- params_to_outputs(l_params_mc, 
                                        l_params_outcome = c(l_params_calib$l_params_outcome, l_params_outcomes_cs_renamed),
                                        l_censor_vars = l_params_calib$l_censor_vars)
      
      # Call item to save
      t(v_outputs_mc)
    }
})
print(stime_mc)
  
# Compare mean, SD, and required cohort size for longitudinal vs. cross-sectional outcomes
df_sample <- data.frame(
  rbind(l_params_calib$df_targets,
        l_params_calib$df_targets[l_params_calib$df_targets$target_groups %in% v_labels_cs, ]),
  mean_mc = colMeans(df_res_mc),
  sd_mc = apply(df_res_mc, 2, sd)
) %>%
  mutate( # Calculate required cohort size
    n_target = n_init * (sd_mc / se)^2,
    label = c(rep("long", nrow(l_params_calib$df_targets)), rep("cs", sum(l_params_calib$df_targets$target_groups %in% v_labels_cs))),
    cs_option = target_groups %in% v_labels_cs
  )

# Get maximum cohort size per target group and calculation type
df_sample_max <- df_sample %>%
  group_by(target_groups, label, cs_option) %>%
  summarise(n_max = max(n_target), .groups = "drop") %>%
  arrange(desc(n_max))

# Set final N as multiple of maximum required cohort size, rounded up to second highest digit
n_final <- mc_multiplier * max(df_sample_max$n_max)
n_digits <- floor(log(n_final, base = 10))
n_final <- ceiling(n_final / 10^(n_digits - 1)) * 10^(n_digits - 1)

# Update calibration parameters
l_params_calib$l_params_model$n_cohort <- n_final

# If the faster cross-sectional calculations do not require a larger cohort size, update calibration outcome parameters
# otherwise, need to choose cohort size by weighing speed increases of calculation against larger cohort size
if (df_sample_max$cs_option[1] == F) {
  l_params_calib$l_params_outcome <- l_params_outcome_cs
}


#### 5. Save parameters  ===========================================

saveRDS(l_params_calib, file = file_params_calib)


#### 6. Coverage analysis  ===========================================

# Run coverage analysis if necessary
if (check_coverage == T) {
  # Generate parameter sample with LHS
  m_param_samp <- lhs_param_samp(prior_map = l_params_calib$prior_map, 
                                 n_samp = n_samp_coverage)
  
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
  
  # Save parameter sample and corresponding outputs
  saveRDS(list(m_param_samp = m_param_samp, 
               m_calib_outputs = m_outputs,
               runtime = stime), 
          file = file_coverage)
  
  # Process targets
  df_targets <- l_params_calib$df_target %>%
    mutate(target_index = factor(target_index)) %>% # Create plot labels
    left_join(df_plot_labels, by = "target_groups")
  df_targets$plot_grps <- factor(df_targets$plot_grps, levels = df_plot_labels$plot_grps)
  
  # Plot and save coverage
  plt_coverage <- plot_coverage(df_targets = df_targets, 
                                m_outputs = m_outputs, 
                                file_fig_coverage = file_fig_coverage,
                                plt_size_text = plt_size_text)
  plt_coverage
}
