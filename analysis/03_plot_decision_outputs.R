###########################  Internal Validation  #########################################
#
#  Objective: Validate BayCANN posteriors by plotting fit of calibration outputs
########################### <<<<<>>>>> ##############################################

rm(list = ls()) # Clean environment
options(scipen = 999) # View data without scientific notation

#### 1.Libraries and functions  ==================================================

###### 1.1 Load packages
library(tidyverse)
library(data.table)
library(patchwork)

###### 1.2 Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)


#### 2. General parameters ========================================================

###### 2.1 Configurations
# Run file to process configurations
source("configs/process_configs.R")

# Extract screening parameters
params_screen <- configs$params_screen

# Get list of output file paths
l_filepaths_imabc <- update_config_paths("files_imabc", configs$paths)
l_filepaths_baycann <- update_config_paths("files_baycann", configs$paths)

###### 2.2 Other parameters
v_methods <- c("imabc", "baycann")

#### 3. Pre-processing actions  ===========================================

# Load IMABC and BayCANN decision outputs if they exist
l_outputs_raw <- list() # Vector for holding saved raw outputs
l_outcomes <- list() # Vector for holding extracted outcomes
for (method in v_methods) {
  if (file.exists(get(paste0("l_filepaths_", method))$file_outputs)) {
    # Load data from file
    l_outputs_raw[[method]] <- readRDS(get(paste0("l_filepaths_", method))$file_outputs)
    l_outputs_raw[[method]] <- l_outputs_raw[[method]]$l_outputs
    
    # Loop over outcome types to extract data
    l_outcomes[[method]] <- list()
    for (outcome in c(names(params_screen$l_outcome_base), names(params_screen$l_outcome_counterfactual))) {
      l_outcomes[[method]][[outcome]] <- list()
      if (outcome %in% names(l_outputs_raw[[method]][[1]][["outputs_base"]])) {
        # Extract base scenario decision outputs to matrix
        m_outputs_base <- do.call(rbind, lapply(l_outputs_raw[[method]], function(u) {
          u[["outputs_base"]][[outcome]]
        }))
        
        # Save data
        l_outcomes[[method]][[outcome]][["base"]] <- m_outputs_base
      }
      
      if (outcome %in% names(l_outputs_raw[[method]][[1]][["outputs_screen"]][[1]])) {
        # Extract screening scenario decision outputs to matrix
        m_outputs_screen <- do.call(rbind, lapply(l_outputs_raw[[method]], function(u) {
          data.frame(scenario = names(u[["outputs_screen"]]),
                     do.call(rbind, lapply(names(u[["outputs_screen"]]), function(nm) {
                       u[["outputs_screen"]][[nm]][[outcome]]
                     })))
        }))
        
        # Save data
        l_outcomes[[method]][[outcome]][["screen"]] <- m_outputs_screen
      }
    }
    
    # Merge outcomes for plotting (LYG vs. test burden)
    l_outcomes[[method]][["plot_data"]] <- l_outcomes[[method]][["lyg"]][["screen"]] %>%
      # Merge base scenario N
      mutate(N = rep(l_outcomes[[method]][["lifeyears"]][["base"]][, "N"],
                     each = length(params_screen$strats))) %>%
      # Merge base scenario screening burden
      mutate(ct_tests_diag_C_base = rep(l_outcomes[[method]][["ntests"]][["base"]],
                                       each = length(params_screen$strats))) %>%
      # Merge screening test burden
      bind_cols(l_outcomes[[method]][["ntests"]][["screen"]] %>%
                  dplyr::select(-scenario))
  }
}

# Merge test and life years gained data

#### 4. Plots ===========================================
# Plot number of tests against life years gained across strategies

