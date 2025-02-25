###########################  Internal Validation  #########################################
#
#  Objective: Validate BayCANN posteriors by plotting fit of calibration outputs
########################### <<<<<>>>>> ##############################################

rm(list = ls()) # Clean environment
options(scipen = 999) # View data without scientific notation

#### 1.Libraries and functions  ==================================================

###### 1.1 Load packages
library(tidyverse)
library(doBy)
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

# Get list of output file paths
l_filepaths_imabc <- update_config_paths("files_imabc", configs$paths)
l_filepaths_baycann <- update_config_paths("files_baycann", configs$paths)


#### 3. Pre-processing actions  ===========================================

# Load IMABC decision outputs if they exist
if (file.exists(l_filepaths_imabc$file_outputs)) {
  l_outputs_imabc <- readRDS(l_filepaths_imabc$file_outputs)
  l_outputs_imabc <- l_outputs_imabc$l_outputs
}

# Load BayCANN decision outputs if they exist
if (file.exists(l_filepaths_baycann$file_outputs)) {
  l_outputs_baycann <- readRDS(l_filepaths_baycann$file_outputs)
  
  # Extract calibration outputs and convert to data frame
  m_outputs <- do.call(rbind, lapply(l_outputs, function(u) {
    reshape_outputs(u[["outputs_base"]])
  }))
}


#### 4. Plots ===========================================
