# Run diagnostics to ensure that default parameters are in search range of true parameters

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

# Load true parameters
true_param_path <- 'disease_simulator/true_param_map_consistent.RData'
load(true_param_path)

# Load default data
l_params_all <- load_default_params()

# Map variables to parameters for tuning - make dataframe of all parameters with "src = unknown"
default_param_map <- make_param_map(l_params_all)

# Merge true and default parameters
compare_param_map <- data.frame(param_map,
                                default_val = default_param_map$param_val)

# Compare multiplied difference
within(compare_param_map, { ratio_of_default = param_val / default_val})

View(compare_param_map)

 
# Load generated sample mapping file and check if true val is within bounds
load(sample_file)
param_map <- param_map %>%
  mutate(within_bounds = (compare_param_map$param_val >=  prior_min) & (compare_param_map$param_val <= prior_max))

# Check generated sample data
sherlock_sample_data <- 'sherlock/calibration_sample.RData'
load(sherlock_sample_data)
nrow(out_calib_targets) 


# Check IMABC data
sherlock_imabc_params_file <- 'sherlock/output/imabc/SimulatedParameters_20240623_2209PDT.csv'
sherlock_imabc_targets_file <- 'sherlock/output/imabc/SimulatedTargets_20240623_2209PDT.csv'
imabc_params <- read.csv(sherlock_imabc_params_file)
imabc_targets <- read.csv(sherlock_imabc_targets_file)
