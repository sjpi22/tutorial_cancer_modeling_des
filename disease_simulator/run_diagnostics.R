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

# Check if directory exists, make if not
dir.create(file.path(data_outpath), showWarnings = FALSE)

# Load default data
l_params_all <- load_default_params()

# Map variables to parameters for tuning - make dataframe of all parameters with "src = unknown"
param_map <- make_param_map(l_params_all)

# Merge true and default parameters
compare_param_map <- data.frame(compare_param_map,
                                default_val = param_map$param_val)

# Compare multiplied difference
within(compare_param_map,
       ratio_of_default = param_val / default_val)

View(compare_param_map)
