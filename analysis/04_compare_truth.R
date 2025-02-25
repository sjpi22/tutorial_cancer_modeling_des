###########################  Ground Truth Comparison   ##########################
#
#  Objective: Compare IMABC and BayCANN posteriors to ground truth parameters
########################### <<<<<>>>>> #########################################

rm(list = ls()) # Clean environment
options(scipen = 999) # View data without scientific notation

#### 1.Libraries and functions  ==================================================

###### 1.1 Load packages
library(readxl)
library(tidyverse)

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
path_figs <- configs$paths$root_dir$figs
file_params_calib <- configs$paths$file_params_calib

# Get list of output file paths
l_filepaths_imabc <- update_config_paths("files_imabc", configs$paths)
l_filepaths_baycann <- update_config_paths("files_baycann", configs$paths)

###### 2.2 Other file paths
file_true_params <- file.path("_ground_truth", "true_params.xlsx")
file_fig_prior <- file.path(path_figs, "true_param_comparison.png")


#### 3. Pre-processing actions  ===========================================

# Load model and calibration parameters
l_params_calib <- readRDS(file_params_calib)

# Load true parameters and remove known parameters
df_true_params <- read_xlsx(file_true_params) %>%
  filter(var_id %in% l_params_calib$prior_map$var_id)

# Load IMABC posteriors
l_posteriors_imabc <- readRDS(l_filepaths_imabc$file_posterior)
m_params_imabc <- l_posteriors_imabc$good_parm_draws %>%
  dplyr::select(l_params_calib$prior_map$var_id)

# Load BayCANN posteriors
m_params_baycann <- read_csv(l_filepaths_baycann$file_posterior) %>%
  dplyr::select(-lp) # Remove last non-parameter column
names(m_params_baycann) <- l_params_calib$prior_map$var_id

# Combine IMABC and BayCANN posteriors in long form
df_params_long <- rbind(
  m_params_imabc %>%
    pivot_longer(everything()) %>%
    mutate(src = "IMABC"),
  m_params_baycann %>%
    pivot_longer(everything()) %>%
    mutate(src = "BayCANN"))


#### 4. Generate plots  ===========================================

# Plot BayCANN parameters against true parameters and priors
plot_params <- ggplot(df_params_long, aes(x = value, fill = src)) + 
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.7) + 
  geom_vline(data = df_true_params %>%
               rename(name = var_id), aes(xintercept = param_val), color = 'red') +
  geom_vline(data = l_params_calib$prior_map %>%
               rename(name = var_id), aes(xintercept = min), color = 'blue') +
  geom_vline(data = l_params_calib$prior_map %>%
               rename(name = var_id), aes(xintercept = max), color = 'blue') +
  facet_wrap(~name, scales = "free") +
  labs(x = "Parameter value", y = "Density", fill = "Method")
plot_params

# Save plot
ggsave(file_fig_prior, plot_params, width = 10, height = 8)
