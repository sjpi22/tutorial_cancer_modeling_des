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
library(data.table)

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
plt_size_text <- configs$params_coverage$plt_size_text # Plot text size

# Get list of output file paths
l_filepaths_imabc <- update_config_paths("files_imabc", configs$paths)
l_filepaths_baycann <- update_config_paths("files_baycann", configs$paths)

###### 2.2 Other file paths
file_true_params <- file.path("_ground_truth", "true_params.xlsx")
file_fig_prior <- file.path(path_figs, "true_param_comparison.png")

###### 2.3 Other parameters
v_methods <- c(imabc = "IMABC", baycann = "BayCANN") # Calibration methods


#### 3. Pre-processing actions  ===========================================

# Load model and calibration parameters
l_params_calib <- readRDS(file_params_calib)

# Load true parameters
df_true_params <- read_xlsx(file_true_params) 

# Update true parameters table if using hazard ratios to parametrize cancer progression variables
if (!is.null(l_params_calib$l_params_model$hr_cancer)) {
  # Update stage distribution parameters
  df_true_params_stage <- df_true_params %>%
    # Subset to stage distribution variables
    filter(var_name == "p_cancer") %>%
    # Replace name of p_cancer variable name and ID
    # and calculate true difference between true and observed stage distribution
    mutate(var_name = "diff_p_cancer",
           param_val = param_val - l_params_calib$l_params_model$p_cancer) %>%
    # Rename variable id
    mutate(var_id = paste(var_name, idx, sep = "."))
  
  # Bind to true parameter table
  df_true_params <- bind_rows(df_true_params, df_true_params_stage)
}

# Remove known parameters
df_true_params <- df_true_params %>%
  filter(var_id %in% l_params_calib$prior_map$var_id)

# Load posteriors
l_m_params <- list()
if ("imabc" %in% tolower(v_methods)) {
  l_posteriors_imabc <- readRDS(l_filepaths_imabc$file_posterior)
  l_m_params[[v_methods["imabc"]]] <- l_posteriors_imabc$good_parm_draws %>%
    dplyr::select(l_params_calib$prior_map$var_id)
}

# Load BayCANN posteriors
if ("baycann" %in% tolower(v_methods)) {
  m_params_baycann <- read_csv(l_filepaths_baycann$file_posterior) %>%
    dplyr::select(-lp) # Remove last non-parameter column
  names(m_params_baycann) <- l_params_calib$prior_map$var_id
  l_m_params[[v_methods["baycann"]]] <- m_params_baycann
}

# Combine IMABC and BayCANN posteriors in long form
df_params_long <- rbindlist(l_m_params, idcol = "method") %>%
  pivot_longer(l_params_calib$prior_map$var_id)
df_params_long$method <- factor(df_params_long$method, levels = v_methods)


#### 4. Generate plots  ===========================================

# Plot BayCANN parameters against true parameters and priors
plt_params <- ggplot(df_params_long, aes(x = value, fill = method)) + 
  geom_histogram(aes(y = after_stat(density)), position = "identity", alpha = 0.7) + 
  geom_vline(data = df_true_params %>%
               rename(name = var_id), aes(xintercept = param_val), color = 'red') +
  geom_vline(data = l_params_calib$prior_map %>%
               rename(name = var_id), aes(xintercept = min), color = 'blue') +
  geom_vline(data = l_params_calib$prior_map %>%
               rename(name = var_id), aes(xintercept = max), color = 'blue') +
  facet_wrap(~name, scales = "free") +
  labs(x = "Parameter value", y = "Density", fill = "Method") +
  scale_x_continuous(breaks = number_ticks(3)) +
  theme_bw(base_size = plt_size_text + 5) +
  theme(plot.title = element_text(size = plt_size_text, face = "bold"),
        axis.text.x = element_text(size = plt_size_text),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = plt_size_text),
        legend.title = element_text(size = plt_size_text),
        legend.text = element_text(size = plt_size_text),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        legend.position = "bottom")
plt_params

# Save plot
ggsave(file_fig_prior, plt_params, width = 12, height = 10)
