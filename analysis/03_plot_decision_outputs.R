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

# Extract relevant parameters from configs
params_screen <- configs$params_screen
unit <- params_screen$l_outcome_counterfactual$lyg$lit_params$unit # Unit for plotting outcomes
plt_size_text <- configs$params_coverage$plt_size_text # Plot text size

# Get list of output file paths
l_filepaths_imabc <- update_config_paths("files_imabc", configs$paths)
l_filepaths_baycann <- update_config_paths("files_baycann", configs$paths)
l_filepaths_decision <- update_config_paths("files_decision", configs$paths)
list2env(l_filepaths_decision, envir = .GlobalEnv)

###### 2.2 Other parameters
v_methods <- c(imabc = "IMABC", baycann = "BayCANN") # Calibration methods
v_screen_multiplier <- 2^seq(0, 3) # Number of screening tests equivalent to cost of one diagnostic test
x_var <- "time_total" # Variable for x-axis
y_var <- "burden_total" # Variable for y-axis
x_int <- 200 # Interval for x-axis
y_int <- 2000 # Interval for y-axis


#### 3. Pre-processing actions  ===========================================

# Get screening interval associated with each strategy
df_intervals <- data.frame(
  scenario = names(params_screen$strats),
  int_test = unname(sapply(names(params_screen$strats), function(u) {
    params_screen$strats[[u]]$int_screen
  }))
)

# Load IMABC and BayCANN decision outputs if they exist
l_outputs_raw <- list() # Vector for holding saved raw outputs
l_outcomes <- list() # Vector for holding extracted outcomes
for (method in names(v_methods)) {
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
                  dplyr::select(ct_tests_screen, ct_tests_diag_total)) %>%
      # Normalize test count by population (LYG already normalized)
      mutate(across(starts_with("ct_"), ~ . / N * unit)) %>%
      mutate(ct_tests_diag_diff = ct_tests_diag_total - ct_tests_diag_C_base)
    
    # Calculate test burden with appending cost fraction
    l_outcomes[[method]][["plot_data"]] <- crossing(l_outcomes[[method]][["plot_data"]],
                                                    multiplier = v_screen_multiplier) %>%
      mutate(burden_total = ct_tests_screen / multiplier + ct_tests_diag_diff) %>%
      # Merge test interval
      merge(df_intervals, by = "scenario")
  }
}


#### 4. Plots ===========================================

# Combine plot dataframes and label methods
df_plot <- data.frame()
for (method in names(l_outcomes)) {
  df_plot <- rbind(
    df_plot,
    l_outcomes[[method]][["plot_data"]] %>%
      mutate(method = v_methods[method])
  )
}
df_plot$method <- factor(df_plot$method, levels = v_methods)

# Get means for plotting lines
df_plot_mean <- df_plot %>%
  group_by(method, scenario, int_test, multiplier) %>%
  summarise(time_total = mean(time_total),
            burden_total = mean(burden_total), 
            .groups = "drop")

# Plot number of tests against life years gained across strategies
plt_outcomes <- ggplot(df_plot, 
                       aes(x = get(x_var), y = get(y_var))) +
  geom_line(data = df_plot_mean, 
            aes(linetype = factor(multiplier)),
            linewidth = 1, color = "gray") +
  geom_point(aes(color = factor(int_test)), 
             alpha = 0.3, size = 1) +
  facet_grid(~method) +
  labs(x = paste0("LYG per ", scales::label_comma()(unit)), 
       y = paste0("Additional diagnostic test burden per ", scales::label_comma()(unit)),
       color = "Screening \ninterval (years)",
       linetype = "Screening tests \nequivalent to cost of \none diagnostic test") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)),
         linetype = guide_legend(ncol = 2)) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma,
                     breaks = number_ticks(5)) +
  theme_bw(base_size = plt_size_text + 5) +
  theme(plot.title = element_text(size = plt_size_text, face = "bold"),
        axis.text.x = element_text(size = plt_size_text),
        axis.text.y = element_text(size = plt_size_text),
        axis.title = element_text(size = plt_size_text),
        legend.title = element_text(size = plt_size_text),
        legend.text = element_text(size = plt_size_text),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        legend.position = "bottom")
plt_outcomes

# Save plot
ggsave(file_fig_decision, plt_outcomes, width = 10, height = 8)
