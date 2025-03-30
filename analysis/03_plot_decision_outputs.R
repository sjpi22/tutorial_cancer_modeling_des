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
test_cost <- sapply(params_screen$test_chars, function(u) { # Number of tests equivalent to cost of one confirmatory test
  u[["cost_ratio"]]
})
plt_size_text <- configs$params_coverage$plt_size_text # Plot text size

# Get list of output file paths
l_filepaths_imabc <- update_config_paths("files_imabc", configs$paths)
l_filepaths_baycann <- update_config_paths("files_baycann", configs$paths)
l_filepaths_decision <- update_config_paths("files_decision", configs$paths)
list2env(l_filepaths_decision, envir = .GlobalEnv)

###### 2.2 Other parameters
v_methods <- c(imabc = "IMABC", baycann = "BayCANN") # Calibration methods
base_test <- "confirm" # Assign name of base test (for diagnosing symptom-detected cases)
x_var <- "time_total" # Variable for x-axis
y_var <- "test_burden" # Variable for y-axis
x_int <- 200 # Interval for x-axis
y_int <- 2000 # Interval for y-axis


#### 3. Pre-processing actions  ===========================================

# Extend base costs
test_cost["base"] <- -test_cost[base_test] # For subtracting cost of base test in non-screening scenario
test_cost["base_cf"] <- test_cost[base_test] # For adding cost of base test in screening scenario

# Get screening interval associated with each strategy
df_intervals <- data.frame(
  scenario = names(params_screen$strats),
  modality = unname(sapply(names(params_screen$strats), function(u) {
    params_screen$strats[[u]]$mod
  })),
  modality_conf = unname(sapply(names(params_screen$strats), function(u) {
    ifelse(is.null(params_screen$strats[[u]]$mod_conf), NA, params_screen$strats[[u]]$mod_conf)
  })),
  int_test = unname(sapply(names(params_screen$strats), function(u) {
    params_screen$strats[[u]]$int_screen
  }))
)

# Load IMABC parameter weights
l_wts <- list()
l_posteriors_imabc <- readRDS(l_filepaths_imabc$file_posterior)
l_wts[["imabc"]] <- l_posteriors_imabc$good_parm_draws$sample_wt

# Use uniform weight for BayCANN
l_wts[["baycann"]] <- 1

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
      # Extract base scenario decision outputs to matrix
      if (outcome %in% names(l_outputs_raw[[method]][[1]][["outputs_base"]])) {
        m_outputs_base <- do.call(rbind, lapply(l_outputs_raw[[method]], function(u) {
          u[["outputs_base"]][[outcome]]
        }))
        
        # Save data
        l_outcomes[[method]][[outcome]][["base"]] <- m_outputs_base
        
        # Truncate length of weights (if unable to produce data for whole posterior)
        if (length(l_wts[[method]]) > 1) {
          l_wts[[method]] <- l_wts[[method]][1:nrow(m_outputs_base)]
        }
      }
      
      # Extract screening scenario decision outputs to matrix
      if (outcome %in% names(l_outputs_raw[[method]][[1]][["outputs_screen"]][[1]])) {
        m_outputs_screen <- do.call(rbind, lapply(l_outputs_raw[[method]], function(u) {
          data.frame(scenario = names(u[["outputs_screen"]]),
                     do.call(rbindlist, list(l = lapply(names(u[["outputs_screen"]]), function(nm) {
                       data.frame(t(u[["outputs_screen"]][[nm]][[outcome]]))
                     }),
                     use.names = TRUE,
                     fill = TRUE)))
        }))
        
        # Save data
        l_outcomes[[method]][[outcome]][["screen"]] <- m_outputs_screen
      }
    }
    
    # Merge outcomes for plotting (LYG vs. test burden)
    l_outcomes[[method]][["plot_data"]] <- l_outcomes[[method]][["lyg"]][["screen"]] %>%
      # Append weights
      mutate(wt = rep(l_wts[[method]], each = n()/length(l_wts[[method]]))) %>%
      # Merge base scenario N
      mutate(N = rep(l_outcomes[[method]][["lifeyears"]][["base"]][, "N"],
                     each = length(params_screen$strats))) %>%
      # Merge base scenario screening burden
      mutate(ct_tests_base = rep(l_outcomes[[method]][["ntests"]][["base"]],
                                        each = length(params_screen$strats))) %>%
      # Merge screening test burden
      bind_cols(l_outcomes[[method]][["ntests"]][["screen"]] %>%
                  dplyr::select(starts_with("ct_tests_")) %>%
                  rename(c(ct_tests_base_cf = "ct_tests_base"))) %>%
      # Normalize test count by population (LYG already normalized)
      mutate(across(starts_with("ct_"), ~ . / N * unit)) %>%
      # Merge test interval
      merge(df_intervals, by = "scenario") %>%
      setDT()
    
    # Calculate final test cost variable
    l_outcomes[[method]][["plot_data"]][, test_burden := rowSums(mapply(`*`, .SD, test_cost), na.rm = TRUE), .SDcols = paste0("ct_tests_", names(test_cost))]
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
  group_by(method, scenario, modality, int_test) %>%
  summarise(time_total = weighted.mean(time_total, wt),
            test_burden = weighted.mean(test_burden, wt), 
            .groups = "drop")

# Plot number of tests against life years gained across strategies
plt_outcomes <- ggplot(df_plot, 
                       aes(x = get(x_var), y = get(y_var))) +
  geom_line(data = df_plot_mean, 
            aes(linetype = factor(modality)),
            linewidth = 1, color = "gray") +
  geom_point(aes(color = factor(int_test)), 
             alpha = 0.2, size = 1) +
  geom_point(data = df_plot_mean,
             aes(shape = factor(modality),
                 fill = factor(int_test)),
             color = "black",
             size = 2) +
  scale_shape_manual(values=21:22) +
  facet_grid(~method) +
  labs(x = paste0("LYG per ", scales::label_comma()(unit)), 
       y = paste0("Additional confirmatory test burden per ", scales::label_comma()(unit)),
       color = "Screening \ninterval (years)") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)),
         linetype = guide_legend(nrow = 2)) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma,
                     breaks = number_ticks(5)) +
  scale_linetype_discrete(name = "Screening modality", labels = c("Confirmatory", "Non-invasive")) +
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
