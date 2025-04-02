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
library(ggdist)
library(dampack)

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
v_quantiles <- configs$params_coverage$v_quantiles

# Get list of output file paths
l_filepaths_imabc <- update_config_paths("files_imabc", configs$paths)
l_filepaths_baycann <- update_config_paths("files_baycann", configs$paths)
l_filepaths_decision <- update_config_paths("files_decision", configs$paths)
list2env(l_filepaths_decision, envir = .GlobalEnv)

###### 2.2 Other parameters
v_methods <- c(truth = "Ground truth", imabc = "IMABC", baycann = "BayCANN") # Calibration methods to evaluate (include "truth" if evaluating against a ground truth)
base_test <- "confirm" # Assign name of base test (for diagnosing symptom-detected cases)
x_var <- "time_total" # Variable for x-axis
y_var <- "test_burden" # Variable for y-axis
x_int <- 200 # Interval for x-axis
y_int <- 2000 # Interval for y-axis
l_filepaths_truth <- list(file_outputs = file.path("_ground_truth", "true_decision_outputs.rds")) # Ground truth outputs


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

# Use uniform weight for BayCANN and truth
l_wts[["baycann"]] <- l_wts[["truth"]] <- 1

# Calculate quantiles and column labels from inner quantile vector
v_quantiles_lb <- (1 - v_quantiles[2]/100)/2
v_quantiles_ub <- (1 + v_quantiles[2]/100)/2
v_quantiles_calc <- sort(c(v_quantiles_lb, v_quantiles_ub))

# Load decision outputs if they exist
l_outputs_raw <- list() # Vector for holding saved raw outputs
l_outcomes <- list() # Vector for holding extracted outcomes
for (method in names(v_methods)) {
  if (file.exists(get(paste0("l_filepaths_", method))$file_outputs)) {
    # Load data from file
    l_outputs_raw[[method]] <- readRDS(get(paste0("l_filepaths_", method))$file_outputs)
    if ("l_outputs" %in% names(l_outputs_raw[[method]])) {
      l_outputs_raw[[method]] <- l_outputs_raw[[method]]$l_outputs
    }
    
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


#### 4. Plots and summary ===========================================

###### 4.1 Decision outcomes
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

# Restructure to plot ground truth average over other methods
# if ("Ground truth" %in% df_plot$method) {
#   # Remove ground truth from plot of all simulations
#   df_plot <- df_plot %>%
#     filter(method != "Ground truth")
#   
#   # Extract ground truth means
#   df_plot_mean_truth <- df_plot_mean %>%
#     filter(method == "Ground truth")
#   
#   # Remove ground truth from plot of all simulations and merge ground truth values
#   df_plot_mean <- df_plot_mean %>%
#     filter(method != "Ground truth") %>%
#     merge(df_plot_mean_truth[, c("scenario", x_var, y_var)], by = "scenario", suffixes = c("", "_truth"))
# }

# Calculate cost-efficiency frontier for plotting
df_plot_icer <- data.frame()
for (method in unique(df_plot$method)) {
  # Filter data to method
  df_plot_mean_method <- df_plot_mean[df_plot_mean$method == method, ]

  # Calculate ICERs for each method
  df_plot_icer_method <- calculate_icers(
    df_plot_mean_method$test_burden,
    df_plot_mean_method$time_total,
    df_plot_mean_method$scenario) %>%
    mutate(method = method)

  # Combine to full data frame
  df_plot_icer <- rbind(df_plot_icer, df_plot_icer_method)
}

# Rename variables in ICER data frame
colnames(df_plot_icer)[c(1:3, ncol(df_plot_icer))] <- c(colnames(df_plot_mean)[2], y_var, x_var, colnames(df_plot_mean)[1])

# Merge test interval and modality
df_plot_icer <- merge(df_plot_icer, df_intervals, by = "scenario")

# Plot number of tests against life years gained across strategies
plt_outcomes <- ggplot(df_plot, 
                       aes(x = get(x_var), y = get(y_var))) +
  geom_point(aes(color = factor(int_test)), # Plot cloud of points simulated from posterior
             alpha = 0.2, size = 1) +
  geom_line(data = df_plot_icer %>% # Plot cost-efficiency frontier
              filter(Status == "ND"), # Keep only non-dominated strategies
            color = "black",
            linewidth = 1) +
  geom_point(data = df_plot_mean, # Plot mean of each strategy
             aes(shape = factor(modality),
                 fill = factor(int_test)),
             color = "black",
             size = 3,
             stroke = 1) +
  # geom_point(data = df_plot_mean, # Plot true value for each strategy
  #            aes(x = get(paste0(x_var, "_truth")),
  #                y = get(paste0(y_var, "_truth"))),
  #            color = "black",
  #            alpha = 1,
  #            shape = 8,
  #            size = 3,
  #            stroke = 1) +
  scale_shape_manual(values = c(21, 24),
                     name = "Screening \nmodality", 
                     labels = c("Gold standard", "Non-invasive")) +
  facet_grid(~method) +
  labs(x = paste0("LYG per ", scales::label_comma()(unit)), 
       y = paste0("Additional confirmatory test burden per ", scales::label_comma()(unit)),
       color = "Screening \ninterval (years)") +
  scale_fill_discrete(guide = "none") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)),
         linetype = guide_legend(nrow = 2),
         shape = guide_legend(nrow = 2)) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma,
                     breaks = number_ticks(5)) +
  scale_linetype_discrete(name = "Screening \nmodality", 
                          labels = c("Gold standard", "Non-invasive")) +
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

###### 4.2 Dwell and sojourn time
l_res_time <- list()
for (method in names(v_methods)) {
  for (outcome in c("dwell_time", "sojourn_time")) {
    # Extract outcome
    v_outcomes <- l_outcomes[[method]][[outcome]][["base"]]
    
    # Calculate mean and quantiles of outputs
    if (length(l_wts[[method]]) == 1) {
      mean_outcome <- mean(v_outcomes)
      ci_outcome <- quantile(v_outcomes, probs = v_quantiles_calc)
    } else {
      mean_outcome <- weighted.mean(v_outcomes, w = l_wts[[method]])
      ci_outcome <- weighted_quantile(
        x = v_outcomes,
        probs = v_quantiles_calc,
        weights = l_wts[[method]]
      )
    }
    
    # Add to results list
    l_res_time <- c(l_res_time, list(
      list(method = method,
           outcome = outcome,
           mean = mean_outcome,
           ci_lb = ci_outcome[1],
           ci_ub = ci_outcome[2])))
  }
}

# Convert list to data frame and save
df_res_time <- rbindlist(l_res_time)
write.csv(df_res_time, file_nathist, row.names = FALSE)
