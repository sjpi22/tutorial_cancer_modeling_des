###########################  Generate BayCANN Sample  ##########################
#
#  Objective: Program to simulate parameter inputs and model outputs for 
#  BayCANN model calibration
########################### <<<<<>>>>> #########################################


#### 1.Libraries and functions  ==================================================
# Clean environment
rm(list = ls())

library(data.table)
library(tidyverse)
library(lhs)
library(doParallel)
library(foreach)
library(assertthat)

###### 1.1 Load functions =================================================
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)


#### 2. General parameters ========================================================

###### 2.1 file paths
file_params <- "data/calibration_params.rds"
outpath <- "output/calibration/BayCANN"
file_sample <- file.path(outpath, "sample_BayCANN.rds")
file_fig_coverage <- file.path(outpath, "plots", "fig_coverage_BayCANN.png")

###### 2.2 BayCANN parameters
n_samp_per_param <- 1000


#### 3. Pre-processing actions  ===========================================

# Load model and calibration parameters
l_params_calib <- readRDS(file_params)

# Set number of cores to use
registerDoParallel(cores=detectCores(logical = TRUE) - l_params_calib$n_cores_reserved_local)


#### 4. Generate random set of inputs  ===========================================

m_param_samp <- with(l_params_calib, {
  # Get number of params to calibrate and number of samples
  n_param <- nrow(prior_map)
  n_samp <- n_samp_per_param * n_param
    
  # Sample unit Latin Hypercube
  m_lhs_unit <- randomLHS(n_samp, n_param)
  
  # Rescale to min/max of each parameter
  m_param_samp <- matrix(nrow = n_samp, ncol = n_param)
  for (i in 1:n_param) {
    m_param_samp[, i] <- qunif(
      m_lhs_unit[, i],
      min = prior_map$min[i],
      max = prior_map$max[i])
  }
  colnames(m_param_samp) <- prior_map$var_id
  
  return(m_param_samp)
}
)


#### 5. Generate corresponding outputs  ===========================================

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
        params_to_calib_outputs(
          l_params_all = l_params_all,
          v_params_update = v_params_update,
          param_map = prior_map,
          l_outcome_params = l_outcome_params,
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
colnames(m_outputs) <- l_params_calib$df_true_targets$target_names

# Check for any NaN
validate_that(
  sum(sapply(m_param_samp, function(x) any(is.nan(x)))) == 0, 
  msg = "Parameters include NaN")
validate_that(
  sum(sapply(m_outputs, function(x) any(is.nan(x)))) == 0, 
  msg = "Outputs include NaN")

# Save parameter sample and corresponding outputs
saveRDS(list(m_param_samp = m_param_samp, 
             m_calib_outputs = m_outputs), 
        file = file_sample)


#### 6. Coverage analysis  ===========================================

# Process targets
df_targets <- l_params_calib$df_true_targets %>%
  mutate(target_index = factor(target_index)) %>% # Create plot labels
  mutate(plot_grp = case_when(target_groups == "stage_distr" ~ "Stage distribution",
                              target_groups == "incidence" ~ "Incidence per 100k by age",
                              target_groups == "prevalence" ~ "Prevalence by age"))

# Convert outputs from wide to long
out_full_bc_cat <- data.frame(m_outputs) %>%
  pivot_longer(
    cols = everything(), 
    names_to = "target_names",
    values_to = "value"
  ) %>%
  mutate(target_index = rep(df_targets$target_index, nrow(m_outputs)),
         plot_grp = rep(df_targets$plot_grp, nrow(m_outputs)))

# Plot distribution of outputs against targets
plot_coverage <- ggplot(data = df_targets) + 
  geom_errorbar(
    aes(x    = target_index, 
        y    = targets, 
        ymin = targets - se, 
        ymax = targets + se),
    width = 0.4, linewidth = 0.9, color="red") +
  theme(legend.position="none") +
  geom_violin(data = out_full_bc_cat,
              aes(x    = target_index,
                  y    = value),
              alpha = 0.4) +
  facet_wrap(~ plot_grp, scales="free") +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(), legend.position="none") +
  scale_fill_manual(values = c("grey10", "grey30"))+
  scale_y_continuous(breaks = number_ticks(5))+
  theme_bw(base_size = 23) +
  theme(plot.title = element_text(size = 22, face = "bold"),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0)) +
  labs(x     = "", y     = "")
plot_coverage

ggsave(file_fig_coverage, plot = plot_coverage,
       width = length(unique(df_targets$plot_grp))*4, height = 6)
