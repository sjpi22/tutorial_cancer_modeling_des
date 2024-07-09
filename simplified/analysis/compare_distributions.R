# Program to compare differences in epidemiological outputs with different distributions

################################################################################
# Setup
################################################################################
# Clear workspace
rm(list = ls())

# Options
options(scipen=999,
        digits=2)

# Load packages
library(readxl)
library(data.table)
library(tidyverse)
library(patchwork)

# Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)

# Common parameters
v_ages_prevalence <- c(seq(30, 80, 10), 100) # Age for lesion prevalence
v_ages <- list(prevalence = v_ages_prevalence,
               incidence = v_ages_prevalence) # Age for cancer incidence 

################################################################################
# Version 1: After onset, all exponential distributions
################################################################################

#### Load and update parameters ####

# Load default data
l_params_exp <- load_default_params()

# Map variables to parameters for tuning - make dataframe of all parameters with "src = unknown"
param_map_exp <- make_param_map(l_params_exp)

# Update
# Mean equals 1/rate
# Variance equals 1/rate
v_update_exp <- c(1/5, 1/3, 1/6, 2, 75)
l_params_exp <- update_param_from_map(l_params_exp, v_update_exp, param_map_exp)
l_params_exp$n_cohort <- 500000

#### Initialize population and disease natural history ####
results_exp <- run_model(l_params_exp)
results_noscreening_exp <- results_exp[['None']]

# Get prevalence, incidence, and stage outputs
l_outputs_exp <- calc_calib_targets(l_params_exp, 
                                    results_noscreening_exp, 
                                    v_ages)

################################################################################
# Version 2: After onset, all gamma distributions
################################################################################

# Change distribution - equivalent to exponential if shape = 1, scale = 1/v_update_exp[i]
# Mean equals shape * scale
# Variance equals shape * scale^2
l_params_gamma <- copy(l_params_exp)

gamma_multiplier <- 2

l_params_gamma$time_1ii_2$distr <- "gamma"
l_params_gamma$time_1ii_2$params <- list(shape=1/v_update_exp[1] * gamma_multiplier, scale=1/gamma_multiplier)

l_params_gamma$time_1i_2$distr <- "gamma"
l_params_gamma$time_1i_2$params <- list(shape=1/v_update_exp[2] * gamma_multiplier, scale=1/gamma_multiplier)

l_params_gamma$time_1i_1ii$distr <- "gamma"
l_params_gamma$time_1i_1ii$params <- list(shape=1/v_update_exp[3] * gamma_multiplier, scale=1/gamma_multiplier)

# Map variables to parameters for tuning - make dataframe of all parameters with "src = unknown"
param_map_gamma <- make_param_map(l_params_gamma)

#### Initialize population and disease natural history ####
results_gamma <- run_model(l_params_gamma)
results_noscreening_gamma <- results_gamma[['None']]

# Get prevalence, incidence, and stage outputs
l_outputs_gamma <- calc_calib_targets(l_params_gamma, 
                                results_noscreening_gamma, 
                                v_ages)


################################################################################
# Version 3: After onset, all weibull distributions
################################################################################

# Change distribution - equivalent to exponential if shape = 1
# mean = scale * gamma(1 + 1/shape)
# variance = scale^2 * (gamma(1+2/shape) + gamma(1+1/shape)^2)
l_params_weibull <- copy(l_params_gamma)

weibull_shape <- 3

l_params_weibull$time_1ii_2$distr <- "weibull"
l_params_weibull$time_1ii_2$params <- list(shape=weibull_shape, scale=1/v_update_exp[1]/gamma(1 + 1/weibull_shape))

l_params_weibull$time_1i_2$distr <- "weibull"
l_params_weibull$time_1i_2$params <- list(shape=weibull_shape, scale=1/v_update_exp[2]/gamma(1 + 1/weibull_shape))

l_params_weibull$time_1i_1ii$distr <- "weibull"
l_params_weibull$time_1i_1ii$params <- list(shape=weibull_shape, scale=1/v_update_exp[3]/gamma(1 + 1/weibull_shape))

# Map variables to parameters for tuning - make dataframe of all parameters with "src = unknown"
param_map_weibull <- make_param_map(l_params_weibull)

#### Initialize population and disease natural history ####
results_weibull <- run_model(l_params_weibull)
results_noscreening_weibull <- results_weibull[['None']]

# Get prevalence, incidence, and stage outputs
l_outputs_weibull <- calc_calib_targets(l_params_weibull, 
                                      results_noscreening_weibull, 
                                      v_ages)

################################################################################
# Plot differences
################################################################################
l_targets <- list(Exponential = l_outputs_exp,
                  Gamma = l_outputs_gamma,
                  Weibull = l_outputs_weibull)
plot_targets <- plot_calib_targets(l_targets)

(plot_targets$prevalence + plot_targets$incidence + plot_targets$stage_distr) + 
  plot_layout(
    ncol = 3,
    guides = "collect"
  ) + plot_annotation(title = 'Calibration targets') &
  theme(legend.position = 'bottom')

# Get means and sample standard deviations of all time distributions
mean_vals <- list()
sd_vals <- list()
for (distr in c("exp", "gamma", "weibull")) {
  time_vars <- get(paste0("results_noscreening_", distr)) %>% 
    dplyr::select(starts_with("time_"))
  
  mean_vals[[distr]] <- time_vars %>%
    colMeans(na.rm=TRUE)
  
  sd_vals[[distr]] <- apply(time_vars, 2, sd, na.rm = TRUE) / sqrt(colSums(!is.na(time_vars)))
}

# Combine means as columns of dataframe
distr_means <- t(as.data.frame(do.call(rbind, mean_vals)))

# Get greatest difference in means
distr_means <- cbind(distr_means, 
                     mean_range = apply(distr_means, 1, function(x) diff(range(x))))

# Combine SDs as columns of dataframe
distr_sds <- t(as.data.frame(do.call(rbind, sd_vals)))

# Get max SD
distr_means <- cbind(distr_means,
                     max_sd = apply(distr_sds, 1, max))

# Calculate ratio of difference in means and SD
distr_means <- as.data.frame(distr_means) %>%
  mutate(sd_ratio = mean_range / max_sd)
View(distr_means)

