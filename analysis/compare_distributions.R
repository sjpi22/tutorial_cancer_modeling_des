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
v_ages_prevalence <- seq(30, 80, 10) # Age for lesion prevalence
v_ages <- list(prevalence = v_ages_prevalence,
               incidence = v_ages_prevalence) # Age for cancer incidence 

################################################################################
# Version 1: After onset, all exponential distributions
# Mean equals 1/rate
# Variance equals 1/rate
################################################################################

#### Load and update parameters ####

# Load default data
l_params_exp <- load_default_params()

# Map variables to parameters for tuning - make dataframe of all parameters with "src = unknown"
param_map_exp <- make_param_map(l_params_exp)

# Update
v_update_exp <- c(1/3, 1/5, 1/4, 2, 75)
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
# Version 2: After onset, all weibull distributions
# equivalent to exponential if shape = 1
# mean = scale * gamma(1 + 1/shape)
# variance = scale^2 * (gamma(1+2/shape) + gamma(1+1/shape)^2)
################################################################################

# Weibull
exp_means <- 1/v_update_exp[1:3]
weibull_shape <- c(3, 3, 3)
weibull_scale <- exp_means / gamma(1 + 1/weibull_shape)

# Change distribution
l_params_weibull <- copy(l_params_exp)

l_params_weibull$time_1ii_2$distr <- "weibull"
l_params_weibull$time_1ii_2$params <- list(shape=weibull_shape[1], scale=weibull_scale[1])

l_params_weibull$time_1i_2$distr <- "weibull"
l_params_weibull$time_1i_2$params <- list(shape=weibull_shape[2], scale=weibull_scale[2])

l_params_weibull$time_1i_1ii$distr <- "weibull"
l_params_weibull$time_1i_1ii$params <- list(shape=weibull_shape[3], scale=weibull_scale[3])

#### Initialize population and disease natural history ####
results_weibull <- run_model(l_params_weibull)
results_noscreening_weibull <- results_weibull[['None']]

# Get prevalence, incidence, and stage outputs
l_outputs_weibull <- calc_calib_targets(l_params_weibull, 
                                        results_noscreening_weibull, 
                                        v_ages)

################################################################################
# Version 3: After onset, all gamma distributions
# Mean equals shape * scale
# Variance equals shape * scale^2
################################################################################

# Change distribution - equivalent to exponential if shape = 1, scale = 1/v_update_exp[i]

l_params_gamma <- copy(l_params_exp)

gamma_shape <- gamma(1 + 1/weibull_shape)^2 / (gamma(1 + 2/weibull_shape) - gamma(1 + 1/weibull_shape)^2)
gamma_scale <- exp_means / gamma_shape

l_params_gamma$time_1ii_2$distr <- "gamma"
l_params_gamma$time_1ii_2$params <- list(shape=gamma_shape[1], scale=gamma_scale[1])

l_params_gamma$time_1i_2$distr <- "gamma"
l_params_gamma$time_1i_2$params <- list(shape=gamma_shape[2], scale=gamma_scale[2])

l_params_gamma$time_1i_1ii$distr <- "gamma"
l_params_gamma$time_1i_1ii$params <- list(shape=gamma_shape[3], scale=gamma_scale[3])


#### Initialize population and disease natural history ####
results_gamma <- run_model(l_params_gamma)
results_noscreening_gamma <- results_gamma[['None']]

# Get prevalence, incidence, and stage outputs
l_outputs_gamma <- calc_calib_targets(l_params_gamma, 
                                results_noscreening_gamma, 
                                v_ages)


################################################################################
# Version 4: After onset, all binary distributions
# mean = probs[1] * x[1] + probs[2] * x[2]
# var = probs[1] * (x[1] - mean)^2 + probs[2] * (x[2] - mean)^2
################################################################################

# Change distribution
l_params_bin <- copy(l_params_exp)

bin_diff <- sqrt(gamma_shape * gamma_scale^2)
bin_min <- exp_means - bin_diff
bin_max <- exp_means + bin_diff

l_params_bin$time_1ii_2$distr <- "empirical"
l_params_bin$time_1ii_2$params <- list(xs = c(bin_min[1], bin_max[1]), probs = c(0.5, 0.5), continuity_correction = NULL)

l_params_bin$time_1i_2$distr <- "empirical"
l_params_bin$time_1i_2$params <- list(xs = c(bin_min[2], bin_max[2]), probs = c(0.5, 0.5), continuity_correction = NULL)

l_params_bin$time_1i_1ii$distr <- "empirical"
l_params_bin$time_1i_1ii$params <- list(xs = c(bin_min[3], bin_max[3]), probs = c(0.5, 0.5), continuity_correction = NULL)

#### Initialize population and disease natural history ####
results_bin <- run_model(l_params_bin)
results_noscreening_bin <- results_bin[['None']]

# Get prevalence, incidence, and stage outputs
l_outputs_bin <- calc_calib_targets(l_params_bin, 
                                        results_noscreening_bin, 
                                        v_ages)


################################################################################
# Version 5: After onset, all uniform distributions with equal variance as weibull/gamma
# mean: 1/2 * (min + max)
# variance: 1/12 * (max - min)^2
################################################################################

# Change distribution
l_params_unif <- copy(l_params_exp)

# Derive parameters
# min + max = 2 * exp_means
# max - min = sqrt(12 * gamma_shape * gamma_scale^2)
# 2 * max = 2 * exp_means + sqrt(12 * gamma_shape * gamma_scale^2)
unif_max = (2 * exp_means + sqrt(12 * gamma_shape * gamma_scale^2)) / 2
unif_min = 2 * exp_means - unif_max
  
l_params_unif$time_1ii_2$distr <- "unif"
l_params_unif$time_1ii_2$params <- list(min = unif_min[1], max = unif_max[1])

l_params_unif$time_1i_2$distr <- "unif"
l_params_unif$time_1i_2$params <- list(min = unif_min[2], max = unif_max[2])

l_params_unif$time_1i_1ii$distr <- "unif"
l_params_unif$time_1i_1ii$params <- list(min = unif_min[3], max = unif_max[3])

#### Initialize population and disease natural history ####
results_unif <- run_model(l_params_unif)
results_noscreening_unif <- results_unif[['None']]

# Get prevalence, incidence, and stage outputs
l_outputs_unif <- calc_calib_targets(l_params_unif, 
                                    results_noscreening_unif, 
                                    v_ages)

# Note: uniform distributions with equal variance as exp do not work because range would need to include negatives





################################################################################
# Plot differences
################################################################################

# Plot distributions
plot_distr_df_wide <- data.frame()
x_diff <- 0.1
x_distr <- seq(0, 10, x_diff)
for (i in 1:length(exp_means)) {
  # Make data frame for continuous distributions
  temp_plot_distr_df <- data.frame(
    var = param_map_exp$var_name[i],
    time = x_distr,
    Exponential = dexp(x_distr, 1/exp_means[i]),
    Weibull = dweibull(x_distr, shape = weibull_shape[i], scale = weibull_scale[i]),
    Gamma = dgamma(x_distr, shape = gamma_shape[i], scale = gamma_scale[i]),
    Uniform = dunif(x_distr, min = unif_min[i], max = unif_max[i]),
    Binary = ifelse(abs(x_distr - bin_min[i]) < x_diff/2 | abs(x_distr - bin_max[i]) < x_diff/2, 0.5, 0)
  )
  
  plot_distr_df_wide <- rbind(plot_distr_df_wide, temp_plot_distr_df)
    
}

plot_distr_df <- plot_distr_df_wide %>%
  pivot_longer(cols = 3:ncol(plot_distr_df_wide))

plot_distr <- ggplot(plot_distr_df, aes(x = time, y = value, color = name)) +
  geom_line() +
  facet_wrap(~var, scales="free", ncol = 3) +
  labs(x = 'Time', y = 'Probability density', color = 'Label')

# Plot targets
l_targets <- list(Exponential = l_outputs_exp,
                  Gamma = l_outputs_gamma,
                  Weibull = l_outputs_weibull,
                  Binary = l_outputs_bin,
                  Uniform = l_outputs_unif)
plot_targets <- plot_calib_targets(l_targets)

((plot_targets$prevalence | plot_targets$incidence | plot_targets$stage_distr) / plot_distr) + 
  plot_layout(
    guides = "collect"
  ) + plot_annotation(title = 'Calibration targets') 

# Get means and sample standard deviations of all time distributions
mean_vals <- list()
sd_vals <- list()
for (distr in c("exp", "gamma", "weibull", "bin", "unif")) {
# for (distr in c("exp",  "bin")) {
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


# Plot means and variances
mean_weibull <- function(shape, scale) {
  return(shape * gamma(1 + 1/scale))
}

var_weibull <- function(shape, scale) {
  return(shape^2 * (gamma(1 + 2/scale) - (gamma(1 + 1/scale))^2))
}

shape_seq <- seq(10)
scale_seq <- seq(1, 5, 0.1)

weibull_df <- data.frame(
  shape = rep(shape_seq, length(scale_seq)),
  scale = rep(scale_seq, length(shape_seq)),
  quantity = "mean"
) %>%
  mutate(val = mean_weibull(shape, scale))

weibull_df <- rbind(weibull_df,
                    data.frame(
                      shape = rep(shape_seq, length(scale_seq)),
                      scale = rep(scale_seq, length(shape_seq)),
                      quantity = "sd"
                    ) %>%
                      mutate(val = sqrt(var_weibull(shape, scale)))
                    )

ggplot(weibull_df, aes(scale, val, color = factor(shape), linetype = quantity)) +
  geom_line()



######### Formulas to get the same mean and variance for weibull and gamma
# Weibull

lambda <- 3
k <- 5
weibull_mean <- lambda * gamma(1 + 1/k)
weibull_var <- lambda^2 * (gamma(1 + 2/k) - gamma(1 + 1/k)^2)
weibull_sample <- rweibull(5000, shape = k, scale = lambda)
print(paste("Weibull mean", round(weibull_mean, 2), "(check with MC simulation:", round(mean(weibull_sample), 2), ")"))
print(paste("Weibull variance", round(weibull_var, 2), "(check with MC simulation:", round(sd(weibull_sample)^2, 2), ")"))


# Gamma
alpha <- gamma(1 + 1/k)^2 / (gamma(1 + 2/k) - gamma(1 + 1/k)^2)
beta <- weibull_mean / alpha
gamma_mean <- alpha * beta
gamma_var <- alpha * beta^2
gamma_sample <- rgamma(5000, shape = alpha, scale = beta)
print(paste("gamma mean", round(gamma_mean, 2), "(check with MC simulation:", round(mean(gamma_sample), 2), ")"))
print(paste("gamma variance", round(gamma_var, 2), "(check with MC simulation:", round(sd(gamma_sample)^2, 2), ")"))

