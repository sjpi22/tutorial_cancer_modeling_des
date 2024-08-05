#* Title: Compare distribution density
#* 
#* Code function: This script compares the density of different distributions 
#*                that share the same mean value (10). This comparison is 
#*                made by calculating their median, 95% and 50% Interquantile
#*                range.
#* 
#* Creation date: July 11 2024
#* Author: David Garibay

# 01 Initial Setup --------------------------------------------------------

## 01.01 Clean environment ------------------------------------------------
remove(list = ls())

#* Refresh environment memory
gc()

## 01.02 Load libraries ----------------------------------------------------
library(dplyr)
library(ggplot2)
library(tidyr)





# 03 Define global variables --------------------------------------------

# Define vector of desired quantiles while manipulating the distributions
v_percentiles <- seq(0, 100, 0.1)

#* Define quantiles to define 95% and 50% Interquantile range (IQR) as well as 
#* the median
v_selected_quantiles <- c(0.025, 0.25, 0.5, 0.75, 0.975)

# 04 Compute distribution densities -------------------------------------

## 04.01 Exponential values ----------------------------------------------

# Exponential distribution parameters
exp_rate <- 0.1

#*  Calculate analytical mean and variance
exp_mean <- 1/exp_rate

# Obtain PDF
v_exp_pdf <- dexp(x = v_percentiles,
                  rate = exp_rate)

# Obtain CDF
v_exp_cdf <- pexp(q = v_percentiles,
                  rate = exp_rate)

v_exp_percentiles <- qexp(p = v_selected_quantiles,
                          rate = exp_rate)

# Create data frame for plotting
df_exp <- data.frame(quantile = v_percentiles,
                     pdf = v_exp_pdf,
                     cdf = v_exp_cdf)

# plt_exp_pdf <- ggplot(data = df_exp,
#                       mapping = aes(x = quantile, y = pdf)) + 
#   theme_bw(base_size = 18) +
#   geom_line() + 
#   geom_vline(xintercept = v_exp_percentiles[1], linetype = "dashed", colour = "dodgerblue2") +
#   geom_vline(xintercept = v_exp_percentiles[5], linetype = "dashed", colour = "dodgerblue2") + 
#   geom_vline(xintercept = v_exp_percentiles[2], linetype = "dashed", colour = "green4") +
#   geom_vline(xintercept = v_exp_percentiles[4], linetype = "dashed", colour = "green4") +
#   geom_vline(xintercept = v_exp_percentiles[3], linetype = "dashed", colour = "red")

## 04.02 Gamma values --------------------------------------------------------

# Gamma distribution paramaters
gamma_shape <- 10
gamma_scale <- 1

#*  Calculate analytical mean and variance
#*  Their definition was found in the documentation using `?rgamma()`
gamma_mean <- gamma_shape*gamma_scale
gamma_var  <- gamma_shape*(gamma_scale^2)

# Obtain PDF
v_gamma_pdf <- dgamma(x = v_percentiles,
                      shape = gamma_shape,
                      scale = gamma_scale)

# Obtain CDF
v_gamma_cdf <- pgamma(q = v_percentiles,
                      shape = gamma_shape,
                      scale = gamma_scale)

v_gamma_percentiles <- qgamma(p = v_selected_quantiles,
                              shape = gamma_shape,
                              scale = gamma_scale)

# Create data frame for plotting
df_gamma <- data.frame(quantile = v_percentiles,
                       pdf = v_gamma_pdf,
                       cdf = v_gamma_cdf)

# plt_gamma_pdf <- ggplot(data = df_gamma,
#                         mapping = aes(x = quantile, y = pdf)) + 
#   theme_bw(base_size = 18) +
#   geom_line() + 
#   geom_vline(xintercept = v_gamma_percentiles[1], linetype = "dashed", colour = "dodgerblue2") +
#   geom_vline(xintercept = v_gamma_percentiles[5], linetype = "dashed", colour = "dodgerblue2") + 
#   geom_vline(xintercept = v_gamma_percentiles[2], linetype = "dashed", colour = "green4") +
#   geom_vline(xintercept = v_gamma_percentiles[4], linetype = "dashed", colour = "green4") +
#   geom_vline(xintercept = v_gamma_percentiles[3], linetype = "dashed", colour = "red")


# 05 Plotting -----------------------------------------------------------

df_plt <- bind_rows(exp   = df_exp,
                    gamma = df_gamma,
                    .id = "distribution")


df_percentiles <- data.frame(
  quantile  = v_selected_quantiles,
  IQR       = c("95% IQR", "50% IQR", "Median", "50% IQR", "95% IQR"),
  exp       = v_exp_percentiles,
  gamma     = v_gamma_percentiles)


df_perc_long <- pivot_longer(data = df_percentiles,
                             cols = c("exp", "gamma"),
                             names_to = "distribution")

plt_distributions_pdf <- ggplot(data = df_plt,
                            mapping = aes(x = quantile, y = pdf)) +
  theme_bw(base_size = 21) +
  theme(legend.position = "bottom") +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  geom_vline(data = df_perc_long, 
             mapping = aes(xintercept = value, color = IQR),
             linetype = "dashed",
             linewidth = 0.8) +
  facet_wrap(~ distribution, nrow = 2) + 
  coord_cartesian(xlim = c(0, 50)) +
  labs(color = "Interquantile range:")
  

plt_distributions_pdf


ggsave(filename = "/Users/cide/Downloads/plt_distributions_pdf.png",
       plot = plt_distributions_pdf,
       width = 10,
       height = 8)
