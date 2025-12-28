###########################  Unit test: Natural history model   ################
#
#  Objective: Run unit tests for natural history model
########################### <<<<<>>>>> ##############################################


#### 1.Libraries and functions  ==================================================
#* Clean environment
rm(list = ls())

library(tidyverse)
library(testthat)

###### 1.1 Load functions =================================================

# Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)


# Define functions to calculate theoretical parameters
f_mean <- function(x, distr_obj) {
  x * query_distr("d", x, distr_obj$distr, distr_obj$params)
}

f_var <- function(x, distr_obj, distr_mean) {
  (x - distr_mean)^2 * query_distr("d", x, distr_obj$distr, distr_obj$params)
}


#### 2. General parameters ========================================================

###### 2.1 File paths
file_model_params <- file.path("_ground_truth", "params_model.rds")

###### 2.2 Other parameters
eps_pct <- 0.01 # tolerated error
p_alpha <- 0.02 # p-value threshold


#### 3. Pre-processing actions  ===========================================

# Set seed
set.seed(2025)

# Load ground truth parameters
l_params_model <- readRDS(file_model_params)

# Simulate natural history model
res <- run_base_model(l_params_model)


#### 4. Tests  ===========================================

# Extract distributions from parameter list
v_distr <- names(l_params_model)
v_distr <- v_distr[grep("^d_", v_distr)]

# Test that time to event distributions are simulated correctly
for (varname in v_distr) {
  # Get variable distribution
  d_test <- l_params_model[[varname]]
  
  # Extract variable removing NAs
  if (varname == "d_time_L_P") { # Special case - lesion-level variable with different name
    v_samp <- res$lesion_level[["time_Lj_Pj"]]
  } else if (varname == "d_n_L") {
    v_samp <- res$patient_level[["n_L_add"]]
  } else if (grepl("^d_time_C._Dc$", varname)) {
    v_samp <- res$patient_level[stage_dx == as.integer(substr(varname, 9, 9))][["time_C_Dc"]]
  } else {
    v_samp <- res$patient_level[[substr(varname, 3, nchar(varname))]]
    if (varname %in% c("d_time_H_Do")) {
      v_male <- res$patient_level$male # Extract sex if needed
      v_male <- v_male[!is.na(v_samp)]
    } else if (varname == "d_time_L_Lj") {
      # Extract scaled time to lesion onset
      v_samp <- res$lesion_level[lesion_id > 1, .(pct = time_L_Lj / (time_H_Do - time_H_L))]
      v_samp <- v_samp$pct
    }
  }
  v_samp <- v_samp[!is.na(v_samp)]
  
  if (varname == "d_n_L") { # Special case - rate multiplied by time from lesion onset to death from other causes for number of additional lesions
    check <- res$patient_level[(time_H_Do > time_H_L) != (!is.na(n_L_add))]
    
    # Recalculate as rate * time from lesion onset to death from other causes
    expected_n <- unlist(res$patient_level[time_H_Do > time_H_L, .(expected_n = d_test$params$lambda * (time_H_Do - time_H_L))])
    
    # Check expected and observed number of lesions 
    test_that("Simulated sample matches distribution", {
      # Perform chi-squared test
      chisq_test <- chisq.test(v_samp, expected_n)
      expect_gt(chisq_test$p.value, p_alpha) # p-value > 0.05 indicates a good fit
    })
  } else {
    # Define theoretical parameters
    if (varname == "d_time_H_Do") {
      # Weights by sex
      v_weights <- c(male = l_params_model$d_male$params$prob,
                     female = 1 - l_params_model$d_male$params$prob)
      
      # Calculated weighted average time to death from other causes
      expected_mean <- c()
      for (sex in names(d_test)) {
        expected_mean[sex] <- integrate(function(x) f_mean(x, d_test[[sex]]), 0, Inf)$value
      }
      expected_mean <- sum(expected_mean*v_weights[names(expected_mean)])
      expected_sd <- NA # Ignore SD for this variable
    } else if (d_test$distr == "binom") {
      expected_mean <- d_test$params$size * d_test$params$prob
      expected_sd <- sqrt(d_test$params$size * d_test$params$prob * (1 - d_test$params$prob))
    } else {
      expected_mean <- integrate(function(x) f_mean(x, d_test), 0, Inf)$value
      expected_sd <- sqrt(integrate(function(x) f_var(x, d_test, expected_mean), 0, Inf)$value)
    }
    
    # Set tolerance as % of expected mean or sampling error: 
    tol_var <- max(qnorm(1 - p_alpha/2) * expected_sd / sqrt(length(v_samp)),
                   expected_mean * eps_pct,
                   na.rm = T)
    
    test_that("Simulated sample matches distribution", {
      # Check mean
      sample_mean <- mean(v_samp)
      expect_equal(sample_mean, expected_mean, tolerance = tol_var)
      
      # Check standard deviation
      if (varname != "d_time_H_Do") {
        sample_sd <- sd(v_samp)
        expect_equal(sample_sd, expected_sd, tolerance = tol_var)
        
        # Kolmogorov-Smirnov test for goodness-of-fit
        if (!d_test$distr %in% c("binom")) {
          ks_test <- do.call(ks.test, c(x = list(v_samp), y = list(paste0("p", d_test$distr)), d_test$params))
          expect_gt(ks_test$p.value, p_alpha) # p-value > 0.05 indicates a good fit # @ output variable and warning about test
        }
      } else {
        # Sex-specific KS test
        l_v_samp <- list(
          male = v_samp[v_male == 1],
          female = v_samp[v_male == 0]
        )
        for (sex in names(d_test)) {
          ks_test <- do.call(ks.test, c(x = list(l_v_samp[[sex]]), y = list(paste0("p", d_test[[sex]]$distr)), d_test[[sex]]$params))
          expect_gt(ks_test$p.value, p_alpha) # p-value > 0.05 indicates a good fit # @ output variable and warning about test
        }
      }
    })
  }
}

# Test correctness of number of lesions
test_that("Test correctness of number of lesions", {
  # Test equality of number of patients with lesions
  expect_equal(nrow(res$patient_level[!is.na(n_L_add)]), length(unique(res$lesion_level$pt_id)))
  
  # Compare number of lesions from patient- and lesion-level data table
  m_n_lesion <- merge(res$patient_level[, .SD, .SDcols = c("pt_id", "n_L_add")], res$lesion_level[, .(n_L = .N), by = pt_id], by = "pt_id")
  expect_equal(m_n_lesion$n_L, m_n_lesion$n_L_add + 1)
})

# Test addition of quantities (assuming 4 stages of cancer)
test_that("Test addition of quantities", {
  # Time to preclinical cancer, lesion-level
  expect_equal(res$lesion_level[, time_H_Pj], res$lesion_level[, time_H_L + time_L_Lj + time_Lj_Pj])
  
  # Time to death from cancer
  expect_equal(res$patient_level[, time_H_Dc], res$patient_level[, time_H_P + time_P_C + time_C_Dc])
  
  # Stage at cancer diagnosis
  res$patient_level[time_H_P < time_H_Do, stage_dx_test := fcase(time_P1_C1 < time_P1_P2, 1,
                                                                 time_P2_C2 < time_P2_P3, 2,
                                                                 time_P3_C3 < time_P3_P4, 3,
                                                                 default = 4)]
  expect_equal(res$patient_level[, stage_dx], res$patient_level[, stage_dx_test])
  
  # Time to cancer diagnosis
  res$patient_level[time_H_P < time_H_Do, time_P_C_test := fcase(time_P1_C1 < time_P1_P2, time_P1_C1,
                                                                 time_P2_C2 < time_P2_P3, time_P1_P2 + time_P2_C2,
                                                                 time_P3_C3 < time_P3_P4, time_P1_P2 + time_P2_P3 + time_P3_C3,
                                                                 time_P3_C3 > time_P3_P4, time_P1_P2 + time_P2_P3 + time_P3_P4 + time_P4_C4)]
  expect_equal(res$patient_level[, time_P_C], res$patient_level[, time_P_C_test])
})

# Test minimum of quantities
test_that("Test minimum of quantities", {
  # Time to preclinical cancer, lesion-level
  expect_equal(res$lesion_level[, time_H_Pj], res$lesion_level[, time_H_L + time_L_Lj + time_Lj_Pj])
  
  # Time to death from cancer
  expect_equal(res$patient_level[time_H_L < time_H_Do, time_H_P], res$lesion_level[, .(min(time_H_Pj)), by = pt_id]$V1)
  
  # Time to death from any cause
  expect_equal(res$patient_level[, time_H_D], res$patient_level[, pmin(time_H_Do, time_H_Dc, na.rm = T)])
})

