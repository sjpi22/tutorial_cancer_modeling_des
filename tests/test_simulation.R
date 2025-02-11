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

set.seed(123)
eps_pct <- 0.01 # tolerated error
alpha <- 0.05 # p-value threshold


#### 3. Pre-processing actions  ===========================================

# Load default parameters
l_params_init <- load_model_params(
  lesion_state = TRUE
)

# Make parameter map and change parameters to be more realistic
df_params <- make_param_map(l_params_init)
df_params$param_val <- c(2, 75, 1/20, 1/10, 1/7, 1/35, 1/6, 1/16, 1/5, 1/10, 1/2)

# Update parameter list
l_params_model <- update_param_from_map(l_params_init, df_params$param_val, df_params)

# Simulate natural history model
res <- run_base_model(l_params_model)


#### 4. Tests  ===========================================

# Test that time to event distributions are simulated correctly
for (varname in unique(df_params$var_name)) {
  # Filter to variable
  df_var <- df_params[df_params$var_name == varname, ]
  
  # Create list of parameters
  l_var_params <- as.list(df_var$param_val)
  names(l_var_params) <- df_var$param_name

  # Get variable distribution
  d_test <- list(
    distr = df_var$var_distr[1],
    params = l_var_params
  )
  
  # Extract variable removing NAs
  if (varname == "d_time_L_P") { # Special case - lesion-level variable with different name
    v_samp <- res$lesion_level[["time_Lj_Pj"]]
  } else if (varname == "d_n_L") {
    v_samp <- res$patient_level[["n_L_add"]]
  } else {
    v_samp <- res$patient_level[[substr(varname, 3, nchar(varname))]]
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
      expect_gt(chisq_test$p.value, alpha) # p-value > 0.05 indicates a good fit
    })
    
    # Change to check distribution of lesion onset time
    d_test <- list(
      distr = "unif",
      params = list(min = 0, max = 1)
    )
    
    # Extract scaled time to lesion onset
    v_samp <- res$lesion_level[lesion_id > 1, .(pct = time_L_Lj / (time_H_Do - time_H_L))]
    v_samp <- v_samp$pct
  } 
  
  # Define theoretical parameters
  expected_mean <- integrate(function(x) f_mean(x, d_test), 0, Inf)$value
  expected_sd <- sqrt(integrate(function(x) f_var(x, d_test, expected_mean), 0, Inf)$value)
  
  # Set tolerance as 1% of expected mean
  # Alternative using sampling error: tol_var <- qnorm(1 - alpha/2) * expected_sd / sqrt(length(v_samp))
  tol_var <- expected_mean * eps_pct
  
  test_that("Simulated sample matches distribution", {
    # Check mean
    sample_mean <- mean(v_samp)
    expect_equal(sample_mean, expected_mean, tolerance = tol_var)
    
    # Check standard deviation
    sample_sd <- sd(v_samp)
    expect_equal(sample_sd, expected_sd, tolerance = tol_var)
    
    # Kolmogorov-Smirnov test for goodness-of-fit
    ks_test <- do.call(ks.test, c(x = list(v_samp), y = list(paste0("p", d_test$distr)), d_test$params))
    expect_gt(ks_test$p.value, alpha) # p-value > 0.05 indicates a good fit # @ output variable and warning about test
  })
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
  res$patient_level[time_H_P < time_H_Do, stage_dx_test := fcase(time_P1_C < time_P1_P2, 1,
                                                                 time_P2_C < time_P2_P3, 2,
                                                                 time_P3_C < time_P3_P4, 3,
                                                                 default = 4)]
  expect_equal(res$patient_level[, stage_dx], res$patient_level[, stage_dx_test])
  
  # Time to cancer diagnosis
  res$patient_level[time_H_P < time_H_Do, time_P_C_test := fcase(time_P1_C < time_P1_P2, time_P1_C,
                                                                 time_P2_C < time_P2_P3, time_P1_P2 + time_P2_C,
                                                                 time_P3_C < time_P3_P4, time_P1_P2 + time_P2_P3 + time_P3_C,
                                                                 time_P3_C > time_P3_P4, time_P1_P2 + time_P2_P3 + time_P3_P4 + time_P4_C)]
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

