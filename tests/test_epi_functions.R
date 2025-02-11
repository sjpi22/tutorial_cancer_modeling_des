# Cross-sectional versus longitudinal prevalence

df_true_params <- readRDS("ground_truth/true_param_map.rds")
l_params <- load_model_params()
l_params <- update_param_from_map(l_params,
                                  df_true_params$param_val,
                                  df_true_params)
l_params$seed <- NULL
set.seed(123)

l_res_long <- list()
l_res_cs <- list()
n_sim <- 30
v_ages <- seq(20, 80, 10)
for (n_cohort in c(10000, 100000, 500000)) {
  # Create matrix placeholders
  m_long <- matrix(nrow = n_sim, ncol = length(v_ages) - 1)
  m_cs <- matrix(nrow = n_sim, ncol = length(v_ages) - 1)
  for (sim in 1:n_sim) {
    # Run model
    res <- run_base_model(l_params)
    
    # Calculate longitudinal prevalence
    temp_prev_long <- calc_prevalence(res, "time_H_P", "time_H_C", "time_H_D",
                                      v_ages = v_ages)
    
    # Append to matrix
    m_long[sim, ] <- temp_prev_long$value
    
    
    # Calculate cross-sectional prevalence by age group, uniform sample
    res[, sample_age := runif(.N, min = 20, max = 80)]
    res[, sample_grp := floor(sample_age/10)*10]
    res[, flg_case := (time_H_P <= sample_age & time_H_C > sample_age)]
    temp_prev_cs <- res[time_H_D > sample_age, mean(flg_case), by = sample_grp]
    temp_prev_cs <- temp_prev_cs[order(sample_grp)]
    m_cs[sim, ] <- temp_prev_cs$V1
  }
  l_res_long <- c(l_res_long, list(m_long))
  l_res_cs <- c(l_res_cs, list(m_cs))
}

# Plot the results
colMeans(l_res_long[[3]])
apply(l_res_long[[3]], 2, sd)
colMeans(l_res_cs[[3]])
apply(l_res_cs[[3]], 2, sd)



