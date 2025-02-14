# Run checks after simulating cohort in simulate_truth.R
# Line: m_cohort <- run_base_model(l_params_all)

#### Consistency checks
# Check that population is ~50% male
assert_that(abs(mean(m_cohort$male) - 0.5) < 0.001) 

# Check time to preclinical cancer
plot(0:l_params_all$max_age, pweibull(0:l_params_all$max_age, shape=v_param_update[1], scale=v_param_update[2]), 
     type='l', xlab="Age", ylab="% with cancer", main="Cancer onset")
lines(sort(m_cohort[time_H_P < l_params_all$max_age, time_H_P]), 1:length(m_cohort[time_H_P < l_params_all$max_age, time_H_P])/n_cohort, col='red')

# Check time to progression distributions
k <- 3
v_cols <- c("blue", "green", "orange", "red")
for (i in head(l_params_all$v_cancer, -1)) {
  var_progress <- paste0("time_P", i, "_P", i+1)
  if (i == 1) {
    plot(0:10, pexp(0:10, rate=v_param_update[k]), type='l', 
         xlab="Time from beginning of stage", ylab="% detected", 
         main="Preclinical cancer progression", lty=i, ylim=c(0, 1))
  } else {
    lines(0:10, pexp(0:10, rate=v_param_update[k]), lty=i)
  }
  lines(sort(m_cohort[!is.na(get(var_progress)), get(var_progress)]), 
        1:length(m_cohort[!is.na(get(var_progress)), get(var_progress)])/nrow(m_cohort[!is.na(get(var_progress))]), 
        col=v_cols[i], lty=i)
  k <- k+2
}

# Check time to detection distributions
k <- 4
for (i in l_params_all$v_cancer) {
  var_detect <- paste0("time_P", i, "_C")
  if (i == 1) {
    plot(0:10, pexp(0:10, rate=v_param_update[k]), type='l', 
         xlab="Time from beginning of stage", ylab="% progressed", 
         main="Cancer detection", lty=i, ylim=c(0, 1))
  } else {
    if (i == 4) k <- length(v_param_update)
    lines(0:10, pexp(0:10, rate=v_param_update[k]), lty=i)
  }
  lines(sort(m_cohort[!is.na(get(var_detect)), get(var_detect)]), 
        1:length(m_cohort[!is.na(get(var_detect)), get(var_detect)])/nrow(m_cohort[!is.na(get(var_detect))]), 
        col=v_cols[i], lty=i)
  k <- k+2
}

# Check time death by stage
for (i in l_params_all$v_cancer) {
  if (i == 1) {
    plot(0:10, pexp(0:10, rate=rate_Cx_Dc[i]), type='l', 
         xlab="Time from detection", ylab="% died", 
         main="Disease-specific survival", lty=i)
  } else {
    lines(0:10, pexp(0:10, rate=rate_Cx_Dc[i]), lty=i)
  }
  lines(sort(m_cohort[stage_dx == i, time_C_Dc]), 
        1:nrow(m_cohort[stage_dx == i])/nrow(m_cohort[stage_dx == i]), 
        col=v_cols[i], lty=i)
}

# Plot time to death by cause
plot(sort(sample(m_cohort[, time_H_D], 1000)), 1:1000/1000, type='l', xlab="Age", ylab="% dead", main="Time to death by cause")
lines(sort(sample(m_cohort[male == 0, time_H_D], 1000)), 1:1000/1000, lty=2)
lines(sort(sample(m_cohort[male == 1, time_H_D], 1000)), 1:1000/1000, lty=3)
lines(sort(sample(m_cohort[, time_H_Do], 1000)), 1:1000/1000, col='blue')
lines(sort(sample(m_cohort[!is.na(time_H_Dc), time_H_Dc], 1000)), 
      1:1000/1000*nrow(m_cohort[!is.na(time_H_Dc)])/n_cohort, col='red')

# Check stage at detection calculations
assert_that(sum(m_cohort[time_P1_C < time_P1_P2, pt_id] != m_cohort[stage_dx == 1, pt_id])==0)
assert_that(sum(m_cohort[time_P1_C > time_P1_P2 & 
                           time_P2_C < time_P2_P3, pt_id] != m_cohort[stage_dx == 2, pt_id])==0)
assert_that(sum(m_cohort[time_P1_C > time_P1_P2 & 
                           time_P2_C > time_P2_P3 &
                           time_P3_C < time_P3_P4, pt_id] != m_cohort[stage_dx == 3, pt_id])==0)
assert_that(sum(m_cohort[time_P1_C > time_P1_P2 & 
                           time_P2_C > time_P2_P3 &
                           time_P3_C > time_P3_P4, pt_id] != m_cohort[stage_dx == 4, pt_id])==0)

# Check time detection calculations
assert_that(sum(m_cohort[time_P1_C < time_P1_P2, time_P_C != time_P1_C])==0)
assert_that(sum(m_cohort[time_P1_C > time_P1_P2 & 
                           time_P2_C < time_P2_P3, 
                         time_P_C != time_P1_P2 + time_P2_C])==0)
assert_that(sum(m_cohort[time_P1_C > time_P1_P2 & 
                           time_P2_C > time_P2_P3 &
                           time_P3_C < time_P3_P4, 
                         time_P_C != time_P1_P2 + time_P2_P3 + time_P3_C])==0)
assert_that(sum(m_cohort[time_P1_C > time_P1_P2 & 
                           time_P2_C > time_P2_P3 &
                           time_P3_C > time_P3_P4,
                         time_P_C != time_P1_P2 + time_P2_P3 + time_P3_P4 + time_P4_C])==0)
