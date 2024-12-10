# # Set screening parameters
# screen_age <- 50
# screen_sens_early <- 0.7
# screen_sens_late <- 0.85
# screen_spec <- 0.8
# 
# # Make copy of no screening data
# m_screening <- copy(m_no_screening)
# 
# # Flag if eligible for screening - have not died or been diagnosed with cancer before screening age
# m_screening[, screen_elig := (time_H_D > screen_age) & (is.na(time_H_C) | time_H_C > screen_age)]
# 
# # Among screen eligible, flag if patient has cancer
# m_screening[screen_elig == TRUE, screen_has_cancer := (time_H_P < screen_age)]
# 
# # Among screen eligible patients with cancer, flag stage at screening age
# m_screening[screen_has_cancer == TRUE, screen_stage_dx := ifelse(time_H_P + time_P1_P2 < screen_age, 2, 1)]
# 
# # Sample probability of detection among people who have developed early-stage cancer (TP = true positive)
# m_screening[screen_stage_dx == 1, screen_TP := rbinom(.N, size = 1, prob = screen_sens_early)]
# 
# # Sample probability of detection among people who have developed late-stage cancer (TP = true positive)
# m_screening[screen_stage_dx == 2, screen_TP := rbinom(.N, size = 1, prob = screen_sens_late)]
# 
# # Update stage and age of detection for people who were screen-detected
# m_screening[screen_TP == TRUE, stage_dx := screen_stage_dx]
# m_screening[screen_TP == TRUE, time_H_C := screen_age]
# 
# # Sample survival among people with screen detected stage 1 cancer
# m_screening[(screen_TP == TRUE) & (screen_stage_dx == 1), time_C_Dc_screen := sample(x = ages_C1_Dc, 
#                                                                                      size = .N, 
#                                                                                      replace = TRUE, 
#                                                                                      prob = prob_C1_Dc) + runif(.N)]
# 
# # Reset to cured if maximum age was sampled
# m_screening[(screen_TP == TRUE) & (screen_stage_dx == 1) & (time_C_Dc_screen > max(ages_C1_Dc)), time_C_Dc_screen := Inf]
# 
# # Sample survival among people with screen detected stage 2 cancer
# m_screening[(screen_TP == TRUE) & (screen_stage_dx == 2), time_C_Dc_screen := sample(x = ages_C2_Dc, 
#                                                                                      size = .N, 
#                                                                                      replace = TRUE, 
#                                                                                      prob = prob_C2_Dc) + runif(.N)]
# 
# # Reset to cured if maximum age was sampled
# m_screening[(screen_TP == TRUE) & (screen_stage_dx == 2) & (time_C_Dc_screen > max(ages_C2_Dc)), time_C_Dc_screen := Inf]
# 
# # Do not allow death in screen-detected scenario to be earlier than death in non-screening scenario
# m_screening[screen_TP == TRUE, time_C_Dc := pmax(time_C_Dc, time_C_Dc_screen)]
# 
# # Calculate time to death from cancer
# m_screening[, time_H_Dc := time_H_C + time_C_Dc]
# 
# # Sample probability of false positive among people who do not have cancer (TN = true negative)
# m_screening[screen_has_cancer == FALSE, screen_FP := rbinom(.N, size = 1, prob = 1 - screen_spec)]
# 
# # Calculate mortality outcomes
# calc_mortality_outcomes(m_screening)