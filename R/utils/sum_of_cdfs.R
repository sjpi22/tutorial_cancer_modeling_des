# Function to calculate CDF of the sum of two distributions at various indices
calc_cdf_sum <- function(v_idx, d_time_1, d_time_2) {
  # Calculate predicted CDF at indices
  p_pred <- sapply(v_idx, function(t) {
    # Calculates sum of the two distributions
    integrate(function(u)
      query_distr("p", t - u, d_time_1$distr, d_time_1$params) * 
        query_distr("d", u, d_time_2$distr, d_time_2$params),
      lower = 0, upper = t)[["value"]]
  })
  return(p_pred)
}

# Create objective function for difference between predicted and observed CDF
obj_fn_cdf <- function(l_params_2, d_time_1, d_time_2, 
                       idx_target, p_target, sd_target = NULL,
                       method = "ss") {
  # Update distribution parameters
  d_time_2$params <- modifyList(d_time_2$params, l_params_2)
  
  # Calculate CDF of sum at indices
  p_pred <- calc_cdf_sum(idx_target, d_time_1, d_time_2)
  
  # Calculate sum of squares or log-likelihood depending on "method" argument
  if (method == "ss") {
    # Calculate sum of squares of observed vs estimated clinical cancer CDF
    if (is.null(sd_target)) sd_target <- 1
    v_diff_cdf <- (p_target - p_pred) / sd_target
    
    # Return negative total sum of squares
    return(-sum(v_diff_cdf^2))
  } else if (method == "ll") {
    # Calculate log likelihood of observed vs estimated clinical cancer CDF
    v_diff_cdf <- dnorm(x = p_target, mean = p_pred, sd = sd_target, log = T) 
    
    # Return total log-likelihood
    return(sum(v_diff_cdf))
  } else {
    stop("Invalid method")
  }
}
