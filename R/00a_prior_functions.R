# prior functions

# Fit spline
fit_spline <- function(x_vals,
                       y_vals,
                       wt = 1,
                       constraints = "none",
                       v_knots = NULL
){
  # Set knots
  if (is.null(v_knots)) {
    # Sample every three values of ages for spline knots
    v_knots <- c(x_vals[seq(1, length(x_vals) - 1, 3)], max(x_vals))
  }
  
  # Set weights
  if (length(wt) != length(x_vals)) {
    wt <- rep(wt, length(x_vals))
  }
  
  # Fit spline
  spline_res <- cobs(
    x = x_vals,
    y = y_vals,
    w = wt,
    constraint = constraints,
    knots = v_knots
  )
  
  return(spline_res)
}

# Convert incidence rate to cumulative incidence
incidence_rate_to_cdf <- function(x_vals,
                                  y_vals,
                                  censored_at_event = TRUE,
                                  distr_censor = NULL, # If not censored at event, time-to-event distribution of censoring from event
                                  wt = 1,
                                  constraints = "none",
                                  v_knots = NULL,
                                  x_pred = NULL) {
  # Fit spline to incidence rate
  spline_incidence <- fit_spline(x_vals = x_vals,
                                 y_vals = y_vals,
                                 wt = wt,
                                 constraints = constraints,
                                 v_knots = v_knots)
  
  # Set x values for prediction
  if (is.null(x_pred)) {
    x_pred <- x_vals
  }
  
  # If censored at time of event, treat rate as hazard
  if (censored_at_event) {
    # Integrate spline to estimate cumulative hazard at ages of interest
    v_chaz <- sapply(x_pred,
                     function(t) integrate(function(x) pmax(0, predict(spline_incidence, x)[, "fit"]), 
                                           lower = 0, upper = t)[["value"]])
    
    # Calculate cumulative probability
    v_cdf <- 1 - exp(-v_chaz)
  } else if (!is.null(distr_censor)) {
    # Integrate spline to estimate cumulative probability
    v_cdf <- sapply(x_vals,
                    function(t) integrate(function(x) pmax(0, predict(spline_incidence, x)[, "fit"]), 
                                          lower = 0, upper = t)[["value"]])
    
    # Estimate percent censored
    v_cdf_censored <- sapply(x_vals,
                             function(t) integrate(function(u) pmax(0, predict(spline_incidence, u)[, "fit"])*query_distr("p", t - u, distr_censor$distr, distr_censor$params), 
                                                   lower = 0, upper = t)[["value"]])
    
    # Recalculate probability density, scaling by percent censored
    v_pdf <- y_vals * (1 - v_cdf_censored)
    
    # Refit spline to probability density
    spline_pdf <- cobs(
      x = x_vals,
      y = v_pdf,
      w = wt,
      constraint = constraints,
      knots = v_knots
    )
    
    # Calculate cumulative probability of by integrating over density
    v_cdf <- sapply(x_pred,
                    function(t) integrate(function(x) pmax(0, predict(spline_pdf, x)[, "fit"]), 
                                          lower = 0, upper = t)[["value"]])
  } else {
    # Integrate spline to estimate cumulative probability at ages of interest
    v_cdf <- sapply(x_pred,
                    function(t) integrate(function(x) pmax(0, predict(spline_incidence, x)[, "fit"]), 
                                          lower = 0, upper = t)[["value"]])
  }
  return(v_cdf)
}

# Fit Weibull parameters to a cumulative distribution function
fit_weibull <- function(x_vals,
                        y_vals,
                        wt = 1) {
  # Calculate Weibull x transformation
  x_transformed = log(x_vals)
  
  # Calculate Weibull y transformation
  y_transformed <- log(-log(1 - y_vals))
  
  # Exclude non-finite values
  nonzero_indices <- which(is.finite(x_transformed) & is.finite(y_transformed))
  x_transformed <- x_transformed[nonzero_indices]
  y_transformed <- y_transformed[nonzero_indices]
  
  # Set weights
  if (length(wt) != length(x_transformed)) {
    wt <- rep(wt, length(x_transformed))
  } else {
    wt <- wt[nonzero_indices]
  }
  
  # Get Weibull estimates using weighted linear regression
  weibull_lm <- lm(y_transformed ~ x_transformed,
                   weights = wt)
  
  # Get shape and scale estimates for time to disease onset distribution
  params <- list(
    shape = weibull_lm$coefficients[2],
    scale = exp(-weibull_lm$coefficients[1]/weibull_lm$coefficients[2])
  )
  return(params)
}


# Fit log-logistic parameters to a cumulative distribution function
fit_llogis <- function(x_vals,
                       y_vals,
                       wt = 1) {
  # Calculate log-logistic x transformation
  x_transformed = log(x_vals)
  
  # Calculate log-logistic y transformation
  y_transformed <- log(y_vals/(1 - y_vals))
  
  # Exclude non-finite values
  nonzero_indices <- which(is.finite(x_transformed) & is.finite(y_transformed))
  x_transformed <- x_transformed[nonzero_indices]
  y_transformed <- y_transformed[nonzero_indices]
  
  # Set weights
  if (length(wt) != length(x_transformed)) {
    wt <- rep(wt, length(x_transformed))
  } else {
    wt <- wt[nonzero_indices]
  }
  
  # Get Weibull estimates using weighted linear regression
  llogis_lm <- lm(y_transformed ~ x_transformed,
                  weights = wt)
  
  # Get shape and scale estimates for time to disease onset distribution
  params <- list(
    shape = llogis_lm$coefficients[2],
    scale = exp(-llogis_lm$coefficients[1]/llogis_lm$coefficients[2])
  )
  return(params)
}


# Perform deconvolution to calculate PDF of addend variable
deconvolve <- function(x_target,
                       y_target,
                       fn_cdf_1, # Function to calculate CDF of time to first state
                       delta = 1,
                       penalty_jump = 0,
                       penalty_flip = 0) {
  # Set times to evaluate time to first state
  v_times_1 <- seq(0, max(x_target), by = delta)
  
  # Get values of CDF of time to first state
  v_cdf_1 <- do.call(fn_cdf_1, list(v_times_1))
  
  # Build matrix for deconvolution with CDF of onset time
  m_cdf <- matrix(0, nrow = length(x_target), ncol = length(v_cdf_1) - 1)
  for (i in 1:nrow(m_cdf)) {
    idx <- x_target[i]/delta # Get maximum index of matrix to fill
    m_cdf[i, 1:idx] <- rev(v_cdf_1[1:idx] + v_cdf_1[2:(idx + 1)])/2 # Trapezoid method
  }
  
  # Convex optimization, least squares penalty
  f_12 <- Variable(ncol(m_cdf))
  objective <- Minimize(sum((y_target - delta*m_cdf %*% f_12)^2)+ 
                          penalty_jump * sum((diff(f_12))^2) + # Penalize large differences in f_12
                          penalty_flip * sum((diff(f_12)[1:(ncol(m_cdf)-2)] - diff(f_12)[2:(ncol(m_cdf)-1)])^2)) # Penalize when direction of f_12 changes
  constraint1 <- f_12 >= 0
  constraint2 <- sum(f_12) <= 1
  problem <- Problem(objective, constraints = list(constraint1, constraint2))
  result <- solve(problem)
  f_12_vals <- result$getValue(f_12)
  
  return(list(x = seq(delta, max(x_target), by = delta) - delta/2,
              pdf = f_12_vals))
}
