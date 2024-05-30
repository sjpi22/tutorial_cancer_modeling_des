# Parametrized functions for hazard, integral of hazard, and cumulative hazard of detection
# x is current level
# a is parameter that increases with rate of early detection
# x_D is fatal level (default 1)
hazard_detect <- function(x, a, x_D = 1) {
  val <- sapply(a/x_D * x, function(u) {
    if(u >= a) {
      return(Inf)
    } else return(a / (a - u) - 1)
  }
  )
  return(val)
}

integral_hazard_detect <- function(x, a, x_D = 1) {
  val <- sapply(a/x_D * x, function(u) {
    if(u >= a) {
      return(Inf)
    } else return(-a * log(a - u) - u)
  }
  )
  return(val)
}

cumulative_hazard_detect <- function(x, a, x_D = 1, x_lower = 0) {
  c_haz <- integral_hazard_detect(x, a, x_D) - 
    integral_hazard_detect(x_lower, a, x_D)
  return(c_haz)
}



#' Density from hazard
#'
#' @param x Vector of quantiles
#' @param haz_rate Hazard rate as function of quantile
#' @param haz_cumulative Cumulative hazard as function of quantile
#' @param log Logical; if TRUE, probabilities p are given as log(p)
#' 
#' @return Density at x
#' 
#' @export
dhazard <- function(x, haz_rate = NULL, haz_cumulative = NULL, log = FALSE, 
                    ...) {
  # Get hazard and cumulative hazard at x
  if(!is.null(haz_cumulative)) {
    H <- match.fun(haz_cumulative)(x)
    if(!is.null(haz_rate)) {
      h <- match.fun(haz_rate)(x)
    } else {
      h <- Deriv(haz_cumulative)(x)
    }
  } else if(!is.null(haz_rate)) {
    H <- integrate(function(u) match.fun(haz_rate)(u), 
                   lower = 0, 
                   upper = x)$value
    h <- match.fun(haz_rate)(x)
  } else stop("Provide function for cumulative hazard or hazard rate")
  
  # Convert to density (derivative of CDF at time point)
  val <- h*exp(-H)
  
  # Log transform density if necessary
  if(log == TRUE) val <- log(val)
  return(val)
}

#' Cumulative distribution function from hazard
#'
#' @param q Vector of quantiles
#' @param haz_rate Hazard rate as function of quantile; either \code{haz_rate} 
#' or \code{haz_cumulative} can be provided but \code{haz_cumulative} is 
#' preferred to prevent integrating, reducing calculation time
#' @param haz_cumulative Cumulative hazard as function of quantile; either 
#' \code{haz_rate} or \code{haz_cumulative} can be provided but 
#' \code{haz_cumulative} is preferred to prevent integrating, reducing 
#' calculation time
#' @param lower.tail logical; if TRUE (default), probabilities are P(X<=x), 
#' otherwise P(X>x)
#' @param log.p Logical; if TRUE, probabilities p are given as log(p)
#' 
#' @return Cumulative distribution function evaluated at q
#' 
#' @export
phazard <- function(q, haz_rate = NULL, haz_cumulative = NULL, 
                    lower.tail = TRUE, log.p = FALSE, ...) {
  # Get cumulative hazard at q
  if(!is.null(haz_cumulative)) {
    H <- match.fun(haz_cumulative)(q)
  } else if(!is.null(haz_rate)) {
    H <- integrate(function(u) match.fun(haz_rate)(u), 
                   lower = 0, 
                   upper = q)$value
  } else stop("Provide function for cumulative hazard or hazard rate")
  
  # Convert to CDF
  val <- 1 - exp(-H)
  
  # Right tail and log transform if necessary
  if(lower.tail == FALSE) val <- 1 - val
  if(log.p) val <- log(val)
  return(val)
}

#' Quantile function of distribution from hazard
#'
#' @param p Vector of probabilities
#' @param haz_rate Hazard rate as function of quantile; either \code{haz_rate} 
#' or \code{haz_cumulative} can be provided but \code{haz_cumulative} is 
#' preferred to prevent integrating, reducing calculation time
#' @param haz_cumulative Cumulative hazard as function of quantile; either 
#' \code{haz_rate} or \code{haz_cumulative} can be provided but 
#' \code{haz_cumulative} is preferred to prevent integrating, reducing 
#' calculation time
#' @param lower.tail logical; if TRUE (default), probabilities are P(X<=x), 
#' otherwise P(X>x)
#' @param log.p Logical; if TRUE, probabilities p are given as log(p)
#' @param quantiles Precalculated quantiles to bound search
#' @param probs Corresponding precalculated probabilities of quantiles to bound 
#' search
#' 
#' @return Quantile of cumulative probability p
#' 
#' @export
qhazard <- function(p, haz_rate = NULL, haz_cumulative = NULL,
                    lower.tail = TRUE, log.p = FALSE, quantiles = NULL,
                    probs = NULL, ...) {
  # Check for hazard functions
  if(is.null(haz_cumulative) & is.null(haz_rate)) {
    stop("Provide function for cumulative hazard or hazard rate")
  } 
  
  # Right tail and log transform if necessary
  if(log.p) p <- exp(p)
  if(lower.tail == FALSE) p <- 1 - p
  
  # Get bounds for search if not provided already
  if(is.null(quantiles)) { 
    res <- bounds_hazard(p, haz_rate, haz_cumulative,
                         process_hazards = FALSE)
    quantiles <- res$quantiles
    probs <- res$probs
  } else { # Check if bounds cover range of inputs
    # Find min and max values of inputted probabilities within boundaries
    min_p <- min(p)
    max_p <- max(p[p < 1])
    
    # Extend quantiles if range of probabilities does not cover range of inputs
    if(max_p > max(probs) | min_p < min(probs)) {
      res <- bounds_hazard(p, haz_rate, haz_cumulative, quantiles = quantiles,
                           probs = probs, process_hazards = FALSE, ...)
      quantiles <- res$quantiles
      probs <- res$probs
    }
  }
  
  # Convert probs to intervals
  prob_index <- findInterval(p, probs)
  quantiles_lb <- c(0, quantiles)
  min_q <- quantiles_lb[prob_index + 1]
  max_q <- quantiles[prob_index + 1]
  
  # Find time value where cumulative distribution function equals probability
  if(!is.null(haz_cumulative)) {
    val <- sapply(seq(length(p)), function(i) {
      if(p[i] == 0) {
        return(0)
      } else if(p[i] == 1) {
        return(Inf)
      } else {
        uniroot(function(q) match.fun(haz_cumulative)(q) + log(1 - p[i]),
                lower = min_q[i],
                upper = max_q[i])$root}
    })
  } else if(!is.null(haz_rate)) {
    val <- sapply(seq(length(p)), function(i) {
      if(p[i] == 0) {
        return(0)
      } else if(p[i] == 1) {
        return(Inf)
      } else {
        uniroot(function(q) integrate(function(u) match.fun(haz_rate)(u), 
                                      lower = 0, 
                                      upper = q)$value + log(1 - p[i]),
                lower = min_q[i],
                upper = max_q[i])$root}
    })
  }
  
  # Return quantile
  return(val)
}

#' Generates n random values from distribution from hazard
#'
#' @param n Number of observations
#' @param haz_rate Hazard rate as function of quantile; either \code{haz_rate} 
#' or \code{haz_cumulative} can be provided but \code{haz_cumulative} is 
#' preferred to prevent integrating, reducing calculation time
#' @param haz_cumulative Cumulative hazard as function of quantile; either 
#' \code{haz_rate} or \code{haz_cumulative} can be provided but 
#' \code{haz_cumulative} is preferred to prevent integrating, reducing 
#' calculation time
#' @param return_prob Logical; if TRUE, return P(X<=x)
#' 
#' @return Vector of randomly generated values
#' 
#' @export
rhazard <- function(n, haz_rate = NULL, haz_cumulative = NULL, 
                    return_prob = FALSE, ...) {
  # Sample probability for inverse CDF transform
  prob <- runif(n)
  
  # Get quantile corresponding to sampled probability
  val <- qhazard(prob, haz_rate, haz_cumulative)
  
  if(return_prob) {
    return(list(val = val, prob = prob))
  } else return(val)
}


# @@@ Calculate vector of quantiles and corresponding probabilities
bounds_hazard <- function(p, haz_rate = NULL, haz_cumulative = NULL,
                    lower.tail = TRUE, log.p = FALSE, 
                    quantiles = NULL, probs = NULL, 
                    process_hazards = TRUE, recalculate_probs = TRUE,
                    n_min_reps = 3, ...) {
  if(process_hazards) {
    # Check for hazard functions
    if(is.null(haz_cumulative) & is.null(haz_rate)) {
      stop("Provide function for cumulative hazard or hazard rate")
    } 
    
    # Right tail and log transform if necessary
    if(log.p) p <- exp(p)
    if(lower.tail == FALSE) p <- 1 - p
  }
  
  # Get upper bound of probability
  if(is.null(quantiles)) { 
    # Generate initial sequence of quantiles and probs
    quantiles <- 2^seq(0, 8)
    probs <- phazard(quantiles, haz_rate, haz_cumulative)
  } else if(recalculate_probs) { # Update probabilities to distribution
    probs <- phazard(quantiles, haz_rate, haz_cumulative)
  }
  
  # Find min and max values of inputted probabilities within boundaries
  min_p <- min(p)
  max_p <- max(p[p < 1])
    
  # Extend quantiles if range of probabilities does not cover range of inputs
  if(max_p > max(probs)) {
    # Test larger quantiles until the probabilities cover the range of inputs
    q_candidates <- max(quantiles) * c(2, 4, 8)
    p_candidate <- max(probs)
    while(max_p > p_candidate) {
      # Calculate associated probabilities
      p_candidates <- phazard(q_candidates, haz_rate, haz_cumulative)
      
      # Append new quantiles/probabilities
      quantiles <- c(quantiles, q_candidates)
      probs <- c(probs, p_candidates)
      
      # Get new quantile and probability candidates
      q_candidates <- max(q_candidates) * c(2, 4, 8)
      p_candidate <- max(p_candidates)
    }
  }
  
  if(min_p < min(probs)) {
    # Test smaller quantiles until the probabilities cover the range of inputs
    q_candidates <- min(quantiles) * c(1/8, 1/4, 1/2)
    p_candidate <- min(probs)
    n_reps <- 1
    while(min_p < p_candidate & n_reps < n_min_reps) {
      # Calculate associated probabilities
      p_candidates <- phazard(q_candidates, haz_rate, haz_cumulative)
      
      # Append new quantiles/probabilities
      quantiles <- c(q_candidates, quantiles)
      probs <- c(p_candidates, probs)
      
      # Get new quantile and probability candidates
      q_candidates <- min(q_candidates) * c(1/8, 1/4, 1/2)
      p_candidate <- min(p_candidates)
      n_reps <- n_reps + 1
    }
    
    # If still does not cover, add 0
    if(min_p < p_candidate) {
      # Append new quantiles/probabilities
      quantiles <- c(0, quantiles)
      probs <- c(0, probs)
    }
  }
  
  # Return quantile
  return(list(quantiles = quantiles, probs = probs))
}
