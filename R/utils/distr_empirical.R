#' Density of empirical distribution
#'
#' @param x Vector of quantiles
#' @param xs Vector of possible values; if continuous, when sorted, each value 
#' of xs represents the lower limit of the interval to which the corresponding 
#' probability applies, while the next value of \code{xs} represents the upper 
#' bound; if \code{NULL}, generate from indices of \code{probs}
#' @param probs Vector of probabilities associated with the corresponding value 
#' in \code{xs}
#' @param warn_sum_probs Logical; if true, print warning if sum of probabilities 
#' is not 1; can be silenced if truncated distribution is used as a conditional 
#' distribution
#' @param max_x Maximum value of \code{xs} for continuous variables as upper bound for 
#' the interval lower bounded by the last value of \code{xs}
#' @param log Logical; if TRUE, probabilities p are given as log(p)
#' 
#' @return Density at xs
#' 
#' @export
dempirical <- function(x, xs = NULL, probs, max_x = NULL, warn_sum_probs = TRUE, 
                       log.p = FALSE, eps = 1e-6) {
  # Sanity checks for empirical distribution
  check_empirical(xs, probs, warn_sum_probs, eps)
  
  # Set xs if NULL to indices of probs, else sort if needed
  if(length(xs) == 0) {
    xs <- seq(length(probs))
  } else if(is.unsorted(xs)) {
    indices <- base::order(xs)
    xs <- xs[indices]
    probs <- probs[indices]
  }
  
  # If maximum value of xs is given, append it for correction
  if (!is.null(max_x)) {
    xs <- c(xs, max_x)
  }
  
  # Find index of largest value of xs for which xs >= x and find probability
  # mass at that point
  lb_index <- sapply(x, function(x) which.max(xs[x >= xs]))
  
  # Set uniform probability density if probability and bounds exist at point
  val <- rep(NA, length(x))
  val[(lb_index < length(xs)) & (lb_index > 0)] <- probs[lb_index] / (xs[lb_index + 1] - xs[lb_index])
  
  # Right tail transform if necessary
  if(log) val <- log(val)
  return(val)
}

#' Cumulative distribution function of empirical distribution
#'
#' @param q Vector of quantiles
#' @param xs Vector of possible values; if continuous, when sorted, each value of
#' xs represents the lower limit of the interval to which the corresponding 
#' probability applies, while the next value of \code{xs} represents the upper 
#' bound; if \code{NULL}, generate from indices of \code{probs}
#' @param probs Vector of probabilities associated with the corresponding value 
#' in \code{xs}
#' @param discrete Logical; default FALSE if \code{xs} is continuous; if TRUE, 
#' indicates that \code{xs} is discrete and overrides 
#' \code{continuity_correction}
#' @param continuity_correction String denoting type of correction for 
#' continuous variables; default is "uniform" for uniformly distributed 
#' correction within intervals, but this is overridden if \code{discrete} = 
#' TRUE; if \code{NULL} is inputted, treats xs as discrete
#' @param warn_sum_probs Logical; if true, print warning if sum of probabilities 
#' is not 1; can be silenced if truncated distribution is used as a conditional 
#' distribution
#' @param max_x Maximum value of \code{xs} for continuous variables as upper bound for 
#' the interval lower bounded by the last value of \code{xs}
#' @param max_ex Remaining expected value conditional on exceeding the maximum 
#' value of \code{xs} for continuous variables
#' @param lower.tail logical; if TRUE (default), probabilities are P(X<=x), 
#' otherwise P(X>x)
#' @param log.p Logical; if TRUE, probabilities p are given as log(p)
#' 
#' @return Cumulative distribution function evaluated at q
#' 
#' @export
pempirical <- function(q, xs = NULL, probs, discrete = FALSE, 
                       continuity_correction = "uniform", max_x = NULL, 
                       max_ex = NULL, warn_sum_probs = TRUE, lower.tail = TRUE, 
                       log.p = FALSE, eps = 1e-6) {
  # Sanity checks for empirical distribution
  check_empirical(xs, probs, warn_sum_probs, eps)
  
  # Set xs if NULL to indices of probs and flag that indices will be returned, 
  # else sort if needed
  return_indices <- FALSE
  if(is.null(xs)) {
    return_indices <- TRUE
    xs <- seq(length(probs))
  } else if(is.unsorted(xs)) {
    indices <- base::order(xs)
    xs <- xs[indices]
    probs <- probs[indices]
  }
  
  # Calculate cumulative mass function (CMF)
  cmf <- c(0, cumsum(probs))
  
  # If maximum value of xs is given, append it for correction
  if (!is.null(max_x)) {
    xs <- c(xs, max_x)
  } else if (!is.null(max_ex)) {
    # Otherwise if remaining expected value conditional on exceeding maximum 
    # value of xs is given, append doubled amount for uniform correction
    xs <- c(xs, max(xs) + 2*max_ex)
  } else {
    # Otherwise append max value
    xs <- c(xs, max(xs))
  }
  
  # Find index of largest probability for which p >= CMF and set lower bound for 
  # quantile
  lb_index <- sapply(q, function(x) which.max(xs[x >= xs]))
  val <- cmf[lb_index]
  
  # Apply uniform correction
  if (!is.null(continuity_correction)) {
    if (tolower(continuity_correction) == "uniform") {
      # Get difference in probability
      v_correction <- ifelse(lb_index < length(cmf),
                             (q - xs[lb_index]) / (xs[lb_index + 1] - xs[lb_index]),
                             0)
      
      # If xs was not given and indices should be returned, return indices and 
      # generated correction
      if(return_indices) {
        return(list(indices = val, correction = v_correction))
      }
      
      # Otherwise update quantile with correction
      val <- val + (cmf[pmin(lb_index + 1, length(probs))] - cmf[lb_index]) * v_correction
    } else stop(paste0(continuity_correction, " continuity correction is not recognized"))
  }
  
  # Right tail and log transform if necessary
  if(lower.tail == FALSE) val <- 1 - val
  if(log.p) val <- log(val)
  return(val)
}

#' Quantile function of empirical distribution
#'
#' @param p Vector of probabilities
#' @param xs Vector of possible values; if continuous, when sorted, each value of
#' xs represents the lower limit of the interval to which the corresponding 
#' probability applies, while the next value of \code{xs} represents the upper 
#' bound; if \code{NULL}, generate from indices of \code{probs}
#' @param probs Vector of probabilities associated with the corresponding value 
#' in \code{xs}
#' @param discrete Logical; default FALSE if \code{xs} is continuous; if TRUE, 
#' indicates that \code{xs} is discrete and overrides 
#' \code{continuity_correction}
#' @param continuity_correction String denoting type of correction for 
#' continuous variables; default is "uniform" for uniformly distributed 
#' correction within intervals, but this is overridden if \code{discrete} = 
#' TRUE; if \code{NULL} is inputted, treats xs as discrete
#' @param warn_sum_probs Logical; if true, print warning if sum of probabilities 
#' is not 1; can be silenced if truncated distribution is used as a conditional 
#' distribution
#' @param max_x Maximum value of \code{xs} for continuous variables as upper bound for 
#' the interval lower bounded by the last value of \code{xs}
#' @param max_ex Remaining expected value conditional on exceeding the maximum 
#' value of \code{xs} for continuous variables
#' @param lower.tail logical; if TRUE (default), probabilities are P(X<=x), 
#' otherwise P(X>x)
#' @param log.p Logical; if TRUE, probabilities p are given as log(p)
#' 
#' @return Quantile of cumulative probability p
#' 
#' @export
qempirical <- function(p, xs = NULL, probs, discrete = FALSE, 
                       continuity_correction = "uniform", max_x = NULL, 
                       max_ex = NULL, warn_sum_probs = TRUE, lower.tail = TRUE, 
                       log.p = FALSE, eps = 1e-6) {
  # Sanity checks for empirical distribution
  check_empirical(xs, probs, warn_sum_probs, eps)
  
  # Right tail and log transform if necessary
  if(log.p) p <- exp(p)
  if(lower.tail == FALSE) p <- 1 - p
  
  # Set xs if NULL to indices of probs and flag that indices will be returned, 
  # else sort if needed
  return_indices <- FALSE
  if(is.null(xs)) {
    return_indices <- TRUE
    xs <- seq(length(probs))
  } else if(is.unsorted(xs)) {
    indices <- base::order(xs)
    xs <- xs[indices]
    probs <- probs[indices]
  }
  
  # Calculate cumulative mass function (CMF)
  cmf <- c(0, cumsum(probs))
  
  # If maximum value of xs is given, append it for correction
  if (!is.null(max_x)) {
    xs <- c(xs, max_x)
  } else if (!is.null(max_ex)) {
    # Otherwise if remaining expected value conditional on exceeding maximum 
    # value of xs is given, append doubled amount for uniform correction
    xs <- c(xs, max(xs) + 2*max_ex)
  } else {
    # Otherwise append max value
    xs <- c(xs, max(xs))
  }
  
  # Find index of largest probability for which p >= CMF and set lower bound for 
  # quantile
  lb_index <- sapply(p, function(x) which.max(cmf[x >= cmf]))
  val <- xs[lb_index]
  
  # Apply uniform correction
  if (!is.null(continuity_correction)) {
    if (tolower(continuity_correction) == "uniform") {
      # Get difference in probability
      v_correction <- ifelse(lb_index < length(cmf),
                             (p - cmf[lb_index]) / (cmf[lb_index + 1] - cmf[lb_index]),
                             0)
      
      # If xs was not given and indices should be returned, return indices and 
      # generated correction
      if(return_indices) {
        return(list(indices = val, correction = v_correction))
      }
      
      # Otherwise update quantile with correction
      val <- val + (xs[pmin(lb_index + 1, length(xs))] - xs[lb_index]) * v_correction
    } else stop(paste0(continuity_correction, " continuity correction is not recognized"))
  }
  
  # Return quantile
  return(val)
}

#' Generates n random values from an empirical distribution
#'
#' @param n Number of observations
#' @param xs Vector of possible values; if continuous, when sorted, each value of
#' xs represents the lower limit of the interval to which the corresponding 
#' probability applies, while the next value of \code{xs} represents the upper 
#' bound; if \code{NULL}, generate from indices of \code{probs}
#' @param probs Vector of probabilities associated with the corresponding value 
#' in \code{xs}
#' @param discrete Logical; default FALSE if \code{xs} is continuous; if TRUE, 
#' indicates that \code{xs} is discrete and overrides 
#' \code{continuity_correction}
#' @param continuity_correction String denoting type of correction for 
#' continuous variables; default is "uniform" for uniformly distributed 
#' correction within intervals, but this is overridden if \code{discrete} = 
#' TRUE; if \code{NULL} is inputted, treats xs as discrete
#' @param max_x Maximum value of \code{xs} for continuous variables as upper bound for 
#' the interval lower bounded by the last value of \code{xs}
#' @param max_ex Remaining expected value conditional on exceeding the maximum 
#' value of \code{xs} for continuous variables
#' @param warn_sum_probs Logical; if true, print warning if sum of probabilities 
#' is not 1; can be silenced if truncated distribution is used as a conditional 
#' distribution
#' 
#' @return Vector of randomly generated values
#' 
#' @export
rempirical <- function(n, xs = NULL, probs, discrete = FALSE, 
                       continuity_correction = "uniform", max_x = NULL, 
                       max_ex = NULL, return_prob = FALSE,
                       warn_sum_probs = TRUE, eps = 1e-6) {
  # Sanity checks for empirical distribution
  check_empirical(xs, probs, warn_sum_probs, eps)
  
  # Set xs if NULL to indices of probs and flag that indices will be returned, 
  # else sort if needed
  return_indices <- FALSE
  if(is.null(xs)) {
    return_indices <- TRUE
    xs <- seq(length(probs))
  } else if(is.unsorted(xs)) {
    index_order <- base::order(xs)
    xs <- xs[index_order]
    probs <- probs[index_order]
  }
  
  # Sample n values of xs weighted by probability distribution with replacement
  val <- sample(x = xs, size = n, prob = probs, replace=TRUE)
  
  # Correction for continuous variables if xs is not discrete and continuous 
  # correction is requested
  if(!discrete) {
    if (!is.null(continuity_correction)) {
      # If maximum value of xs is given, append it for correction
      if (!is.null(max_x)) xs <- c(xs, max_x)
      
      if (tolower(continuity_correction) == "uniform") {
        # Sample uniformly within intervals
        v_correction <- runif(n)
        
        # If xs was not given and indices should be returned, return indices and 
        # generated correction
        if(return_indices) {
          return(list(indices = val, correction = v_correction))
        }
        
        # If remaining expected value conditional on exceeding maximum value of  
        # xs is given, append doubled amount for uniform correction
        if (!is.null(max_ex)) xs <- c(xs, max(xs) + 2*max_ex)
        
        # Get interval width for each generated value and add correction scaled 
        # by interval width
        v_interval <- xs[pmin(length(xs), match(val, xs) + 1)] - val
        val <- val + v_correction * v_interval
      } else stop(paste0(continuity_correction, " continuity correction is not recognized"))
    }
  }
  return(val)
}

#' Sanity checks for empirical distribution; checks if probabilities sum to 1, 
#' checks that xs and probs are the same length
#'
#' @param xs Vector of possible values
#' @param probs Vector of probabilities associated with the corresponding value 
#' in \code{xs}
#' @param warn_sum_probs Logical; if true, print warning if sum of probabilities 
#' is not 1; can be silenced if truncated distribution is used as a conditional 
#' distribution
#' 
#' @return Vector of randomly generated values
#' 
#' @export
check_empirical <- function(xs, probs, warn_sum_probs = TRUE, eps = 1e-6) {
  # Check if probabilities sum to 1, otherwise output warning
  if(warn_sum_probs) {
    if(abs(sum(probs) - 1 ) > eps) {
      print("Probabilities do not sum to 1")
    }
  }
  
  # Check that xs and probs are the same length
  if(length(xs) == 0) {
    xs <- seq(length(probs))
  } else if(length(xs) != length(probs)) stop("xs and probs have different lengths")
}

