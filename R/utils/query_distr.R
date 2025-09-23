#' General function to query a distribution's density, cumulative distribution 
#' function, quantile function, or random generation function following the 
#' format in the R stats package
#'
#' \code{query_distr} is a flexible method to generate values from an inputted 
#' distribution
#'
#' @param target Type of output required, consistent with the R stats package: 
#' \code{d} for density, \code{p} for cumulative probability, \code{q} for 
#' quantile, and \code{r} for randomly generated value
#' @param x For \code{target} = \code{d} or \code{p}, the vector of quantiles 
#' at which to evaluate; for \code{target} = \code{q}, the vector of 
#' probabilities; for \code{target} = \code{r}, the number of observations to 
#' generate
#' @param distr String with the distribution name
#' @param params List of distribution parameters named correspondingly to the 
#' distribution's family of functions
#' @param ... List of optional arguments to the target distribution function
#' 
#' @return Target value from the distribution
#' 
#' @export
query_distr <- function(target, x, distr, params, ...) {
  # Get distribution name for function
  if(tolower(distr) == "exponential") {
    distr_name <- "exp"
  } else if(tolower(distr) == "binomial") {
      distr_name <- "binom"
  } else if(tolower(distr) == "uniform") {
    distr_name <- "unif"
  } else {
    distr_name <- tolower(distr)
  }
  
  # Call function based on target and distribution name if function exists
  function_name <- paste0(target, distr_name)
  if(is.function(match.fun(function_name))) {
    # If cure_max parameter exists, save cure parameter
    if ("cure_max" %in% names(params)) {
      cure_max <- params$cure_max
      params$cure_max <- NULL
      
      # If target is quantile, adjust probability for unscaled model
      if (target == "q") {
        x <- x / cure_max
        x[x > 1] <- NA
      }
    }
    
    # Calculate distribution results
    val <- do.call(match.fun(function_name), c(list(x), params, ...))
    
    # If cure_max parameter exists, adjust for cure model
    if (exists("cure_max")) {
      val <- cure_adjust(target, val, cure_max)
    }
    return(val)
  } else {
    stop(paste0("Function ", function_name, " not found"))
  }
}


# Function to adjust outputs for cure model (where CDF_final = cure_max*CDF_initial)
cure_adjust <- function(target,
                        val,
                        cure_max) {
  # If sampling from distribution, randomly select values to keep or set NA
  if (target == "r") {
    v_keep <- rbinom(length(val), 1, cure_max)
    val[v_keep == 0] <- NA
  } else if (target %in% c("p", "d")) {
    # If calculating PDF or CDF, scale by cure_max
    val <- val * cure_max
  } else if (!target %in% "q") {
    # Adjustment assumed to already be made for quantile
    stop("target must be 'r', 'q', 'p', or 'd'")
  }
  return(val)
}