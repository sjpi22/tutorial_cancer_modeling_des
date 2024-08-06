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
    val <- do.call(match.fun(function_name), c(list(x), params, ...))
    return(val)
  } else {
    stop(paste0("Function ", function_name, " not found"))
  }
}