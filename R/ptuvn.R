#' Distribution function of truncated univariate normal distribution
#'
#' \code{ptuvn} calculates the cumulative distribution function (cdf) of truncated univariate normal distribution.
#'
#' @param x value at which cdf is desired.
#' @param mean mean of the underlying univariate normal distribution.
#' @param sd standard deviation of the underlying univariate normal distribution.
#' @param lower lower bound for truncation.
#' @param upper upper bound for truncation.
#'
#' @return \code{ptuvn} returns the cumulative distribution function (with same dimension and type as \code{x}) of truncated univariate normal distribution.
#'
#' @examples
#' ptuvn(x= -3:3, mean=0, sd=1 ,lower= -2, upper=2)
#'
ptuvn <- function(x, mean, sd, lower, upper){
  result <- (stats::pnorm(x,mean,sd)- stats::pnorm(lower,mean,sd))/(stats::pnorm(upper,mean,sd) - stats::pnorm(lower,mean,sd))
  result[x <= lower] <- 0 ; result[x >= upper] <- 1

  result
}

