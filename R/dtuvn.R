#' Density function of truncated univariate normal distribution
#'
#' \code{dtuvn} calculates the probability density function (pdf) of truncated univariate normal distribution.
#'
#' @param x value at which density is desired.
#' @param mean mean of the underlying univariate normal distribution.
#' @param sd standard deviation of the underlying univariate normal distribution.
#' @param lower lower bound for truncation.
#' @param upper upper bound for truncation.
#'
#' @return \code{dtuvn} returns the density (with same dimension and type as \code{x}) of truncated univariate normal distribution.
#'
#' @examples
#' dtuvn(x= -3:3, mean=0, sd=1 ,lower= -2, upper=2)
#'
dtuvn <- function(x,mean,sd,lower,upper){
  value <- stats::dnorm(x,mean,sd)/(stats::pnorm(upper,mean,sd) - stats::pnorm(lower,mean,sd))
  return(ifelse(x<lower | x>upper,0,value))
}
