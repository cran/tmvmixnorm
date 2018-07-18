#' Random number generation for truncated univariate normal distribution
#'
#' \code{rtuvn} simulate truncated univariate normal distribution within the interval.
#'
#' @param n number of random number generating.
#' @param mean mean of the underlying univariate normal distribution.
#' @param sd standard deviation of the underlying univariate normal distribution.
#' @param lower lower bound for truncation.
#' @param upper upper bound for truncation.
#'
#' @return \code{rtuvn} returns a vector of random number follows truncated univariate normal distribution.
#'
#' @examples
#' set.seed(1203)
#' ans <- rtuvn(n=1000, mean=1, sd=2, lower=-2, upper=3)
#' summary(ans)
#'
#' # Check if the sample matches with CDF by KS test
#' ks.test(ans,"ptuvn",1,2,-2,3)
#'
rtuvn <- function(n=1, mean=0, sd=1, lower, upper){
  # transform the boundaries
  a <- (lower - mean)/sd
  b <- (upper - mean)/sd

  # generate n samples from TN(0,1;a,b)
  Z <- rep(0,n)
  for (i in 1:n){
    temp <- imp(a,b)
    Z[i] <- temp$x
  }

  # transform the data back
  samp <- sd*Z + mean

  return(samp)
}
