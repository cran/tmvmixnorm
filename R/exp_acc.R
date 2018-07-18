#' Acceptance rate of translated-exponential rejection sampling
#'
#' \code{exp_acc_opt} calculate the acceptance rate of translated-exponential rejection sampling for the truncation interval (a,b).
#'
#' @param a lower bound for truncation.
#' @param b upper bound for truncation.
#'
#' @examples
#' set.seed(1203)
#' exp_acc_opt(1,2)
#'
exp_acc_opt <- function(a,b){
  lambda <- ( a+sqrt(a^2+4) )/2
  rate <- sqrt(2*pi)*(pnorm(b)-pnorm(a))*lambda*exp(-lambda^2/2+lambda*a)
  return(rate)
}
