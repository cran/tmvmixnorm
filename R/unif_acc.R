#' Acceptance rate of uniform rejection sampling
#'
#' \code{unif_acc} calculate the acceptance rate of uniform rejection sampling for the truncation interval (a,b).
#'
#' @param a lower bound for truncation.
#' @param b upper bound for truncation.
#'
#' @examples
#' set.seed(1203)
#' unif_acc(1,2)
#'
unif_acc <- function(a,b){
  cons <- sqrt(2*pi)/(b-a)
  pab <- pnorm(b)-pnorm(a)
  cp <- cons*pab
  if (a<=0 & b>=0) {rate <- cp}
  if (a>0) {rate <- cp*exp(a^2/2)}
  return(rate)
}
