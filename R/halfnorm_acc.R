#' Acceptance rate of half-normal rejection sampling
#'
#' \code{halfnorm_acc} calculates the acceptance rate of half-normal rejection sampling for the truncation interval (a,b).
#'
#' @param a lower bound for truncation.
#' @param b upper bound for truncation.
#'
#' @examples
#' set.seed(1203)
#' halfnorm_acc(1,2)
#'
halfnorm_acc <- function(a,b){
  rate <- 2*(pnorm(b)-pnorm(a))
  return(rate)
}


