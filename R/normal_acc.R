#' Acceptance rate of normal rejection sampling
#'
#' \code{norm_acc} calculates the acceptance rate of normal rejection sampling for the truncation interval (a,b).
#'
#' @param a lower bound for truncation.
#' @param b upper bound for truncation.
#'
#' @examples
#' set.seed(1203)
#' norm_acc(1,2)
#'
norm_acc <- function(a,b){
  rate <- pnorm(b)-pnorm(a)
  return(rate)
}

