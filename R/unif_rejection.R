#' Uniform rejection sampling
#'
#' \code{unif_rej} is used for uniform rejection sampling.
#'
#' @param a lower bound
#' @param b upper bound
#'
#' @return \code{unif_rej} returns a list
#' \code{x}: sampled value; and
#' \code{acc}: total number of draw used.
#'
#' @examples
#' set.seed(1)
#' unif_rej(a=1, b=2)
#'
unif_rej <- function(a,b){
  acc <- 0
  repeat{
    x <- stats::runif(1,a,b)
    u <- stats::runif(1)
    # different cases for the ratio
    if (0>=a & 0<=b){rho <- exp(-x^2/2)}
    if (a>0){rho <- exp(-(x^2-a^2)/2)}
    if (b<0){rho <- exp(-(x^2-b^2)/2)}
    acc <- acc+1 # accept x at acc-th draw
    if (u <= rho)
      return(list(x=x,acc=acc))
  }
}

