#' Half-normal rejection sampling
#'
#' \code{halfnorm_rej} is used for half-normal rejection sampling.
#'
#' @param a lower bound
#' @param b upper bound
#'
#' @return \code{halfnorm_rej} returns a list
#' \code{x}: sampled value; and
#' \code{acc}: total number of draw used.
#'
#' @examples
#' set.seed(1)
#' halfnorm_rej(a=1, b=Inf)
#'
halfnorm_rej <- function(a,b){
  acc <- 0
  repeat{
    x <- stats::rnorm(1)
    acc <- acc+1 # accept x at acc-th draw
    if (abs(x)>=a & abs(x)<=b)
      return(list(x=abs(x),acc=acc))
  }
}
