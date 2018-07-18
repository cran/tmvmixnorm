#' Normal rejection sampling
#'
#' \code{norm_rej()} is used for normal rejection sampling.
#'
#' @param a lower bound
#' @param b upper bound
#'
#' @return \code{norm_rej} returns a list
#' \code{x}: sampled value; and
#' \code{acc}: total number of draw used.
#'
#' @examples
#' set.seed(1)
#' norm_rej(a=1, b=Inf)
#'
norm_rej <- function(a,b=Inf){
  acc <- 0
  repeat{
    x <- stats::rnorm(1)
    acc <- acc+1 # accept x at acc-th draw
    if (x>=a&x<=b)
      return(list(x=x,acc=acc))
  }
}
