#' Translated-exponential rejection sampling
#'
#' \code{exp_rej} is used for translated-exponential rejection sampling.
#'
#' @param a lower bound
#' @param b upper bound
#' @param lam lambda for translated-exponential only
#'
#' @return \code{exp_rej} returns a list
#' \code{x}: sampled value; and
#' \code{acc}: total number of draw used.
#'
#' @examples
#' set.seed(1)
#' exp_rej(a=1, b=Inf)
#'
exp_rej <- function(a,b=Inf,lam='default'){

  if (lam=='opt'){lambda <- ( a+sqrt(a^2+4) )/2} else {lambda=a}
  acc <- 0
  # rho is continuous at lambda=a, so don't need to change the expression of rho
  repeat{
    x <- stats::rweibull(1,shape=1,scale=1/lambda)+a
    u <- stats::runif(1)
    rho <- exp(-(x-lambda)^2/2)
    acc <- acc+1 # accept x at acc-th draw
    if (u <= rho & x<b)
      return(list(x=x,acc=acc))
  }
}
