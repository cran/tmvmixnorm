#' Rejection sampling of standardized truncated univariate normal distribution
#'
#' \code{imp} contains a general function for rejection sampling of standardized truncated univariate normal distribution in (a,b).
#'
#' @param a lower bound for truncation.
#' @param b upper bound for truncation.
#'
#' @return \code{imp} returns a list
#' \code{x}: sampled value; and
#' \code{acc}: total number of draw used.
#'
#' @examples
#' imp(1,Inf) # Case 1: [a,infty)
#' imp(-1,1) # Case 2: 0 in [a,b], a<0<b
#' imp(1,2) # Case 3: [a,b], a>=0
#'
imp <- function(a,b){
  lowerb1 <- function(a){
    lowerb1 <- sqrt(pi/2)*exp(a^2/2)+a
    return(lowerb1)
  }

  # lower bound of b for normal vs unif in 0 in [a,b]
  lowerb2 <- function(a){
    lowerb2 <- sqrt(2*pi)+a
    return(lowerb2)
  }


  # lower bound of b for exp vs unif in [a,b]>=0
  lowerb <- function(a){
    lambda <- ( a+sqrt(a^2+4) )/2
    lowerb <- a+ exp(1/2)/lambda * exp(  ( a^2-a*sqrt(a^2+4) )/4  )
    return(lowerb)
  }


  # Case 1: [a,infty)
  imp_case1 <- function(a,b=Inf){
    if (a<0) {
      samp <- norm_rej(a=a,b=b)
    } else {
      if (a<0.25696) {
        samp <- halfnorm_rej(a=a,b=b)
      } else {
        samp <- exp_rej(a=a,b=b,lam='opt')
      }
    }
    return(samp)
  }

  # Case 2: 0 in [a,b], a<0<b
  imp_case2 <- function(a,b){
    if (b>lowerb2(a)) {
      samp <- norm_rej(a=a,b=b)
    } else {
      samp <- unif_rej(a=a,b=b)}
    return(samp)
  }

  # Case 3: [a,b], a>=0
  imp_case3 <- function(a,b){
    if (a<=0.25696) {
      blower1 <- lowerb1(a)
      if (b<= blower1) {
        samp <- unif_rej(a=a,b=b)
      } else {
        samp <- halfnorm_rej(a=a,b=b)
      }
    } else {
      blower2 <- lowerb(a)
      if (b<=blower2) {
        samp <- unif_rej(a=a,b=b)
      } else {
        samp <- exp_rej(a=a,b=b,lam='opt')
      }
    }
    return(samp)
  }

  # Case 4: (-infty,b] (symmetric to Case 1)
  imp_case4 <- function(a,b){
    temp <- imp_case1(a=-b,b=-a)
    x <- -temp$x
    samp <- list(x=x,acc=temp$acc)
    return(samp)
  }

  # Case 5: [a,b], b<=0 (symmetric to Case 3)
  imp_case5 <- function(a,b){
    temp <- imp_case3(a=-b,b=-a)
    x <- -temp$x
    samp <- list(x=x,acc=temp$acc)
    return(samp)
  }

  if (a==-Inf | b==Inf) {
    if (b==Inf) {
      samp <- imp_case1(a=a,b=b)} else {
        samp <- imp_case4(a=a,b=b) }
  } else {
    if (a<0 & b>0) {samp <- imp_case2(a=a,b=b)}
    if (a>=0) {samp <- imp_case3(a=a,b=b)}
    if (b<=0) {samp <- imp_case5(a=a,b=b)}
  }
  return(samp)
}
