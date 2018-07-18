#' Acceptance rate of truncated univariate normal distribution rejection sampling
#'
#' \code{imp_acc} calculates the acceptance rate of truncated univariate standardized normal distribution rejection sampling for the truncation interval (a,b).
#'
#' @param a lower bound for truncation.
#' @param b upper bound for truncation.
#'
#' @examples
#' imp_acc(1,Inf) # Case 1: [a,infty)
#' imp_acc(-1,1) # Case 2: 0 in [a,b], a<0<b
#' imp_acc(1,2) # Case 3: [a,b], a>=0
#'
imp_acc <- function(a,b){
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
  imp_acc_case1 <- function(a,b=Inf){
    if (a<0) {
      rate <- norm_acc(a=a,b=b)
    } else {
      if (a<0.25696) {
        rate <- halfnorm_acc(a=a,b=b)
      } else {
        rate <- exp_acc_opt(a=a,b=b)
      }
    }
    return(rate)
  }

  # Case 2: 0 in [a,b], a<0<b
  imp_acc_case2 <- function(a,b){
    if (b>lowerb2(a)) {
      rate <- norm_acc(a=a,b=b)
    } else {
      rate <- unif_acc(a=a,b=b)}
    return(rate)
  }

  # Case 3: [a,b], a>=0
  imp_acc_case3 <- function(a,b){
    if (a<=0.25696) {
      blower1 <- lowerb1(a)
      if (b<= blower1) {
        rate <- unif_acc(a=a,b=b)
      } else {
        rate <- halfnorm_acc(a=a,b=b)
      }
    } else {
      blower2 <- lowerb(a)
      if (b<=blower2) {
        rate <- unif_acc(a=a,b=b)
      } else {
        rate <- exp_acc_opt(a=a,b=b)
      }
    }
    return(rate)
  }

  if (a==-Inf | b==Inf) {
    if (b==Inf) {
      rate <- imp_acc_case1(a=a)} else {
        rate <- imp_acc_case1(a=-b) }
  } else {
    if (a<0 & b>0) {rate <- imp_acc_case2(a=a,b=b)}
    if (a>=0) {rate <- imp_acc_case3(a=a,b=b)}
    if (b<=0) {rate <- imp_acc_case3(a=-b,b=-a)}
  }
  return(rate)
}

