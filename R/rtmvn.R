#' Random number generation for truncated multivariate normal distribution subject to linear inequality constraints
#'
#' \code{rtmvn} simulates truncated multivariate (p-dimensional) normal distribution subject to linear inequality constraints. The constraints should be written as a matrix (\code{D}) with \code{lower} and \code{upper} as the lower and upper bounds for those constraints respectively. Note that \code{D} can be non-full rank, which generalize many traditional methods.
#'
#' @param n number of random samples desired (sample size).
#' @param Mean mean vector of the underlying multivariate normal distribution.
#' @param Sigma positive definite covariance matrix of the underlying multivariate normal distribution.
#' @param D matrix or vector of coefficients of linear inequality constraints.
#' @param lower vector of lower bounds for truncation.
#' @param upper vector of upper bounds for truncation.
#' @param int initial value vector for Gibbs sampler (satisfying truncation), if \code{NULL} then determine automatically.
#' @param burn burn-in iterations discarded (default as \code{10}).
#' @param thin thinning lag (default as \code{1}).
#'
#' @return \code{rtmvn} returns a (\code{n*p}) matrix (or vector when \code{n=1}) containing random numbers which approximately follows truncated multivariate  normal distribution.
#'
#' @examples
#' # Example for full rank with strong dependence
#' d <- 3
#' rho <- 0.9
#' Sigma <- matrix(0, nrow=d, ncol=d)
#' Sigma <- rho^abs(row(Sigma) - col(Sigma))
#'
#' D1 <- diag(1,d) # Full rank
#'
#' set.seed(1203)
#' ans.1 <- rtmvn(n=1000, Mean=1:d, Sigma, D=D1, lower=rep(-1,d), upper=rep(1,d),
#' int=rep(0,d), burn=50)
#'
#' apply(ans.1, 2, summary)
#'
#' # Example for non-full rank
#' d <- 3
#' rho <- 0.5
#' Sigma <- matrix(0, nrow=d, ncol=d)
#' Sigma <- rho^abs(row(Sigma) - col(Sigma))
#'
#' D2 <- matrix(c(1,1,1,0,1,0,1,0,1),ncol=d)
#' qr(D2)$rank # 2
#'
#' set.seed(1228)
#' ans.2 <- rtmvn(n=100, Mean=1:d, Sigma, D=D2, lower=rep(-1,d), upper=rep(1,d), burn=10)
#'
#' apply(ans.2, 2, summary)
#'
rtmvn <- function(n, Mean, Sigma, D=diag(1,length(Mean)), lower, upper, int=NULL, burn=10, thin=1){
  if (length(Mean) == 1) {
    result <- rtuvn(n=n, mean=Mean, sd=c(Sigma), lower=lower, upper=upper)
  } else{

  if ( any(lower >= upper)) stop("lower bound must be smaller than upper bound\n")

  bound.check <- 0

  if (!is.null(int)) {
    inits_test <- D%*%int
    lower.log <- inits_test >= lower + 1e-8  # small tol for get away from bound
    upper.log <- inits_test <= upper - 1e-8  # small tol for get away from bound
    bound.check <- prod(lower.log*upper.log)
    if(bound.check == 0) cat("initial is outside or too close from boundary, will be auto-correlated!\n")
  } else if (bound.check == 0) {

    D.inv <- MASS::ginv(D)
    int <- D.inv%*%(lower + upper)/2
  }


  if( any (c(burn,thin,n) %% 1 != 0))  stop("burn, thin and n must be integer\n")
  if ( any(c(burn, thin, n -1) < 0) ) stop("burn, thin must be  non-negative interger, n must be positive integer\n")


  if(is.vector(D)==TRUE){
    Rtilde <- t(as.matrix(D))
    lower <- as.vector(lower)
    upper <- as.vector(upper)
  } else {
    Rtilde <- D
  }

  a <- lower - Rtilde%*%Mean
  b <- upper - Rtilde%*%Mean
  Sigma.chol <- t(chol(Sigma))
  R <- Rtilde%*%Sigma.chol

  p <- ncol(R) # number of parameters, i.e. length of beta vector
  #   m <- nrow(R) # number of constraints

  z <- solve(Sigma.chol)%*%(int-Mean)  # int is the initial value for the original problem

  keep.x <- matrix(NA, ncol=p, nrow=(thin+1)*n+burn)

  for (i in 1:((thin+1)*n+burn)){
    for (j in 1:p){  # p is the number of the degree of the multinormal
      rj <- as.vector(R[,j])   # the jth column of R
      Rj <- as.matrix(R[,-j])  # m by p-1 matrix by removing the jth column of R
      zj <- as.vector(z[-j])   # p-1 vector by removing the jth element
      a.temp <- a - Rj%*%zj
      b.temp <- b - Rj%*%zj

      pos <- rj>0   # if r_jk>0
      neg <- rj<0   # if r_jk<0

      if(sum(pos)==0){
        lower.pos <- -Inf
        upper.pos <- Inf
      } else {
        lower.pos <- max(a.temp[pos]/rj[pos])  # when r_jk>0, a associated with lower bound
        upper.pos <- min(b.temp[pos]/rj[pos])  #              b                 upper
      }

      if(sum(neg)==0){
        upper.neg <- Inf
        lower.neg <- -Inf
      } else {
        upper.neg <- min(a.temp[neg]/rj[neg])  # when r_jk<0, a                 upper
        lower.neg <- max(b.temp[neg]/rj[neg])  #              b                 lower
      }

      lower.j <- max(lower.pos,lower.neg)   # lower for z_j
      upper.j <- min(upper.pos,upper.neg)   # upper

      z[j] <- rtuvn(lower=lower.j,upper=upper.j)
    }
    ############################################################

    # for the original results
    x <- Sigma.chol%*%z+Mean
    keep.x[i,] <- x
  }

  final.ind <- 1:((thin+1)*n+burn)
  final.ind <- final.ind[(burn+1):length(final.ind)]
  final.ind <- seq(1,length(final.ind),by=thin+1) + thin + burn

  if (n == 1) {result <- c(keep.x[final.ind,])  } else{ result <- keep.x[final.ind,]}
  }

  return( result )
}
