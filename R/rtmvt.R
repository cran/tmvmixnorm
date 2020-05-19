#' Random number generation for truncated multivariate Student's t distribution subject to linear inequality constraints
#'
#' \code{rtmvt} simulates truncated multivariate (p-dimensional) Student's t distribution subject to linear inequality constraints. The constraints should be written as a matrix (\code{D}) with \code{lower} and \code{upper} as the lower and upper bounds for those constraints respectively. Note that \code{D} can be non-full rank, which generalizes many traditional methods.
#'
#' @param n number of random samples desired (sample size).
#' @param Mean location vector of the multivariate Student's t distribution.
#' @param Sigma positive definite dispersion matrix of the multivariate t distribution.
#' @param nu degrees of freedom for Student-t distribution.
#' @param D matrix or vector of coefficients of linear inequality constraints.
#' @param lower lower bound vector for truncation.
#' @param upper upper bound vector for truncation.
#' @param int initial value vector for Gibbs sampler (satisfying truncation), if \code{NULL} then determine automatically.
#' @param burn burn-in iterations discarded (default as \code{10}).
#' @param thin thinning lag (default as \code{1}).
#'
#' @return \code{rtmvt} returns a (\code{n*p}) matrix (or vector when \code{n=1}) containing random numbers which follows truncated multivariate Student-t distribution.
#'
#' @examples
#' # Example for full rank
#' d <- 3
#' rho <- 0.5
#' nu <- 10
#' Sigma <- matrix(0, nrow=d, ncol=d)
#' Sigma <- rho^abs(row(Sigma) - col(Sigma))
#'
#' D1 <- diag(1,d) # Full rank
#'
#' set.seed(1203)
#' ans.t <- rtmvt(n=1000, Mean=1:d, Sigma, nu=nu, D=D1, lower=rep(-1,d), upper=rep(1,d),
#' burn=50, thin=0)
#'
#' apply(ans.t, 2, summary)
#'
rtmvt <- function(n, Mean, Sigma, nu, D, lower, upper, int=NULL, burn=10, thin=1){

  if ( any(lower >= upper)) stop("lower bound must be smaller than upper bound\n")
  bound.check <- 0

  if (!is.null(int)) {
    inits_test <- D%*%int
    lower.log <- inits_test >= lower + 1e-8  # small tol for get away from bound
    upper.log <- inits_test <= upper - 1e-8  # small tol for get away from bound
    bound.check <- prod(lower.log*upper.log)

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

  x <- solve(Sigma.chol)%*%(int-Mean) # initial value for the transformed tmvt
  # int is the initial value for the original problem

  keep.t <- matrix(NA, ncol=p, nrow=(thin+1)*n+burn)

  for (i in 1:((thin+1)*n+burn)){
    u <- stats::rchisq(1,df=nu) # sample from chisq(nu)
    denom <- sqrt(u/nu)
    lw <- c(a*denom)
    up <- c(b*denom)
    z0 <- c(x*denom)  # initial value for the tmvn for this step
    z <- c( rtmvn(n=1,Mean=rep(0,p),Sigma=diag(1,p),D=R,lower=lw,upper=up,int=z0, burn=0))
    # sample from standard tmvn

    x <- z/denom  # sample from standard tmvt

    ############################################################
    # for the original results
    w <- Sigma.chol%*%x+Mean
    keep.t[i,] <- w
  }

  final.ind <- 1:((thin+1)*n+burn)
  final.ind <- final.ind[(burn+1):length(final.ind)]
  final.ind <- seq(1,length(final.ind),by=thin+1) + thin + burn

  if (n == 1) {result <- c(keep.t[final.ind,])  } else{ result <- keep.t[final.ind,]}
  return( result )
}



