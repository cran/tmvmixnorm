% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/halfnorm_rejction.R
\name{halfnorm_rej}
\alias{halfnorm_rej}
\title{Half-normal rejection sampling}
\usage{
halfnorm_rej(a, b)
}
\arguments{
\item{a}{lower bound}

\item{b}{upper bound}
}
\value{
\code{halfnorm_rej} returns a list
\code{x}: sampled value; and
\code{acc}: total number of draw used.
}
\description{
\code{halfnorm_rej} is used for half-normal rejection sampling.
}
\examples{
set.seed(1)
halfnorm_rej(a=1, b=Inf)

}
