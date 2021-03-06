% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unif_acc.R
\name{unif_acc}
\alias{unif_acc}
\title{Acceptance rate of uniform rejection sampling}
\usage{
unif_acc(a, b)
}
\arguments{
\item{a}{lower bound for truncation.}

\item{b}{upper bound for truncation.}
}
\description{
\code{unif_acc} calculates the acceptance rate of uniform rejection sampling for the truncation interval (a,b).
}
\examples{
set.seed(1203)
unif_acc(1,2)

}
