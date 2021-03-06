% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normal_acc.R
\name{norm_acc}
\alias{norm_acc}
\title{Acceptance rate of normal rejection sampling}
\usage{
norm_acc(a, b)
}
\arguments{
\item{a}{lower bound for truncation.}

\item{b}{upper bound for truncation.}
}
\description{
\code{norm_acc} calculates the acceptance rate of normal rejection sampling for the truncation interval (a,b).
}
\examples{
set.seed(1203)
norm_acc(1,2)

}
