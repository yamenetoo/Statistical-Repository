\name{negloglik_KwL}
\alias{negloglik_KwL}
\title{Negative log-likelihood of the KwL distribution}
\description{Negative log-likelihood used internally by \code{\link{fit_KwL}}.}
\usage{negloglik_KwL(par, y)}
\arguments{
  \item{par}{Numeric vector \code{c(theta, a, b)}.}
  \item{y}{Numeric vector of observations.}
}
\value{Negative log-likelihood (a large sentinel value if \code{par} is infeasible).}
