\name{rKwL}
\alias{rKwL}
\title{Random generation from the KwL distribution}
\description{Generates variates from the Kumaraswamy-Lindley distribution by the inversion method.}
\usage{rKwL(n, theta, a, b)}
\arguments{
  \item{n}{Number of variates to generate.}
  \item{theta}{Lindley parameter, \eqn{\theta > 0}.}
  \item{a}{Kumaraswamy shape parameter, \eqn{a > 0}.}
  \item{b}{Kumaraswamy shape parameter, \eqn{b > 0}.}
}
\value{Numeric vector of length \code{n}.}
