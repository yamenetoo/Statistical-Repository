\name{QLindley}
\alias{QLindley}
\title{Lindley quantile function via the Lambert W function}
\description{Closed-form inverse of \code{\link{GLindley}}, obtained as the p = 0 special case of the Lindley-Geometric quantile function of Al-Lami, Jabbari Khamnei and Heydari (2025).}
\usage{QLindley(v, theta)}
\arguments{
  \item{v}{Numeric vector of probabilities in (0,1).}
  \item{theta}{Lindley parameter, \eqn{\theta > 0}.}
}
\value{Numeric vector of Lindley quantiles.}
