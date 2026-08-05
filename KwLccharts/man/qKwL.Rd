\name{qKwL}
\alias{qKwL}
\title{Quantile function of the KwL distribution}
\description{Closed-form quantile function of the Kumaraswamy-Lindley distribution, obtained by inverting the Kumaraswamy-G cdf and then the Lindley cdf via the Lambert W function (see \code{\link{QLindley}}).}
\usage{qKwL(u, theta, a, b)}
\arguments{
  \item{u}{Numeric vector of probabilities in (0,1).}
  \item{theta}{Lindley parameter, \eqn{\theta > 0}.}
  \item{a}{Kumaraswamy shape parameter, \eqn{a > 0}.}
  \item{b}{Kumaraswamy shape parameter, \eqn{b > 0}.}
}
\value{Numeric vector of quantiles.}
\examples{
qKwL(0.5, theta = 1, a = 1.5, b = 1.5)
pKwL(qKwL(0.5, 1, 1.5, 1.5), 1, 1.5, 1.5)
}
