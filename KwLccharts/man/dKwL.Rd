\name{dKwL}
\alias{dKwL}
\title{Density of the Kumaraswamy-Lindley (KwL) distribution}
\description{Density of the three-parameter Kumaraswamy-Lindley distribution, obtained from the Lindley baseline via the Kumaraswamy-G construction of Cordeiro and de Castro (2011). Reduces to the Lindley density when \code{a = b = 1}.}
\usage{dKwL(x, theta, a, b, log = FALSE)}
\arguments{
  \item{x}{Non-negative numeric vector of quantiles.}
  \item{theta}{Lindley parameter, \eqn{\theta > 0}.}
  \item{a}{Kumaraswamy shape parameter, \eqn{a > 0}.}
  \item{b}{Kumaraswamy shape parameter, \eqn{b > 0}.}
  \item{log}{Logical; if \code{TRUE}, the log-density is returned.}
}
\value{Numeric vector of (log-)densities.}
\examples{
dKwL(1.5, theta = 1, a = 1.5, b = 1.5)
integrate(dKwL, 0, Inf, theta = 1, a = 1.5, b = 1.5)$value
}
