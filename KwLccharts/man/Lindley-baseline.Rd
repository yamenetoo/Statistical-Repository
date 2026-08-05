\name{Lindley-baseline}
\alias{gLindley}
\alias{GLindley}
\title{Lindley density and cumulative distribution function}
\description{
Density (\code{gLindley}) and cumulative distribution function
(\code{GLindley}) of the one-parameter Lindley distribution, the baseline
of the Kumaraswamy-Lindley family.
}
\usage{
gLindley(x, theta)
GLindley(x, theta)
}
\arguments{
  \item{x}{Non-negative numeric vector.}
  \item{theta}{Lindley parameter, \eqn{\theta > 0}.}
}
\value{
\code{gLindley} returns a numeric vector of densities;
\code{GLindley} returns a numeric vector of cumulative probabilities.
}
