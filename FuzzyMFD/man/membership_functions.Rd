\name{mf_triangular}
\alias{mf_triangular}
\alias{mf_trapezoidal}
\alias{mf_gaussian}
\alias{mf_bell}
\alias{mf_sigmoid}
\alias{mf_indicator}
\title{Membership Function Shapes}
\description{
Vectorised membership functions: triangular, trapezoidal, Gaussian,
generalised bell, sigmoid and the indicator function of a crisp interval.
}
\usage{
mf_triangular(x, a, b, c, height = 1)
mf_trapezoidal(x, a, b, c, d, height = 1)
mf_gaussian(x, mean, sd)
mf_bell(x, a, b, c)
mf_sigmoid(x, a, c)
mf_indicator(x, lower, upper)
}
\arguments{
  \item{x}{numeric vector of elements of the universe.}
  \item{a, b, c, d}{shape parameters: corners of the triangle or trapezoid;
    for \code{mf_bell} the width \code{a}, slope \code{b} and centre \code{c};
    for \code{mf_sigmoid} the slope \code{a} and crossover point \code{c}.}
  \item{height}{maximum membership.}
  \item{mean, sd}{centre and spread of the Gaussian shape.}
  \item{lower, upper}{limits of the crisp interval.}
}
\details{
\code{mf_indicator} is the indicator function of Equation (1) of Ross (2000),
equal to 1 for \code{lower <= x <= upper} and 0 otherwise. The generalised
bell is \eqn{1 / (1 + |(x - c)/a|^{2b})} and the sigmoid is
\eqn{1 / (1 + \exp(-a(x - c)))}.
}
\value{A numeric vector of membership values in [0, 1].}
\seealso{\code{\link{fuzzy_set_fun}} to turn any of them into a fuzzy set.}
\examples{
x <- seq(0, 80, by = 10)
mf_triangular(x, 0, 20, 40)
mf_trapezoidal(x, 10, 20, 40, 60)
mf_indicator(c(4.99, 6), 5, 7)
curve(mf_bell(x, 10, 2, 40), 0, 80, ylab = "membership")
}
\keyword{math}
