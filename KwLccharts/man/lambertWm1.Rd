\name{lambertWm1}
\alias{lambertWm1}
\title{Lambert W function, lower real branch}
\description{Numerically evaluates the lower (negative) real branch \eqn{W_{-1}(x)} of the Lambert W function on \eqn{[-1/e,0)}, via Halley's iteration.}
\usage{lambertWm1(x)}
\arguments{\item{x}{Numeric vector with \eqn{-1/e \le x < 0}.}}
\value{Numeric vector of \eqn{W_{-1}(x)} values (\code{NA} outside the valid domain).}
\examples{
lambertWm1(-0.1)
lambertWm1(-1/exp(1))
}
