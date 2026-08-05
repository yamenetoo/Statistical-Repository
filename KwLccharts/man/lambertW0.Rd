\name{lambertW0}
\alias{lambertW0}
\title{Lambert W function, principal branch}
\description{Numerically evaluates the principal branch \eqn{W_0(x)} of the Lambert W function via Halley's iteration, with no external dependencies.}
\usage{lambertW0(x)}
\arguments{\item{x}{Numeric vector, \eqn{x \ge -1/e}.}}
\value{Numeric vector of \eqn{W_0(x)} values.}
\examples{
lambertW0(1)
lambertW0(exp(1))
}
