\name{fuzzify}
\alias{fuzzify}
\title{Fuzzification of Crisp Readings}
\description{
Represent crisp readings that carry imprecision as fuzzy sets (Section 3 of
Ross, 2000), e.g. a voltmeter reading that is accurate to +/- 1 percent
(Figure 1.6).
}
\usage{
fuzzify(value, spread, relative = FALSE,
        shape = c("triangular", "trapezoidal", "gaussian"),
        core = 0, name = NULL)
}
\arguments{
  \item{value}{numeric vector of crisp readings.}
  \item{spread}{half-width of the support around each reading (recycled).}
  \item{relative}{if \code{TRUE}, \code{spread} and \code{core} are fractions
    of \code{abs(value)} (e.g. \code{0.01} for +/- 1 percent).}
  \item{shape}{shape of the membership function.}
  \item{core}{half-width of the plateau of full membership
    (trapezoidal shape).}
  \item{name}{optional labels.}
}
\details{
The triangular shape has support \code{[value - spread, value + spread]} and
peak at \code{value}; the trapezoidal shape adds a core of half-width
\code{core}; the Gaussian shape uses standard deviation \code{spread / 3} on
the same support.
}
\value{
A \code{"fuzzyset"} for a single reading, otherwise a list of them.
}
\seealso{\code{\link{agreement}}}
\examples{
## a 10 V reading with +/- 1 percent imprecision (Figure 1.6)
R <- fuzzify(10, spread = 0.01, relative = TRUE)
fuzzy_support(R)
fuzzify(c(3.1, 4.7), spread = 0.5, shape = "trapezoidal", core = 0.1)
}
\keyword{math}
