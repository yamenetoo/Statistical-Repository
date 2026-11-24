\name{fuzzy_union}
\alias{fuzzy_union}
\alias{fuzzy_intersection}
\alias{fuzzy_complement}
\alias{agreement}
\title{Union, Intersection, Complement and Agreement of Fuzzy Sets}
\description{
Standard (max / min) fuzzy set operations, computed exactly for
piecewise-linear membership functions, and the degree of agreement between a
fuzzy set and a crisp or fuzzy reading.
}
\usage{
fuzzy_union(..., name = NULL)
fuzzy_intersection(..., name = NULL)
fuzzy_complement(A, universe = range(A$x), name = NULL)
agreement(A, reading)
}
\arguments{
  \item{\dots}{fuzzy sets, or a single list of fuzzy sets.}
  \item{name}{optional label of the result.}
  \item{A}{a \code{"fuzzyset"} object.}
  \item{universe}{\code{c(lower, upper)}: the universe on which the complement
    \eqn{1 - \mu_A(x)} is taken.}
  \item{reading}{a crisp reading (numeric) or a fuzzy reading
    (\code{"fuzzyset"}).}
}
\details{
The union \eqn{\mu_{A \cup B}(x) = \max(\mu_A(x), \mu_B(x))} and the
intersection \eqn{\mu_{A \cap B}(x) = \min(\mu_A(x), \mu_B(x))} of two
piecewise-linear functions are piecewise linear with knots at the knots of
both sets and at the points where the two functions cross; these points are
added, so the result is exact.

\code{agreement} follows Figure 1.7 of Ross (2000): for a crisp reading it is
the membership of the reading in \code{A}; for a fuzzy reading it is the
height of the intersection of the two sets.
}
\value{
A \code{"fuzzyset"} object, or for \code{agreement} a number in [0, 1].
}
\seealso{\code{\link{defuzzify}}, \code{\link{fuzzify}}}
\examples{
out <- ross2000_sets("output")
C <- fuzzy_union(out)          # Figure 1.15c
C
fuzzy_intersection(out$C1, out$C2)

## convex sets have a convex intersection (Figure 1.5)
is_convex(fuzzy_intersection(fs_triangular(0, 2, 5), fs_trapezoidal(2, 4, 6, 8)))

v <- ross2000_sets("voltage")
agreement(v$low, 5.8)                         # crisp reading
agreement(v$medium, fuzzify(7.2, spread = 1)) # fuzzy reading
}
\keyword{math}
