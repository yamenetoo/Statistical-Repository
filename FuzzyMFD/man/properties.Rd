\name{fuzzy_core}
\alias{fuzzy_core}
\alias{fuzzy_support}
\alias{fuzzy_boundary}
\alias{fuzzy_height}
\alias{crossover_points}
\alias{fuzzy_prototype}
\alias{is_normal}
\alias{is_subnormal}
\alias{is_convex}
\alias{is_fuzzy_number}
\title{Features of a Membership Function}
\description{
The terms of Section 2.1 of Ross (2000) that describe the shape of a
membership function.
}
\usage{
fuzzy_core(A, tol = 1e-9)
fuzzy_support(A, tol = 1e-9)
fuzzy_boundary(A, tol = 1e-9)
fuzzy_height(A)
crossover_points(A, level = 0.5, tol = 1e-9)
fuzzy_prototype(A, tol = 1e-9)
is_normal(A, tol = 1e-9)
is_subnormal(A, tol = 1e-9)
is_convex(A, tol = 1e-9)
is_fuzzy_number(A, tol = 1e-9)
}
\arguments{
  \item{A}{a \code{"fuzzyset"} object.}
  \item{tol}{numerical tolerance.}
  \item{level}{membership level of the crossover points (0.5 in the chapter).}
}
\details{
\describe{
  \item{core}{elements with full membership, \eqn{\mu_A(x) = 1}.}
  \item{support}{elements with non-zero membership, \eqn{\mu_A(x) > 0}
    (returned as closed intervals).}
  \item{boundaries}{elements with partial membership, \eqn{0 < \mu_A(x) < 1}.}
  \item{height}{\eqn{\max_x \mu_A(x)}; a set with height 1 is normal, otherwise
    sub-normal.}
  \item{prototype}{the element of a set whose core is a single point,
    \eqn{norm(A)}; \code{NA} otherwise.}
  \item{convex}{Equation (3): \eqn{\mu_A(y) \ge \min[\mu_A(x), \mu_A(z)]} for all
    \eqn{x < y < z}.}
  \item{fuzzy number}{a convex normal fuzzy set on the real line.}
  \item{crossover points}{elements with \eqn{\mu_A(x) = 0.5}.}
}
For piecewise-linear membership functions all of these are computed exactly
from the knots.
}
\value{
\code{fuzzy_core}, \code{fuzzy_support} and \code{fuzzy_boundary} return a
two-column matrix of intervals (\code{lower}, \code{upper}); a single point is
an interval with equal ends. \code{fuzzy_height} and \code{fuzzy_prototype} return a
number, \code{crossover_points} a numeric vector and the \code{is_*} functions
a logical value.
}
\seealso{\code{\link{fuzzy_set}}}
\examples{
## convex normal set and non-convex normal set (Figure 1.4)
A <- fuzzy_set(c(0, 2, 5), c(0, 1, 0))
B <- fuzzy_set(c(0, 2, 3, 4, 6), c(0, 1, 0.7, 0.9, 0))
is_convex(A)
is_convex(B)
is_fuzzy_number(A)

## core, support and boundaries (Figure 1.2)
Tz <- fs_trapezoidal(1, 3, 5, 8)
fuzzy_core(Tz)
fuzzy_support(Tz)
fuzzy_boundary(Tz)
crossover_points(Tz)

## a sub-normal set
S <- fs_triangular(0, 1, 2, height = 0.6)
fuzzy_height(S)
is_subnormal(S)
fuzzy_prototype(fs_triangular(0, 1, 2))
}
\keyword{math}
