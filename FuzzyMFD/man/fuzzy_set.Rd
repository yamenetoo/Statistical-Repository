\name{fuzzy_set}
\alias{fuzzy_set}
\alias{fuzzy_set_fun}
\alias{fs_triangular}
\alias{fs_trapezoidal}
\alias{fs_rectangular}
\alias{fs_singleton}
\alias{fs_gaussian}
\alias{membership}
\alias{print.fuzzyset}
\alias{summary.fuzzyset}
\alias{print.summary.fuzzyset}
\title{Fuzzy Sets on the Real Line}
\description{
Create fuzzy sets stored as piecewise-linear membership functions and
evaluate their membership values.
}
\usage{
fuzzy_set(x, mu, name = NULL)
fuzzy_set_fun(fun, lower, upper, n = 1001L, name = NULL)
fs_triangular(a, b, c, height = 1, name = NULL)
fs_trapezoidal(a, b, c, d, height = 1, name = NULL)
fs_rectangular(a, b, name = NULL)
fs_singleton(x0, height = 1, name = NULL)
fs_gaussian(mean, sd, lower = mean - 4 * sd, upper = mean + 4 * sd,
            n = 401L, name = NULL)
membership(A, z)

\method{print}{fuzzyset}(x, digits = 4, \dots)
\method{summary}{fuzzyset}(object, \dots)
\method{print}{summary.fuzzyset}(x, digits = 4, \dots)
}
\arguments{
  \item{x}{for \code{fuzzy_set}: distinct knots; for the methods: the object.}
  \item{mu}{membership values in [0, 1] at the knots.}
  \item{name}{optional label of the set.}
  \item{fun}{a vectorised membership function, evaluated on a grid.}
  \item{lower, upper}{range of the grid.}
  \item{n}{number of grid points.}
  \item{a, b, c, d}{corner points of the triangle (\code{a <= b <= c}) or
    trapezoid (\code{a <= b <= c <= d}); for \code{fs_rectangular} the crisp
    interval \code{[a, b]}.}
  \item{height}{maximum membership (1 for a normal set).}
  \item{x0}{location of a singleton.}
  \item{mean, sd}{centre and spread of a Gaussian membership function.}
  \item{A}{a \code{"fuzzyset"} object.}
  \item{z}{points at which to evaluate the membership.}
  \item{digits}{significant digits to print.}
  \item{object}{a \code{"fuzzyset"} object.}
  \item{\dots}{unused.}
}
\details{
The membership is interpolated linearly between the knots and is zero
outside \code{[min(x), max(x)]}; a set whose first or last knot has positive
membership therefore jumps to zero there (as a crisp interval does).
Triangles, trapezoids and crisp intervals are represented exactly; smooth
shapes are represented on a grid of \code{n} points.

\code{summary} reports the features of Section 2 of Ross (2000): height,
normality, convexity, core, support, boundaries, crossover points, prototype
and centroid.
}
\value{
An object of class \code{"fuzzyset"}: a list with components \code{x},
\code{mu} and \code{name}. \code{membership} returns a numeric vector.
}
\references{
Ross, T. J. (2000). Membership functions, fuzzification and defuzzification.
In \emph{Fuzzy Systems in Medicine}, 49--78. Springer-Verlag.
}
\seealso{\code{\link{mf_triangular}}, \code{\link{fuzzy_core}},
  \code{\link{plot.fuzzyset}}}
\examples{
## the fuzzy set H of heights "near 6 feet" (Figure 1.1b)
H <- fs_triangular(4.5, 6, 7.5, name = "near 6 feet")
H
membership(H, c(5, 5.917, 6, 7))
summary(H)

## the crisp set A = [5, 7] (Figure 1.1a)
A <- fs_rectangular(5, 7, name = "A")
membership(A, c(4.99, 6))

## a smooth shape on a grid
G <- fs_gaussian(20, 5, name = "normal temperature")
membership(G, 25)
}
\keyword{math}
