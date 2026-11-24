\name{defuzzify}
\alias{defuzzify}
\alias{defuzz_max}
\alias{defuzz_centroid}
\alias{defuzz_weighted}
\alias{defuzz_meanmax}
\title{Defuzzification}
\description{
Reduce a fuzzy output to a single crisp value (Section 5 of Ross, 2000) by the
max-membership (height), centroid, weighted-average or mean-max method.
}
\usage{
defuzzify(sets, method = c("centroid", "max", "meanmax", "weighted"), \dots)
defuzz_max(A, which = c("first", "last", "all"))
defuzz_centroid(A)
defuzz_weighted(sets, centers = NULL, weights = NULL)
defuzz_meanmax(A)
}
\arguments{
  \item{sets}{a fuzzy set or a list of fuzzy sets, the parts
    \eqn{C_1, \dots, C_k} of the output (Equation 14).}
  \item{method}{defuzzification method.}
  \item{\dots}{passed to \code{defuzz_max} or \code{defuzz_weighted}.}
  \item{A}{a \code{"fuzzyset"} object.}
  \item{which}{for \code{defuzz_max} when the maximum is not unique: the
    first or last maximising element, or the intervals of all of them.}
  \item{centers}{centre of each part; by default the middle of its plateau of
    maximum membership (its mean when the part is symmetric).}
  \item{weights}{weight of each part; by default its maximum membership.}
}
\details{
\describe{
  \item{max-membership (15)}{\eqn{z^*} with \eqn{\mu_C(z^*) \ge \mu_C(z)} for
    all \eqn{z}.}
  \item{centroid (16)}{\eqn{z^* = \int \mu_C(z) z\,dz / \int \mu_C(z)\,dz},
    computed in closed form for the piecewise-linear membership function.}
  \item{weighted average (17)}{\eqn{z^* = \sum \mu_C(\bar z) \bar z / \sum \mu_C(\bar z)},
    each part weighted by its maximum membership; meant for symmetric parts.}
  \item{mean-max (18)}{\eqn{z^* = (a + b)/2} where \eqn{[a, b]} spans the
    elements of maximum membership.}
}
\code{defuzzify} applies the centroid, max-membership and mean-max methods to
the union \eqn{C = \bigcup_i C_i} of the parts, and the weighted-average
method to the parts themselves.
}
\value{
A number (\code{defuzz_max(which = "all")} returns a matrix of intervals).
}
\note{
For the mean-arterial-pressure example of the chapter the exact centroid is
4.953 mm Hg; the chapter reports 4.9 mm Hg.
}
\seealso{\code{\link{fuzzy_union}}, \code{\link{ross2000_sets}}}
\examples{
B <- ross2000_sets("map")              # Figures 1.20-1.22
defuzzify(B, "centroid")               # Figure 1.23
defuzzify(B, "weighted")               # Figure 1.24: 5.41
defuzzify(B, "meanmax")                # Figure 1.25: 6.5
defuzzify(B, "max", which = "all")

plot_fuzzy_sets(list(union = fuzzy_union(B)), legend = FALSE)
abline(v = defuzzify(B), lty = 2)
}
\keyword{math}
