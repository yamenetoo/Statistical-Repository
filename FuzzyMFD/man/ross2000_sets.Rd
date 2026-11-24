\name{ross2000_sets}
\alias{ross2000_sets}
\title{Example Fuzzy Sets of Ross (2000)}
\description{Fuzzy sets used in the figures and examples of the chapter.}
\usage{
ross2000_sets(which = c("map", "output", "height", "voltage"))
}
\arguments{
  \item{which}{\code{"map"}: the three measurements \eqn{B_1, B_2, B_3} of
    mean arterial pressure in mm Hg (Figures 1.20-1.22); \code{"output"}: the
    two parts \eqn{C_1, C_2} of a fuzzy output (Figure 1.15);
    \code{"height"}: the crisp set \eqn{A = [5, 7]} and the fuzzy set
    \eqn{H} "near 6 feet" (Figure 1.1); \code{"voltage"}: illustrative
    "low voltage" and "medium voltage" sets (Figure 1.7).}
}
\value{A named list of \code{"fuzzyset"} objects.}
\examples{
B <- ross2000_sets("map")
plot_fuzzy_sets(B, xlab = "MAP (mm Hg)")
}
\keyword{datasets}
