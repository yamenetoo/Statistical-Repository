\name{FuzzyMFD-package}
\alias{FuzzyMFD-package}
\alias{FuzzyMFD}
\docType{package}
\title{Membership Functions, Fuzzification and Defuzzification}
\description{
The methods of Ross (2000) for fuzzy sets on the real line.

\itemize{
  \item Fuzzy sets and shapes: \code{\link{fuzzy_set}}, \code{\link{fs_triangular}},
    \code{\link{fs_trapezoidal}}, \code{\link{fs_rectangular}}, \code{\link{fs_gaussian}},
    \code{\link{mf_triangular}} and the other \code{mf_*} functions.
  \item Features of the membership function (Section 2): \code{\link{fuzzy_core}},
    \code{\link{fuzzy_support}}, \code{\link{fuzzy_boundary}}, \code{\link{fuzzy_height}},
    \code{\link{is_normal}}, \code{\link{is_convex}}, \code{\link{is_fuzzy_number}},
    \code{\link{crossover_points}}, \code{\link{fuzzy_prototype}}.
  \item Operations: \code{\link{fuzzy_union}}, \code{\link{fuzzy_intersection}},
    \code{\link{fuzzy_complement}}, \code{\link{agreement}}.
  \item Fuzzification (Section 3): \code{\link{fuzzify}}.
  \item Membership value assignment (Section 4): \code{\link{rank_ordering}},
    \code{\link{entropy_threshold}}, \code{\link{entropy_partition}},
    \code{\link{partition_sets}}, \code{\link{inductive_mf}}.
  \item Defuzzification (Section 5): \code{\link{defuzzify}}, \code{\link{defuzz_max}},
    \code{\link{defuzz_centroid}}, \code{\link{defuzz_weighted}}, \code{\link{defuzz_meanmax}}.
  \item Worked examples: \code{\link{ross_colors}}, \code{\link{ross_ellipse}},
    \code{\link{ross2000_sets}}; the script
    \code{system.file("examples", "ross2000_examples.R", package = "FuzzyMFD")}
    reproduces the results of the chapter.
}

A fuzzy set is stored as a piecewise-linear membership function (knots and
membership values; zero outside the knots), so triangles, trapezoids, crisp
intervals and all their unions and intersections are represented exactly and
the centroid is computed in closed form.
}
\references{
Ross, T. J. (2000). Membership functions, fuzzification and defuzzification.
In P. S. Szczepaniak, P. J. G. Lisboa and J. Kacprzyk (eds.),
\emph{Fuzzy Systems in Medicine}, 49--78. Springer-Verlag, Berlin Heidelberg.

Zadeh, L. A. (1965). Fuzzy sets. \emph{Information and Control}, 8, 338--353.
}
\author{Mohamad Yamen AL Mohamad \email{yamenmohamad@tabrizu.ac.ir}}
\keyword{package}
