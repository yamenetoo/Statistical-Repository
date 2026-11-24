\name{partition_sets}
\alias{partition_sets}
\alias{inductive_mf}
\title{Membership Functions from Partition Thresholds}
\description{
Build overlapping membership functions whose toes are the thresholds found by
entropy minimisation (Figures 1.11, 1.13 and 1.14 of Ross, 2000).
}
\usage{
partition_sets(thresholds, range, labels = NULL, overlap = 0.5)
inductive_mf(x, class, levels = 2, range = base::range(x), labels = NULL,
             candidates = NULL, pure = c("middle", "entropy", "none"),
             overlap = 0.5)
}
\arguments{
  \item{thresholds}{numeric thresholds, or an \code{"entropy_partition"}
    object.}
  \item{range}{\code{c(lower, upper)}: the universe of the variable.}
  \item{labels}{labels of the sets; by default NG, PO (two sets), NG, ZE, PO
    (three), NB, NS, ZE, PS, PB (five) or NB, NM, NS, ZE, PS, PM, PB (seven).}
  \item{overlap}{for a single threshold: half-width of the overlap as a
    fraction of the distance from the threshold to the nearer end of the
    range.}
  \item{x, class, levels, candidates, pure}{passed to
    \code{\link{entropy_partition}}.}
}
\details{
With \eqn{k \ge 2} sorted thresholds \eqn{t_1 < \dots < t_k} on
\eqn{[lower, upper]}, the first set equals 1 on \eqn{[lower, t_1]} and falls to
0 at \eqn{t_2}; set \eqn{i} (\eqn{1 < i < k}) is the triangle with peak
\eqn{t_i} and feet \eqn{t_{i-1}, t_{i+1}}; the last set rises from 0 at
\eqn{t_{k-1}} to 1 at \eqn{t_k} and stays 1 up to \eqn{upper}. Neighbouring
sets therefore cross at membership 0.5. A single threshold (the primary
threshold alone) gives two sets that cross at membership 0.5 over it.
}
\value{
\code{partition_sets} returns a named list of \code{"fuzzyset"} objects;
\code{inductive_mf} returns a list with the \code{partition} and the
\code{sets}.
}
\seealso{\code{\link{entropy_partition}}, \code{\link{plot_fuzzy_sets}}}
\examples{
data(ross_ellipse)
res <- inductive_mf(ross_ellipse$x, ross_ellipse$class, levels = 2,
                    range = c(0, 100))
res$partition
plot_fuzzy_sets(res$sets, xlab = "x = a/b")   # Figure 1.14

## the first partition alone (Figure 1.13)
plot_fuzzy_sets(partition_sets(1.5, range = c(0, 100)))
}
\keyword{math}
