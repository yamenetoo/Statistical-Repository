\name{entropy_threshold}
\alias{entropy_threshold}
\alias{entropy_partition}
\alias{print.entropy_threshold}
\alias{print.entropy_partition}
\title{Inductive Reasoning by Entropy Minimisation}
\description{
Find thresholds between classes of sample data by minimising the entropy of
the split (Section 4.3 of Ross, 2000), and partition the data repeatedly into
primary, secondary, tertiary, \dots thresholds.
}
\usage{
entropy_threshold(x, class, candidates = NULL,
                  pure = c("middle", "entropy", "none"), k = 1)
entropy_partition(x, class, levels = 2, candidates = NULL,
                  pure = c("middle", "entropy", "none"), k = 1)

\method{print}{entropy_threshold}(x, digits = 3, \dots)
\method{print}{entropy_partition}(x, digits = 4, \dots)
}
\arguments{
  \item{x}{numeric sample values; for the print methods, the object.}
  \item{class}{class label of each sample (two or more classes).}
  \item{candidates}{candidate thresholds; by default the mid-points between
    adjacent distinct sample values.}
  \item{pure}{what to do in a region that contains a single class, where the
    entropy carries no class information: \code{"middle"} splits at the
    middle candidate (as in the chapter), \code{"entropy"} at the minimum of
    \eqn{S}, \code{"none"} does not split.}
  \item{k}{the positive constant of Equation (6).}
  \item{levels}{number of partitioning levels: 1 gives the primary threshold
    PRI, 2 adds the secondary thresholds SEC1 and SEC2, 3 the tertiary
    thresholds TER1 to TER4, and so on.}
  \item{digits}{digits to print.}
  \item{\dots}{unused.}
}
\details{
For a candidate threshold \eqn{x} the samples with value \eqn{\le x} form the
\eqn{p} region and the others the \eqn{q} region. With Equations (7)-(13),
\deqn{S(x) = p(x) S_p(x) + q(x) S_q(x),}
\deqn{S_p(x) = -k \sum_j p_j(x) \ln p_j(x), \quad S_q(x) = -k \sum_j q_j(x) \ln q_j(x),}
\deqn{p_j(x) = \frac{n_j(x) + 1}{n(x) + 1}, \quad q_j(x) = \frac{N_j(x) + 1}{N(x) + 1}, \quad p(x) = \frac{n(x)}{n}, \quad q(x) = 1 - p(x),}
where \eqn{n_j(x)} and \eqn{N_j(x)} count the class-\eqn{j} samples in the two
regions, \eqn{n(x)} and \eqn{N(x)} all samples, and \eqn{n} the samples of the
region being split. The threshold is the candidate with minimum entropy.

\code{entropy_partition} applies \code{entropy_threshold} to the whole data
(PRI), then to the samples on each side of every threshold found.
}
\value{
\code{entropy_threshold} returns an object of class
\code{"entropy_threshold"} with the table of \eqn{p_j, q_j, p, q, S_p, S_q, S}
for every candidate (\code{table}), the chosen \code{threshold} and its
entropy \code{S}. \code{entropy_partition} returns an object of class
\code{"entropy_partition"} with a data frame \code{thresholds} (level, label,
threshold, entropy, number of samples) sorted by threshold, and the
individual steps.
}
\seealso{\code{\link{partition_sets}}, \code{\link{inductive_mf}},
  \code{\link{ross_ellipse}}}
\examples{
data(ross_ellipse)
## Table 1.3: the candidates of the chapter
entropy_threshold(ross_ellipse$x, ross_ellipse$class,
                  candidates = c(0.7, 1.0, 1.5, 3.45))
## all mid-points
et <- entropy_threshold(ross_ellipse$x, ross_ellipse$class)
et$threshold

## primary and secondary thresholds (Tables 1.3-1.5)
entropy_partition(ross_ellipse$x, ross_ellipse$class, levels = 2)
}
\keyword{math}
