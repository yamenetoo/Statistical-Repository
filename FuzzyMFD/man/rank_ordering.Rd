\name{rank_ordering}
\alias{rank_ordering}
\alias{print.rank_ordering}
\alias{plot.rank_ordering}
\title{Membership Values by Rank Ordering}
\description{
Assign membership values from pairwise preferences (Section 4.2 of Ross,
2000): preferences of an individual, a committee or a poll are summed per item
and normalised.
}
\usage{
rank_ordering(P)

\method{print}{rank_ordering}(x, digits = 3, \dots)
\method{plot}{rank_ordering}(x, ylab = "membership", xlab = "",
     main = "Rank-ordering membership", \dots)
}
\arguments{
  \item{P}{square matrix of pairwise preferences: \code{P[i, j]} is the number
    of respondents preferring item \code{i} to item \code{j}; the diagonal is
    ignored.}
  \item{x}{a \code{"rank_ordering"} object.}
  \item{digits}{digits to print.}
  \item{ylab, xlab, main}{plot labels.}
  \item{\dots}{further arguments passed to \code{plot}.}
}
\details{
The row sums give the number of times each item was preferred; their share of
all responses (in percent) and the rank are reported as in Table 1.1 of the
chapter. The membership of each item in the fuzzy set "best" is its percentage
divided by the largest percentage (Figure 1.9). A warning is given when
\code{P[i, j] + P[j, i]} is not the same for all pairs.
}
\value{
An object of class \code{"rank_ordering"}: a list with \code{table}
(\code{item}, \code{total}, \code{percent}, \code{rank}, \code{membership}),
\code{preferences}, \code{responses} and the named vector \code{membership}.
}
\seealso{\code{\link{ross_colors}}}
\examples{
data(ross_colors)
ro <- rank_ordering(ross_colors)   # Table 1.1
ro
plot(ro)                           # Figure 1.9
}
\keyword{math}
