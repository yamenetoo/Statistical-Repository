\name{KwLBootCC}
\alias{KwLBootCC}
\title{Bootstrap percentile control chart for the KwL distribution}
\description{Implements the parametric-bootstrap percentile control-chart algorithm of Al-Lami, Jabbari Khamnei and Heydari (2025, PLOS ONE, for the Lindley-Geometric distribution), adapted to the three-parameter Kumaraswamy-Lindley distribution.}
\usage{
KwLBootCC(u, alpha, theta, a, b, n = 5, B = 2000, MC = 50,
          start = c(theta = 1, a = 1, b = 1), verbose = FALSE)
}
\arguments{
  \item{u}{Percentile of interest, \eqn{0 < u < 1}.}
  \item{alpha}{Nominal false-alarm probability (e.g. 0.0027).}
  \item{theta}{KwL Lindley parameter used to generate bootstrap subgroups.}
  \item{a}{KwL Kumaraswamy shape parameter used to generate bootstrap subgroups.}
  \item{b}{KwL Kumaraswamy shape parameter used to generate bootstrap subgroups.}
  \item{n}{Subgroup size.}
  \item{B}{Number of bootstrap iterations per Monte-Carlo repetition.}
  \item{MC}{Number of outer Monte-Carlo repetitions.}
  \item{start}{Starting values passed to \code{\link{fit_KwL}}.}
  \item{verbose}{If \code{TRUE}, prints progress for each Monte-Carlo repetition.}
}
\value{A list with \code{MLCL}, \code{MUCL} (mean control limits), \code{SDLCL}, \code{SDUCL} (their Monte-Carlo standard deviations), and \code{LCL_all}, \code{UCL_all} (the per-repetition limits).}
\examples{
\donttest{
set.seed(1)
KwLBootCC(u = 0.05, alpha = 0.0027, theta = 1, a = 1.5, b = 1.5, n = 5, B = 200, MC = 5)
}
}
