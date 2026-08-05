\name{KwL_ARL}
\alias{KwL_ARL}
\title{Average run length simulation for the KwL bootstrap control chart}
\description{Simulates Phase-II monitoring under the KwL bootstrap control chart and returns the average run length (ARL) and its standard deviation, either in-control (Phase-II parameters equal to Phase-I parameters) or out-of-control (a specified parameter shift).}
\usage{
KwL_ARL(u, LCL, UCL, theta1, a1, b1, n = 5, n_runs = 200, max_run = 2000,
        start = c(theta = 1, a = 1, b = 1))
}
\arguments{
  \item{u}{Percentile of interest.}
  \item{LCL}{Lower control limit (typically from \code{\link{KwLBootCC}}).}
  \item{UCL}{Upper control limit (typically from \code{\link{KwLBootCC}}).}
  \item{theta1}{KwL Lindley parameter generating the Phase-II subgroups.}
  \item{a1}{KwL Kumaraswamy shape parameter generating the Phase-II subgroups.}
  \item{b1}{KwL Kumaraswamy shape parameter generating the Phase-II subgroups.}
  \item{n}{Subgroup size.}
  \item{n_runs}{Number of independent run lengths to simulate.}
  \item{max_run}{Cap on the run length (to guarantee termination).}
  \item{start}{Starting values passed to \code{\link{fit_KwL}}.}
}
\value{A list with \code{ARL} (mean run length), \code{SDRL} (its standard deviation), and \code{RL_all} (the individual run lengths).}
\examples{
\donttest{
set.seed(1)
lim <- KwLBootCC(0.05, 0.0027, 1, 1.5, 1.5, n = 5, B = 200, MC = 5)
KwL_ARL(0.05, lim$MLCL, lim$MUCL, 1, 1.5, 1.5, n_runs = 10, max_run = 300)
}
}
