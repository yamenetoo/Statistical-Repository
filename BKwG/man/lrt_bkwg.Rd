\name{lrt_bkwg}
\alias{lrt_bkwg}
\alias{wald_bkwg}
\title{Likelihood Ratio and Wald Tests for BKw-G Models}
\description{
\code{lrt_bkwg} tests a nested sub-model against a larger fitted model by
the likelihood ratio statistic. \code{wald_bkwg} performs a joint Wald test
of \eqn{H_0: \theta_j = \theta_{j0}}{H0: theta_j = theta_j0} for a set of
free parameters.
}
\usage{
lrt_bkwg(fit, restricted = "KwG")
wald_bkwg(fit, H0)
}
\arguments{
  \item{fit}{the larger model, an object returned by \code{\link{fit_bkwg}}
    with \code{method = "mle"}.}
  \item{restricted}{either the name of a sub-model (see
    \code{\link{bkwg_submodels}}), which is then fitted to the same data, or
    an already fitted restricted model.}
  \item{H0}{named vector or list of hypothesised values, e.g.
    \code{c(alpha = 1, beta = 1)}.}
}
\details{
Under \eqn{H_0}, the likelihood ratio statistic
\eqn{LR = 2(\ell_1 - \ell_0)}{LR = 2 (l1 - l0)} and the Wald statistic are
asymptotically chi-squared with degrees of freedom equal to the number of
restrictions. All restrictions of the sub-models lie in the interior of the
parameter space, so the usual chi-squared limit applies.
}
\value{An object of class \code{"htest"}.}
\seealso{\code{\link{fit_bkwg}}, \code{\link{compare_bkwg}}}
\examples{
set.seed(3)
x <- rbkwg(300, a = 3, b = 0.5, alpha = 2, beta = 0.6, par = c(1, 1))
fit <- fit_bkwg(x, baseline = "weibull")
lrt_bkwg(fit, "KwG")     # H0: alpha = beta = 1
lrt_bkwg(fit, "G")       # H0: the baseline Weibull is enough
wald_bkwg(fit, c(alpha = 1, beta = 1))
}
\keyword{htest}
