\name{compare_bkwg}
\alias{compare_bkwg}
\title{Compare BKw-G Sub-Models and Baselines on a Data Set}
\description{
Fits several sub-models of the BKw-G family, for one or more baselines, and
returns a table with the log-likelihood, AIC, CAIC, BIC, HQIC, \eqn{W^*},
\eqn{A^*} and Kolmogorov-Smirnov statistics, sorted by AIC.
}
\usage{
compare_bkwg(x, baseline = "weibull",
             submodels = c("BKwG", "KwG", "BetaG", "EKwG", "BEG", "EG", "G"),
             method = "mle", ...)
}
\arguments{
  \item{x}{numeric vector of observations.}
  \item{baseline}{a character vector of built-in baseline names, or a single
    \code{\link{bkwg_baseline}} object.}
  \item{submodels}{the sub-models to fit (see \code{\link{bkwg_submodels}}).}
  \item{method}{estimation method passed to \code{\link{fit_bkwg}}.}
  \item{\dots}{further arguments passed to \code{\link{fit_bkwg}}.}
}
\value{
A data frame with one row per fitted model; the list of fitted models is
stored in the attribute \code{"fits"}.
}
\seealso{\code{\link{fit_bkwg}}, \code{\link{gof_bkwg}},
\code{\link{lrt_bkwg}}}
\examples{
set.seed(7)
x <- rbkwg(200, a = 2, b = 0.7, alpha = 3, beta = 0.8, par = c(1.2, 1))
tab <- compare_bkwg(x, baseline = "weibull",
                    submodels = c("BKwG", "KwG", "BetaG", "G"))
tab
\donttest{
## several baselines at once
compare_bkwg(x, baseline = c("weibull", "gamma", "lnorm"),
             submodels = c("BKwG", "G"))
}
}
\keyword{models}
