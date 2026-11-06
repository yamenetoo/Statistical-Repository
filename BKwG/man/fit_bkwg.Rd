\name{fit_bkwg}
\alias{fit_bkwg}
\alias{bkwg_submodels}
\title{Fit BKw-G Models and Their Sub-Models}
\description{
Estimates the parameters of a BKw-G distribution, or of one of its nested
sub-models, for any baseline, by maximum likelihood, maximum product of
spacings, minimum Anderson-Darling or Cramer-von Mises distance, ordinary or
weighted least squares.
}
\usage{
fit_bkwg(x, baseline = "weibull", submodel = "BKwG", fixed = NULL,
         start = NULL, method = c("mle", "mps", "ad", "cvm", "lse", "wlse"),
         nstart = 3, hessian = TRUE, control = list())

bkwg_submodels()
}
\arguments{
  \item{x}{numeric vector of observations.}
  \item{baseline}{baseline distribution: a built-in name (see
    \code{\link{bkwg_baselines}}) or a \code{\link{bkwg_baseline}} object.}
  \item{submodel}{one of \code{"BKwG"} (full model), \code{"KwG"},
    \code{"BetaG"}, \code{"EKwG"}, \code{"BEG"}, \code{"EG"} or \code{"G"};
    see \code{bkwg_submodels()}.}
  \item{fixed}{optional named vector or list of parameters held fixed at the
    given values (shape parameters \code{a}, \code{b}, \code{alpha},
    \code{beta} or baseline parameters). Overrides the restrictions of
    \code{submodel}.}
  \item{start}{optional named vector or list of starting values for free
    parameters. By default starting values for the baseline come from simple
    moment or quantile estimators, combined with a grid of shape values.}
  \item{method}{estimation method: \code{"mle"} (maximum likelihood),
    \code{"mps"} (maximum product of spacings), \code{"ad"} (minimum
    Anderson-Darling distance), \code{"cvm"} (minimum Cramer-von Mises
    distance), \code{"lse"} (least squares) or \code{"wlse"} (weighted least
    squares).}
  \item{nstart}{number of the best grid starting points from which the
    optimiser is run; the best solution is kept.}
  \item{hessian}{logical; compute the observed information matrix and the
    standard errors (methods \code{"mle"} and \code{"mps"} only).}
  \item{control}{list of control parameters for \code{\link{optim}}.}
}
\details{
Optimisation is carried out on an unconstrained scale (log for positive
parameters, logit for parameters with two finite bounds) with Nelder-Mead
followed by BFGS, from several starting points. Standard errors are based on
the inverse of the observed information matrix of the log-likelihood,
computed with \pkg{numDeriv} when it is installed and with
\code{\link{optimHess}} otherwise.

Sub-models of the BKw-G family:
\tabular{lll}{
  \code{"BKwG"} \tab Beta Kumaraswamy-G \tab no restriction \cr
  \code{"KwG"} \tab Kumaraswamy-G \tab \eqn{\alpha = \beta = 1}{alpha = beta = 1} \cr
  \code{"BetaG"} \tab Beta-G \tab \eqn{a = b = 1} \cr
  \code{"EKwG"} \tab Exponentiated Kumaraswamy-G \tab \eqn{\beta = 1}{beta = 1} \cr
  \code{"BEG"} \tab Beta exponentiated-G \tab \eqn{b = 1} \cr
  \code{"EG"} \tab Exponentiated-G \tab \eqn{b = \alpha = \beta = 1}{b = alpha = beta = 1} \cr
  \code{"G"} \tab Baseline \tab \eqn{a = b = \alpha = \beta = 1}{a = b = alpha = beta = 1}
}
Note that when \eqn{\alpha = 1}{alpha = 1} the parameters \eqn{b} and
\eqn{\beta}{beta} enter the model only through \eqn{b\beta}{b beta}, so the
full model is weakly identified when the data point to
\eqn{\alpha \approx 1}{alpha close to 1}; in that case a sub-model is
preferable. Scale-type baseline parameters can also be confounded with the
shape parameters (for instance the Weibull scale with \eqn{b}); always check
the standard errors and compare with the sub-models.
}
\value{
An object of class \code{"bkwgfit"}, a list with components
\item{estimate}{all parameters (free and fixed).}
\item{coefficients}{the estimated (free) parameters.}
\item{vcov, se}{covariance matrix and standard errors of the free
  parameters.}
\item{loglik, AIC, CAIC, BIC, HQIC}{log-likelihood at the estimate and
  information criteria.}
\item{fixed}{the fixed parameters.}
\item{method, submodel, baseline, n, k, x, convergence, call}{other
  information.}
\code{bkwg_submodels} returns a data frame describing the sub-models.
}
\seealso{
\code{\link{summary.bkwgfit}}, \code{\link{gof_bkwg}},
\code{\link{lrt_bkwg}}, \code{\link{compare_bkwg}}
}
\examples{
bkwg_submodels()

set.seed(123)
x <- rbkwg(200, a = 2, b = 1.5, baseline = "weibull",
           par = c(shape = 1.5, scale = 2))

## Kumaraswamy-Weibull sub-model with a known scale
fit_kw <- fit_bkwg(x, baseline = "weibull", submodel = "KwG",
                   fixed = c(scale = 2))
fit_kw
summary(fit_kw)
confint(fit_kw, type = "log")
\donttest{
## the full model and other estimation methods
fit <- fit_bkwg(x, baseline = "weibull")
summary(fit)
fit_mps <- fit_bkwg(x, baseline = "weibull", submodel = "KwG", method = "mps")
}
}
\keyword{models}
\keyword{distribution}
