#' unitTransformDist: Generalized Unit-Interval Distributions via the Odds Transform
#'
#' Constructs, fits, and diagnoses continuous probability distributions on
#' (0,1) obtained from an arbitrary base distribution on (0, Inf) via
#' \eqn{Y = 1/(1+X)}. The base distribution is supplied directly as
#' density/CDF/quantile/random-generation functions, so any distribution
#' the user can write down is supported -- not only a fixed registry of
#' presets.
#'
#' @section Main entry points:
#' \itemize{
#'   \item \code{\link{utd}} construct a distribution from any base d/p/q/r functions
#'   \item \code{\link{fit_utd_mle}} maximum likelihood estimation
#'   \item \code{\link{profile_likelihood_ci}} corrected profile-likelihood interval
#'   \item \code{\link{run_monte_carlo}} Monte Carlo study of estimator behaviour
#'   \item \code{\link{bias_correct}} unified bias-correction dispatcher
#' }
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats optim uniroot qchisq runif quantile sd var median qnorm ks.test density
#' @importFrom graphics plot abline barplot
## usethis namespace: end
NULL
