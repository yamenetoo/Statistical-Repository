\name{bkwgfit-methods}
\alias{bkwgfit-methods}
\alias{print.bkwgfit}
\alias{summary.bkwgfit}
\alias{print.summary.bkwgfit}
\alias{coef.bkwgfit}
\alias{vcov.bkwgfit}
\alias{logLik.bkwgfit}
\alias{nobs.bkwgfit}
\alias{confint.bkwgfit}
\alias{residuals.bkwgfit}
\alias{quantile.bkwgfit}
\alias{plot.bkwgfit}
\title{Methods for Fitted BKw-G Models}
\description{
Printing, summaries, coefficients, covariance matrix, log-likelihood,
confidence intervals, residuals, fitted quantiles and diagnostic plots for
objects returned by \code{\link{fit_bkwg}}.
}
\usage{
\method{print}{bkwgfit}(x, digits = 4, ...)
\method{summary}{bkwgfit}(object, ...)
\method{print}{summary.bkwgfit}(x, digits = 4, ...)
\method{coef}{bkwgfit}(object, complete = FALSE, ...)
\method{vcov}{bkwgfit}(object, ...)
\method{logLik}{bkwgfit}(object, ...)
\method{nobs}{bkwgfit}(object, ...)
\method{confint}{bkwgfit}(object, parm, level = 0.95,
        type = c("wald", "log"), ...)
\method{residuals}{bkwgfit}(object, type = c("quantile", "pit", "cox-snell"),
          ...)
\method{quantile}{bkwgfit}(x, probs = c(0.25, 0.5, 0.75), ...)
\method{plot}{bkwgfit}(x, y, which = 1:4, ...)
}
\arguments{
  \item{x, object}{an object of class \code{"bkwgfit"} (or
    \code{"summary.bkwgfit"} for its print method).}
  \item{digits}{number of significant digits.}
  \item{complete}{logical; if \code{TRUE} return fixed parameters as well.}
  \item{parm}{names or indices of parameters; all free parameters by
    default.}
  \item{level}{confidence level.}
  \item{type}{for \code{confint}: \code{"wald"} intervals on the original
    scale or \code{"log"} intervals (delta method on the log scale, which
    keep positive parameters positive). For \code{residuals}: normal
    quantile residuals \eqn{\Phi^{-1}(F(x))}{qnorm(F(x))},
    probability integral transform \eqn{F(x)} or Cox-Snell residuals
    \eqn{-\log(1 - F(x))}{-log(1 - F(x))}.}
  \item{probs}{probabilities for the fitted quantiles.}
  \item{y}{not used.}
  \item{which}{which plots to draw: 1 histogram with fitted density, 2
    empirical and fitted cdf, 3 P-P plot, 4 Q-Q plot.}
  \item{\dots}{further arguments (currently ignored).}
}
\value{
\code{summary} returns an object of class \code{"summary.bkwgfit"} with the
coefficient table (estimates, standard errors, Wald z statistics and
p-values for \eqn{H_0: \theta = 0}{H0: theta = 0}) and the information
criteria. The other methods return what their generics return.
}
\seealso{\code{\link{fit_bkwg}}, \code{\link{wald_bkwg}}}
\examples{
set.seed(1)
x <- rbkwg(200, a = 1.5, b = 2, alpha = 0.6, beta = 1, baseline = "exp",
           par = c(rate = 0.5))
fit <- fit_bkwg(x, baseline = "exp", submodel = "EKwG")
summary(fit)
coef(fit, complete = TRUE)
AIC(fit); BIC(fit)
confint(fit)
quantile(fit, c(0.1, 0.5, 0.9))
plot(fit)
}
\keyword{models}
