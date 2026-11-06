\name{gof_bkwg}
\alias{gof_bkwg}
\alias{print.bkwggof}
\title{Goodness-of-Fit Statistics for Fitted BKw-G Models}
\description{
Computes the Cramer-von Mises (\eqn{W^*}) and Anderson-Darling (\eqn{A^*})
statistics of Chen and Balakrishnan (1995) and the Kolmogorov-Smirnov
statistic for a fitted model, optionally with parametric bootstrap p-values
that account for parameter estimation.
}
\usage{
gof_bkwg(fit, B = 0, seed = NULL)

\method{print}{bkwggof}(x, digits = 4, ...)
}
\arguments{
  \item{fit}{an object returned by \code{\link{fit_bkwg}}.}
  \item{B}{number of parametric bootstrap replicates; with \code{B = 0} only
    the asymptotic Kolmogorov-Smirnov p-value is reported.}
  \item{seed}{optional random seed for the bootstrap.}
  \item{x}{an object of class \code{"bkwggof"}.}
  \item{digits}{number of significant digits.}
  \item{\dots}{ignored.}
}
\details{
In the parametric bootstrap, samples are drawn from the fitted model and the
model is re-estimated on each sample with the same method, so the
p-values take the estimation of the parameters into account. Smaller
statistics indicate a better fit.
}
\value{
An object of class \code{"bkwggof"}: a list with a data frame
\code{table} (statistics and p-values), the number \code{B} of successful
bootstrap replicates and the fit.
}
\references{
Chen, G. and Balakrishnan, N. (1995). A general purpose approximate
goodness-of-fit test. \emph{Journal of Quality Technology}, \bold{27},
154--161.
}
\seealso{\code{\link{compare_bkwg}}, \code{\link{fit_bkwg}}}
\examples{
set.seed(10)
x <- rbkwg(150, a = 2, b = 1, alpha = 0.5, beta = 2, par = c(1.5, 1))
fit <- fit_bkwg(x, baseline = "weibull", submodel = "BetaG")
gof_bkwg(fit)
\donttest{
gof_bkwg(fit, B = 50, seed = 1)
}
}
\keyword{htest}
