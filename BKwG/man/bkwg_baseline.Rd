\name{bkwg_baseline}
\alias{bkwg_baseline}
\alias{bkwg_baselines}
\alias{print.bkwg_baseline}
\title{Baseline Distributions for the BKw-G Family}
\description{
\code{bkwg_baselines} lists the built-in baseline distributions \eqn{G}.
\code{bkwg_baseline} builds a user-defined baseline that can be used in every
function of the package.
}
\usage{
bkwg_baselines()

bkwg_baseline(name, pG, dG, qG = NULL, par = numeric(0),
              lower = NULL, upper = NULL, support = c(-Inf, Inf),
              start = NULL)

\method{print}{bkwg_baseline}(x, ...)
}
\arguments{
  \item{name}{a short name for the baseline.}
  \item{pG}{the baseline cdf, a function \code{function(q, par)} or
    \code{function(q, par, lower.tail = TRUE, log.p = FALSE)}; \code{par} is
    a named numeric vector. Supplying \code{log.p} improves tail accuracy.}
  \item{dG}{the baseline pdf, \code{function(x, par)} or
    \code{function(x, par, log = FALSE)}.}
  \item{qG}{optional quantile function \code{function(p, par)}. If
    \code{NULL}, it is obtained by numerical inversion of \code{pG}.}
  \item{par}{named numeric vector of default parameter values (its names
    define the baseline parameters).}
  \item{lower, upper}{lower and upper bounds of the parameters (used for
    validation and by the optimiser). Default: unbounded.}
  \item{support}{the support \code{c(lower, upper)} of the baseline.}
  \item{start}{optional function \code{function(x)} returning named starting
    values computed from data, or a fixed named vector.}
  \item{x}{a \code{"bkwg_baseline"} object.}
  \item{\dots}{ignored.}
}
\details{
Built-in baselines (parameter names in brackets): \code{"exp"} (rate),
\code{"weibull"} (shape, scale), \code{"norm"} (mean, sd), \code{"lnorm"}
(meanlog, sdlog), \code{"gamma"} (shape, rate), \code{"logis"} (location,
scale), \code{"gumbel"} (location, scale; maximum Gumbel),
\code{"lomax"} (shape, scale; \eqn{G(x) = 1 - (1 + x/scale)^{-shape}}{G(x) = 1 - (1 + x/scale)^(-shape)}),
\code{"llogis"} (shape, scale; log-logistic),
\code{"frechet"} (shape, scale; \eqn{G(x) = \exp\{-(x/scale)^{-shape}\}}{G(x) = exp(-(x/scale)^(-shape))}),
\code{"burrxii"} (c, k, scale; \eqn{G(x) = 1 - [1 + (x/scale)^c]^{-k}}{G(x) = 1 - [1 + (x/scale)^c]^(-k)})
and \code{"unif"} (no parameters, support \eqn{(0, 1)}; gives the
beta-Kumaraswamy distribution on the unit interval).
}
\value{
\code{bkwg_baselines} returns a data frame. \code{bkwg_baseline} returns an
object of class \code{"bkwg_baseline"}.
}
\seealso{\code{\link{dbkwg}}, \code{\link{fit_bkwg}}}
\examples{
bkwg_baselines()

## a Rayleigh baseline written by the user
rayleigh <- bkwg_baseline("rayleigh",
  pG = function(q, par) ifelse(q > 0, 1 - exp(-q^2 / (2 * par[["sigma"]]^2)), 0),
  dG = function(x, par) ifelse(x > 0, x / par[["sigma"]]^2 *
                               exp(-x^2 / (2 * par[["sigma"]]^2)), 0),
  qG = function(p, par) par[["sigma"]] * sqrt(-2 * log1p(-p)),
  par = c(sigma = 1), lower = 0, upper = Inf, support = c(0, Inf),
  start = function(x) c(sigma = sqrt(mean(x^2) / 2)))
rayleigh
dbkwg(1:3, a = 2, b = 1, alpha = 1.5, beta = 0.8, baseline = rayleigh,
      par = c(sigma = 2))
}
\keyword{distribution}
