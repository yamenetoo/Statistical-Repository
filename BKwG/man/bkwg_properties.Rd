\name{bkwg_properties}
\alias{bkwg_properties}
\alias{mbkwg}
\alias{mgf_bkwg}
\alias{bkwg_stats}
\alias{entropy_bkwg}
\alias{mrl_bkwg}
\title{Moments, Entropies and Other Properties of the BKw-G Distribution}
\description{
Numerical computation of raw moments, the moment generating function,
summary measures (mean, variance, skewness, kurtosis, Bowley skewness and
Moors kurtosis), Shannon and Renyi entropies and the mean residual life.
}
\usage{
mbkwg(r = 1, a = 1, b = 1, alpha = 1, beta = 1,
      baseline = "weibull", par = NULL, ...)
mgf_bkwg(t, a = 1, b = 1, alpha = 1, beta = 1,
         baseline = "weibull", par = NULL, ...)
bkwg_stats(a = 1, b = 1, alpha = 1, beta = 1,
           baseline = "weibull", par = NULL, ...)
entropy_bkwg(a = 1, b = 1, alpha = 1, beta = 1,
             baseline = "weibull", par = NULL,
             type = c("shannon", "renyi"), delta = 2, ...)
mrl_bkwg(t, a = 1, b = 1, alpha = 1, beta = 1,
         baseline = "weibull", par = NULL, ...)
}
\arguments{
  \item{r}{vector of moment orders.}
  \item{t}{vector of arguments (of the generating function or of the mean
    residual life).}
  \item{a, b, alpha, beta}{positive shape parameters.}
  \item{baseline, par}{baseline distribution and its parameters, as in
    \code{\link{dbkwg}}.}
  \item{type}{\code{"shannon"} or \code{"renyi"} entropy.}
  \item{delta}{order of the Renyi entropy (\eqn{\delta > 0}, \eqn{\delta \neq 1}{delta != 1}).}
  \item{\dots}{further arguments passed to \code{\link{integrate}}, such as
    \code{rel.tol}.}
}
\details{
All quantities are computed through the quantile function, e.g.
\eqn{E[X^r] = \int_0^1 Q(u)^r du}{E[X^r] = int_0^1 Q(u)^r du},
\eqn{H = -\int_0^1 \log f(Q(u)) du}{H = -int_0^1 log f(Q(u)) du},
\eqn{I_R(\delta) = (1-\delta)^{-1} \log \int_0^1 f(Q(u))^{\delta - 1} du}{I_R(delta) = log( int_0^1 f(Q(u))^(delta - 1) du ) / (1 - delta)}
and \eqn{m(t) = \int_{F(t)}^1 Q(u) du / (1 - F(t)) - t}{m(t) = int_F(t)^1 Q(u) du / (1 - F(t)) - t}.
A value \code{NA} (with a warning) is returned when the integral does not
converge, for example when the moment does not exist.
}
\value{
Numeric vectors. \code{bkwg_stats} returns a named vector with
\code{mean}, \code{variance}, \code{sd}, \code{skewness}, \code{kurtosis}
(not excess), \code{median}, \code{bowley} and \code{moors}.
}
\seealso{\code{\link{dbkwg}}, \code{\link{bkwg_mixture}}}
\examples{
mbkwg(1:2, a = 2, b = 1, alpha = 1.5, beta = 2, par = c(2, 1))
bkwg_stats(a = 2, b = 1, alpha = 1.5, beta = 2, par = c(2, 1))
entropy_bkwg(a = 2, b = 1, alpha = 1.5, beta = 2, par = c(2, 1))
entropy_bkwg(a = 2, b = 1, alpha = 1.5, beta = 2, par = c(2, 1),
             type = "renyi", delta = 0.5)
mgf_bkwg(c(-1, 0.5), a = 2, b = 1, alpha = 1.5, beta = 2, par = c(2, 1))
mrl_bkwg(c(0.5, 1), a = 2, b = 1, alpha = 1.5, beta = 2, par = c(2, 1))
}
\keyword{distribution}
