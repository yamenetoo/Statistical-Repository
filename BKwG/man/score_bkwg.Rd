\name{score_bkwg}
\alias{score_bkwg}
\title{Analytical Score Vector of the BKw-G Shape Parameters}
\description{
Computes the analytical score (gradient of the log-likelihood) with respect
to the shape parameters \eqn{(a, b, \alpha, \beta)}{(a, b, alpha, beta)} for
given baseline parameters.
}
\usage{
score_bkwg(x, a = 1, b = 1, alpha = 1, beta = 1,
           baseline = "weibull", par = NULL, sum = TRUE)
}
\arguments{
  \item{x}{vector of observations.}
  \item{a, b, alpha, beta}{positive shape parameters.}
  \item{baseline, par}{baseline distribution and its parameters, as in
    \code{\link{dbkwg}}.}
  \item{sum}{logical; if \code{TRUE} return the total score, otherwise the
    \code{length(x)} by 4 matrix of individual contributions.}
}
\details{
With \eqn{T = G^a}, \eqn{L = \log(1 - T)} and \eqn{K = 1 - (1 - T)^b}, the
individual score components are
\deqn{U_a = 1/a + \log G - (b\beta - 1) T \log G / (1 - T) + (\alpha - 1) b (1-T)^{b-1} T \log G / K,}{U_a = 1/a + log G - (b beta - 1) T log G / (1 - T) + (alpha - 1) b (1-T)^(b-1) T log G / K,}
\deqn{U_b = 1/b + \beta L - (\alpha - 1)(1 - K) L / K,}{U_b = 1/b + beta L - (alpha - 1)(1 - K) L / K,}
\deqn{U_\alpha = \log K - \psi(\alpha) + \psi(\alpha + \beta),}{U_alpha = log K - psi(alpha) + psi(alpha + beta),}
\deqn{U_\beta = b L - \psi(\beta) + \psi(\alpha + \beta),}{U_beta = b L - psi(beta) + psi(alpha + beta),}
where \eqn{\psi} is the digamma function. Setting the total score to zero
gives the likelihood equations for the shape parameters.
}
\value{
A named numeric vector of length 4 (\code{sum = TRUE}) or a matrix.
}
\seealso{\code{\link{fit_bkwg}}}
\examples{
set.seed(2)
x <- rbkwg(200, a = 2, b = 1.5, alpha = 0.8, beta = 1.2, par = c(1.5, 1))
score_bkwg(x, a = 2, b = 1.5, alpha = 0.8, beta = 1.2, par = c(1.5, 1))
}
\keyword{distribution}
