\name{fit_KwL}
\alias{fit_KwL}
\title{Maximum-likelihood fit of the KwL distribution}
\description{Maximizes the Kumaraswamy-Lindley log-likelihood via \code{\link[stats]{optim}} (L-BFGS-B), restarting from five dispersed starting values and keeping the best converged solution, to reduce sensitivity to local optima.}
\usage{fit_KwL(y, start = c(theta = 1, a = 1, b = 1))}
\arguments{
  \item{y}{Numeric vector of observations (\eqn{y_i > 0}).}
  \item{start}{Named numeric vector, initial values for \code{c(theta,a,b)} used as the first starting point.}
}
\value{A list with elements \code{theta}, \code{a}, \code{b} (the MLEs) and \code{value} (minimized negative log-likelihood), or \code{NULL} if no starting value converged.}
\examples{
set.seed(1)
y <- rKwL(200, theta = 1, a = 1.5, b = 2)
fit_KwL(y)
}
