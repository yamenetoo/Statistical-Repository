# =========================================================================
# Kumaraswamy-Lindley (KwL) distribution
# Baseline: Lindley(theta) with
#   g(x;theta) = theta^2/(theta+1) * (1+x) * exp(-theta*x),         x>0
#   G(x;theta) = 1 - (1 + theta*x/(theta+1)) * exp(-theta*x)
#
# Kumaraswamy-G construction (Cordeiro & de Castro, 2011):
#   F(x) = 1 - (1 - G(x)^a)^b
#   f(x) = a*b*g(x)*G(x)^(a-1)*(1-G(x)^a)^(b-1)
#
# Parameters: theta > 0 (Lindley scale/shape), a > 0, b > 0 (Kumaraswamy shapes)
# =========================================================================

# NOTE: requires lambertW.R (lambertWm1) to be sourced beforehand.

# ---- Lindley baseline ----------------------------------------------------
# Lindley quantile function via Lambert W (lower branch), p = 0 case of the
# Lindley-Geometric quantile in Al-Lami, Jabbari Khamnei & Heydari (2025).
#' Lindley density (baseline of the KwL family)
#' @param x Non-negative numeric vector.
#' @param theta Lindley parameter, \eqn{\theta>0}.
#' @return Numeric vector of densities.
#' @export
gLindley <- function(x, theta) {
  (theta^2 / (theta + 1)) * (1 + x) * exp(-theta * x)
}

#' Lindley cumulative distribution function (baseline of the KwL family)
#' @param x Non-negative numeric vector.
#' @param theta Lindley parameter, \eqn{\theta>0}.
#' @return Numeric vector of cumulative probabilities.
#' @export
GLindley <- function(x, theta) {
  1 - (1 + theta * x / (theta + 1)) * exp(-theta * x)
}

#' Lindley quantile function via the Lambert W function
#'
#' Closed-form inverse of \code{\link{GLindley}}, obtained as the \eqn{p=0}
#' special case of the Lindley-Geometric quantile function of Al-Lami,
#' Jabbari Khamnei & Heydari (2025).
#'
#' @param v Numeric vector of probabilities in \eqn{(0,1)}.
#' @param theta Lindley parameter, \eqn{\theta>0}.
#' @return Numeric vector of Lindley quantiles.
#' @export
QLindley <- function(v, theta) {
  # With p = 0 (pure Lindley case of the Lindley-Geometric quantile),
  # the general LG quantile W_{-1}(-(v-1)(theta+1)e^{-(theta+1)}/(p*v-1))
  # reduces to W_{-1}((v-1)(theta+1)e^{-(theta+1)}).
  arg <- (v - 1) * (theta + 1) * exp(-(theta + 1))
  Wm1 <- lambertWm1(arg)
  -(Wm1 + theta + 1) / theta
}

# ---- Kumaraswamy-Lindley pdf / cdf / quantile / rng ----------------------

#' Density of the Kumaraswamy-Lindley (KwL) distribution
#'
#' \eqn{f(x;\theta,a,b) = a\,b\,g(x;\theta)\,G(x;\theta)^{a-1}\,[1-G(x;\theta)^a]^{b-1}},
#' where \eqn{g(\cdot;\theta)}, \eqn{G(\cdot;\theta)} are the Lindley density
#' and cdf. Reduces to the Lindley density when \code{a = b = 1}.
#'
#' @param x Non-negative numeric vector of quantiles.
#' @param theta Lindley parameter, \eqn{\theta>0}.
#' @param a,b Kumaraswamy shape parameters, \eqn{a,b>0}.
#' @param log Logical; if \code{TRUE}, the log-density is returned.
#' @return Numeric vector of (log-)densities.
#' @examples
#' dKwL(1.5, theta = 1, a = 1.5, b = 1.5)
#' integrate(dKwL, 0, Inf, theta = 1, a = 1.5, b = 1.5)$value  # ~= 1
#' @export
dKwL <- function(x, theta, a, b, log = FALSE) {
  Gx <- GLindley(x, theta)
  gx <- gLindley(x, theta)
  Gx <- pmin(pmax(Gx, 1e-300), 1 - 1e-15)
  logf <- log(a) + log(b) + log(gx) + (a - 1) * log(Gx) + (b - 1) * log1p(-Gx^a)
  if (log) return(logf)
  exp(logf)
}

#' Cumulative distribution function of the KwL distribution
#'
#' \eqn{F(x;\theta,a,b) = 1-[1-G(x;\theta)^a]^b}.
#'
#' @inheritParams dKwL
#' @return Numeric vector of cumulative probabilities.
#' @export
pKwL <- function(x, theta, a, b) {
  Gx <- GLindley(x, theta)
  Gx <- pmin(pmax(Gx, 0), 1 - 1e-15)
  1 - (1 - Gx^a)^b
}

# Quantile function of KwL: solve F(x) = u
#   (1 - G(x)^a)^b = 1-u  =>  G(x) = [1 - (1-u)^(1/b)]^(1/a)  =>  x = Q_Lindley(G(x))
#' Quantile function of the KwL distribution
#'
#' Closed-form quantile function obtained by inverting the Kumaraswamy-G cdf
#' and then the Lindley cdf via the Lambert W function (see
#' \code{\link{QLindley}}).
#'
#' @param u Numeric vector of probabilities in \eqn{(0,1)}.
#' @param theta,a,b KwL parameters as in \code{\link{dKwL}}.
#' @return Numeric vector of quantiles.
#' @examples
#' qKwL(0.5, theta = 1, a = 1.5, b = 1.5)
#' pKwL(qKwL(0.5, 1, 1.5, 1.5), 1, 1.5, 1.5)  # returns 0.5
#' @export
qKwL <- function(u, theta, a, b) {
  v <- (1 - (1 - u)^(1 / b))^(1 / a)
  v <- pmin(pmax(v, 1e-14), 1 - 1e-14)
  QLindley(v, theta)
}

#' Random generation from the KwL distribution
#'
#' Generates variates by the inversion method: draws \eqn{U\sim
#' \mathrm{Unif}(0,1)} and returns \code{qKwL(U, theta, a, b)}.
#'
#' @param n Number of variates to generate.
#' @param theta,a,b KwL parameters as in \code{\link{dKwL}}.
#' @return Numeric vector of length \code{n}.
#' @export
rKwL <- function(n, theta, a, b) {
  u <- stats::runif(n, 0, 1)
  qKwL(u, theta, a, b)
}
