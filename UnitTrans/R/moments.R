#' Raw moments, mean, variance, skewness and kurtosis
#'
#' Computed by numerical integration of the quantile function,
#' \eqn{E[V^r] = \int_0^1 Q(p)^r dp}.
#' @param dist,par,trans As in [dtrans()].
#' @param r Order(s) of the raw moments.
#' @return `moments_trans` returns a named vector of raw moments;
#'   `summary_trans` returns mean, variance, skewness and excess kurtosis.
#' @examples
#' moments_trans("gamma", c(2, 1), "rational", r = 1:2)
#' summary_trans("norm", c(0, 1), "logit")
#' @export
moments_trans <- function(dist, par, trans, r = 1:4) {
  out <- vapply(r, function(k)
    stats::integrate(function(p) qtrans(p, dist, par, trans)^k, 0, 1,
                     rel.tol = 1e-8, subdivisions = 500L)$value, 0)
  stats::setNames(out, paste0("m", r))
}

#' @rdname moments_trans
#' @export
summary_trans <- function(dist, par, trans) {
  m <- moments_trans(dist, par, trans, 1:4)
  mu <- m[[1]]; v <- m[[2]] - mu^2
  c3 <- m[[3]] - 3 * mu * m[[2]] + 2 * mu^3
  c4 <- m[[4]] - 4 * mu * m[[3]] + 6 * mu^2 * m[[2]] - 3 * mu^4
  c(mean = mu, variance = v, skewness = c3 / v^1.5, ex.kurtosis = c4 / v^2 - 3)
}

#' Survival and hazard functions
#'
#' @param v Values in (0,1).
#' @param dist,par,trans As in [dtrans()].
#' @return Numeric vector.
#' @examples
#' hazard_trans(c(.2, .5), "gamma", c(2, 1), "rational")
#' @export
hazard_trans <- function(v, dist, par, trans)
  dtrans(v, dist, par, trans) / survival_trans(v, dist, par, trans)

#' @rdname hazard_trans
#' @export
survival_trans <- function(v, dist, par, trans)
  ptrans(v, dist, par, trans, lower.tail = FALSE)
