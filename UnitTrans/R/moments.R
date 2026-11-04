#' Raw moments, mean, variance, skewness and kurtosis
#'
#' Computed by numerical integration of the quantile function,
#' \eqn{E[V^r] = \int_0^1 Q(p)^r dp}. Central moments in `summary_trans` are
#' integrated directly (not obtained from raw moments), which avoids
#' cancellation errors.
#' @param dist,par,trans As in [dtrans()].
#' @param r Order(s) of the raw moments.
#' @return `moments_trans` returns a named vector of raw moments;
#'   `summary_trans` returns mean, variance, skewness and excess kurtosis.
#' @examples
#' moments_trans("gamma", c(2, 1), "rational", r = 1:2)
#' summary_trans("norm", c(0, 1), "logit")
#' @export
moments_trans <- function(dist, par, trans, r = 1:4) {
  tr <- get_transformation(trans); bl <- .as_base(dist); par <- .name_par(par, bl)
  out <- vapply(r, function(k)
    stats::integrate(function(p) qtrans(p, bl, par, tr)^k, 0, 1,
                     rel.tol = 1e-8, subdivisions = 500L)$value, 0)
  stats::setNames(out, paste0("m", r))
}

#' @rdname moments_trans
#' @export
summary_trans <- function(dist, par, trans) {
  tr <- get_transformation(trans); bl <- .as_base(dist); par <- .name_par(par, bl)
  mu <- stats::integrate(function(p) qtrans(p, bl, par, tr), 0, 1,
                         rel.tol = 1e-10, subdivisions = 500L)$value
  cm <- function(k)
    stats::integrate(function(p) (qtrans(p, bl, par, tr) - mu)^k, 0, 1,
                     rel.tol = 1e-8, subdivisions = 500L)$value
  v <- cm(2)
  c(mean = mu, variance = v, skewness = cm(3) / v^1.5, ex.kurtosis = cm(4) / v^2 - 3)
}

#' Survival and hazard functions
#'
#' The survival function is the upper tail computed directly from the baseline
#' and the hazard is evaluated on the log scale, so both stay accurate near 1.
#' @param v Values in (0,1).
#' @param dist,par,trans As in [dtrans()].
#' @return Numeric vector.
#' @examples
#' hazard_trans(c(.2, .5), "gamma", c(2, 1), "rational")
#' @export
hazard_trans <- function(v, dist, par, trans) {
  bl <- .as_base(dist); par <- .name_par(par, bl)
  exp(dtrans(v, bl, par, trans, log = TRUE) -
        ptrans(v, bl, par, trans, lower.tail = FALSE, log.p = TRUE))
}

#' @rdname hazard_trans
#' @export
survival_trans <- function(v, dist, par, trans)
  ptrans(v, dist, par, trans, lower.tail = FALSE)
