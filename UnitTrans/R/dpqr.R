#' Density, distribution, quantile and random generation
#'
#' Functions for the distribution of \eqn{V = g(Y)} on (0,1), where \eqn{Y}
#' follows a baseline distribution `dist` with parameters `par`.
#'
#' Tail probabilities are computed directly from the baseline (for a
#' decreasing transformation the upper tail of \eqn{Y} is the lower tail of
#' \eqn{V}), never as `1 - p`, so `lower.tail = FALSE` and `log.p = TRUE` stay
#' accurate far into the tails.
#'
#' @param v,q Vector of quantiles in (0,1).
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param dist Baseline distribution name, e.g. `"gamma"` (uses `dgamma`,
#'   `pgamma`, `qgamma`).
#' @param par Baseline parameters: named vector/list, or unnamed in the order of
#'   the baseline's formals.
#' @param trans Transformation name (see [trans_list()]) or object.
#' @param log Logical; return log-density.
#' @param lower.tail Logical; if `FALSE`, return \eqn{P(V > q)}.
#' @param log.p Logical; probabilities are given/returned as `log(p)`.
#' @return Numeric vector. `rtrans` may return values numerically equal to 0
#'   or 1 for extreme draws.
#' @examples
#' dtrans(c(.2, .5, .8), "norm", c(0, 1), "logit")
#' ptrans(0.5, "gamma", c(2, 1), "exponential")
#' qtrans(c(.25, .5), "norm", c(0, 1), "probit")
#' set.seed(1); rtrans(5, "gamma", c(2, 1), "rational")
#' @name dpqr_trans
NULL

.clip01 <- function(x, eps = 1e-10) pmin(pmax(x, eps), 1 - eps)

#' @rdname dpqr_trans
#' @export
dtrans <- function(v, dist, par, trans, log = FALSE) {
  tr <- get_transformation(trans)
  bl <- .as_base(dist)
  par <- .name_par(par, bl)
  out <- rep(-Inf, length(v))
  ok <- is.finite(v) & v > 0 & v < 1
  vv <- v[ok]
  if (length(vv)) {
    y <- tr$ginv(vv)
    out[ok] <- .call_d(bl, y, par) + tr$logjac(vv)
  }
  out[is.na(v)] <- NA
  if (log) out else exp(out)
}

#' @rdname dpqr_trans
#' @export
ptrans <- function(q, dist, par, trans, lower.tail = TRUE, log.p = FALSE) {
  tr <- get_transformation(trans)
  bl <- .as_base(dist)
  par <- .name_par(par, bl)
  out <- numeric(length(q))
  out[is.na(q)] <- NA
  lo <- !is.na(q) & q <= 0
  hi <- !is.na(q) & q >= 1
  mid <- !is.na(q) & !lo & !hi
  p_lo <- if (lower.tail) 0 else 1
  out[lo] <- p_lo
  out[hi] <- 1 - p_lo
  if (any(mid)) {
    y <- tr$ginv(q[mid])
    # P(V <= q) = P(Y <= y) if g increases, P(Y >= y) if g decreases
    lt_b <- if (tr$decreasing) !lower.tail else lower.tail
    out[mid] <- .call_p(bl, y, par, lower.tail = lt_b, log.p = log.p)
  }
  if (log.p) out[lo | hi] <- log(out[lo | hi])
  out
}

#' @rdname dpqr_trans
#' @export
qtrans <- function(p, dist, par, trans, lower.tail = TRUE, log.p = FALSE) {
  tr <- get_transformation(trans)
  bl <- .as_base(dist)
  par <- .name_par(par, bl)
  lt_b <- if (tr$decreasing) !lower.tail else lower.tail
  tr$g(.call_q(bl, p, par, lower.tail = lt_b, log.p = log.p))
}

#' @rdname dpqr_trans
#' @export
rtrans <- function(n, dist, par, trans) {
  qtrans(stats::runif(n), dist, par, trans)
}
