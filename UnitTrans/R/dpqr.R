#' Density, distribution, quantile and random generation
#'
#' Functions for the distribution of \eqn{V = g(Y)} on (0,1), where \eqn{Y}
#' follows a baseline distribution `dist` with parameters `par`.
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
#' @return Numeric vector.
#' @examples
#' dtrans(c(.2, .5, .8), "norm", c(0, 1), "logit")
#' ptrans(0.5, "gamma", c(2, 1), "exponential")
#' qtrans(c(.25, .5), "norm", c(0, 1), "probit")
#' set.seed(1); rtrans(5, "gamma", c(2, 1), "rational")
#' @name dpqr_trans
NULL

#' @rdname dpqr_trans
#' @export
dtrans <- function(v, dist, par, trans, log = FALSE) {
  tr <- get_transformation(trans)
  out <- rep(-Inf, length(v))
  ok <- is.finite(v) & v > 0 & v < 1
  vv <- v[ok]
  if (length(vv)) {
    y <- tr$ginv(vv)
    ly <- .call_base("d", dist, y, par, log = TRUE)
    out[ok] <- ly + tr$logjac(vv)
  }
  out[is.na(v)] <- NA
  if (log) out else exp(out)
}

#' @rdname dpqr_trans
#' @export
ptrans <- function(q, dist, par, trans, lower.tail = TRUE) {
  tr <- get_transformation(trans)
  out <- numeric(length(q))
  out[is.na(q)] <- NA
  lo <- !is.na(q) & q <= 0
  hi <- !is.na(q) & q >= 1
  mid <- !is.na(q) & !lo & !hi
  out[lo] <- 0
  out[hi] <- 1
  if (any(mid)) {
    y <- tr$ginv(q[mid])
    out[mid] <- if (tr$decreasing)
      .call_base("p", dist, y, par, lower.tail = FALSE)
    else .call_base("p", dist, y, par)
  }
  if (lower.tail) out else 1 - out
}

#' @rdname dpqr_trans
#' @export
qtrans <- function(p, dist, par, trans, lower.tail = TRUE) {
  tr <- get_transformation(trans)
  if (!lower.tail) p <- 1 - p
  pp <- if (tr$decreasing) 1 - p else p
  tr$g(.call_base("q", dist, pp, par))
}

#' @rdname dpqr_trans
#' @export
rtrans <- function(n, dist, par, trans) {
  qtrans(stats::runif(n), dist, par, trans)
}
