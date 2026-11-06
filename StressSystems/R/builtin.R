#' Topp-Leone distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the Topp-Leone distribution on (0, 1) with
#' \eqn{F(x) = [x(2 - x)]^{\nu}}, \eqn{\nu > 0}.
#'
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param shape Shape parameter \eqn{\nu > 0} (a single value).
#' @param log,log.p Logical; return log values.
#' @param lower.tail Logical; if \code{TRUE} probabilities are \eqn{P(X \le x)}.
#' @name toppleone
#' @export
dtoppleone <- function(x, shape, log = FALSE) {
  .check_pos(shape, "shape")
  out <- rep(-Inf, length(x))
  ok <- !is.na(x) & x > 0 & x < 1
  xo <- x[ok]
  out[ok] <- log(2 * shape) + log1p(-xo) + (shape - 1) * log(xo * (2 - xo))
  out[is.na(x)] <- NA_real_
  if (log) out else exp(out)
}

#' @rdname toppleone
#' @export
ptoppleone <- function(q, shape, lower.tail = TRUE, log.p = FALSE) {
  .check_pos(shape, "shape")
  qq <- pmin(pmax(q, 0), 1)
  out <- (qq * (2 - qq))^shape
  if (!lower.tail) out <- 1 - out
  if (log.p) log(out) else out
}

#' @rdname toppleone
#' @export
qtoppleone <- function(p, shape) {
  .check_pos(shape, "shape")
  1 - sqrt(1 - p^(1 / shape))
}

#' @rdname toppleone
#' @export
rtoppleone <- function(n, shape) qtoppleone(runif(n), shape)

#' Generalized Rayleigh (Burr type X) distribution
#'
#' Two-parameter generalized Rayleigh distribution of Kundu and Raqab with
#' \eqn{F(x) = (1 - e^{-(\lambda x)^2})^{\alpha}}, \eqn{x > 0}.
#' Here \code{shape} is \eqn{\alpha} and \code{rate} is \eqn{\lambda}.
#' Note that \code{VGAM::dgenray} uses \code{scale = 1 / rate}.
#'
#' @inheritParams toppleone
#' @param rate Rate (inverse scale) parameter \eqn{\lambda > 0}.
#' @name genrayleigh
#' @export
dgenrayleigh <- function(x, shape, rate = 1, log = FALSE) {
  .check_pos(shape, "shape")
  .check_pos(rate, "rate")
  out <- rep(-Inf, length(x))
  ok <- !is.na(x) & x > 0
  xo <- x[ok]
  z <- (rate * xo)^2
  out[ok] <- log(2) + log(shape) + 2 * log(rate) + log(xo) - z +
    (shape - 1) * log(-expm1(-z))
  out[is.na(x)] <- NA_real_
  if (log) out else exp(out)
}

#' @rdname genrayleigh
#' @export
pgenrayleigh <- function(q, shape, rate = 1, lower.tail = TRUE, log.p = FALSE) {
  .check_pos(shape, "shape")
  .check_pos(rate, "rate")
  out <- (-expm1(-(rate * pmax(q, 0))^2))^shape
  if (!lower.tail) out <- 1 - out
  if (log.p) log(out) else out
}

#' @rdname genrayleigh
#' @export
qgenrayleigh <- function(p, shape, rate = 1) {
  .check_pos(shape, "shape")
  .check_pos(rate, "rate")
  sqrt(-log1p(-p^(1 / shape))) / rate
}

#' @rdname genrayleigh
#' @export
rgenrayleigh <- function(n, shape, rate = 1) qgenrayleigh(runif(n), shape, rate)

.check_pos <- function(v, nm) {
  if (length(v) != 1 || !is.finite(v) || v <= 0)
    stop("'", nm, "' must be a single positive number")
  invisible(TRUE)
}

# Families shipped with the package (support and constraints are known).
.ss_builtin <- function() {
  list(
    toppleone = list(d = dtoppleone, p = ptoppleone, q = qtoppleone,
                     r = rtoppleone, lower = 0, upper = 1,
                     constraints = list(shape = c(0, Inf))),
    genrayleigh = list(d = dgenrayleigh, p = pgenrayleigh, q = qgenrayleigh,
                       r = rgenrayleigh, lower = 0, upper = Inf,
                       constraints = list(shape = c(0, Inf), rate = c(0, Inf)))
  )
}
