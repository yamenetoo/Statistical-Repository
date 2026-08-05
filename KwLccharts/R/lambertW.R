# Numerically robust Lambert W function (principal branch W0 and lower branch W_{-1})
# Implemented via series initialization + Halley's iteration (no external deps).

#' Lambert W function, principal branch (W0)
#'
#' Numerically evaluates the principal (upper) branch of the Lambert W
#' function, the solution W(x) of \eqn{W(x) e^{W(x)} = x}, via Halley's
#' iteration. Implemented without external dependencies.
#'
#' @param x Numeric vector, \eqn{x \ge -1/e}.
#' @return Numeric vector of W0(x) values.
#' @examples
#' lambertW0(1)   # 0.5671433
#' lambertW0(exp(1))  # 1
#' @export
lambertW0 <- function(x) {
  # valid for x >= -1/e
  out <- numeric(length(x))
  for (i in seq_along(x)) {
    xi <- x[i]
    if (is.na(xi)) { out[i] <- NA; next }
    if (xi == 0) { out[i] <- 0; next }
    # initial guess
    if (xi <= 1) {
      w <- xi * (1 - xi + 1.5*xi^2)  # series near 0
      if (xi > -1/exp(1) + 1e-8 && xi < -1/exp(1)*0.5) {
        p <- sqrt(2*(exp(1)*xi + 1))
        w <- -1 + p - p^2/3
      }
    } else {
      L1 <- log(xi); L2 <- log(L1)
      w <- L1 - L2 + L2/L1
    }
    for (k in 1:100) {
      ew <- exp(w)
      f  <- w*ew - xi
      wp1 <- w + 1
      denom <- ew*wp1 - ((w+2)*f)/(2*wp1)
      if (denom == 0 || !is.finite(denom)) break
      wnew <- w - f/denom
      if (abs(wnew - w) < 1e-14) { w <- wnew; break }
      w <- wnew
    }
    out[i] <- w
  }
  out
}

#' Lambert W function, lower real branch (W_{-1})
#'
#' Numerically evaluates the lower (negative) real branch of the Lambert W
#' function on \eqn{[-1/e, 0)}, via Halley's iteration.
#'
#' @param x Numeric vector with \eqn{-1/e \le x < 0}.
#' @return Numeric vector of \eqn{W_{-1}(x)} values (\code{NA} outside the
#'   valid domain).
#' @examples
#' lambertWm1(-0.1)      # about -3.5772
#' lambertWm1(-1/exp(1)) # -1
#' @export
lambertWm1 <- function(x) {
  # lower (negative) real branch, valid for -1/e <= x < 0
  out <- numeric(length(x))
  for (i in seq_along(x)) {
    xi <- x[i]
    if (is.na(xi)) { out[i] <- NA; next }
    if (xi >= 0) { out[i] <- NA; next }
    if (xi < -exp(-1) - 1e-12) { out[i] <- NA; next }
    # initial guess near the branch point
    if (xi > -0.3) {
      p <- -sqrt(2*(exp(1)*xi + 1))
      w <- -1 + p - p^2/3 - 11*p^3/72
    } else {
      L1 <- log(-xi); L2 <- log(-L1)
      w <- L1 - L2 + L2/L1
    }
    for (k in 1:200) {
      ew <- exp(w)
      f  <- w*ew - xi
      wp1 <- w + 1
      denom <- ew*wp1 - ((w+2)*f)/(2*wp1)
      if (denom == 0 || !is.finite(denom)) break
      wnew <- w - f/denom
      if (!is.finite(wnew)) break
      if (abs(wnew - w) < 1e-14) { w <- wnew; break }
      w <- wnew
    }
    out[i] <- w
  }
  out
}
