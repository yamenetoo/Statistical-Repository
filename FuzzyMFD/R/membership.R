## Vectorised membership functions -----------------------------------------

mf_indicator <- function(x, lower, upper) {
  ## Equation (1): the indicator (characteristic) function of a crisp set
  as.numeric(x >= lower & x <= upper)
}

mf_triangular <- function(x, a, b, c, height = 1) {
  if (!(a <= b && b <= c && a < c)) stop("need a <= b <= c with a < c")
  y <- numeric(length(x))
  up <- x >= a & x <= b
  dn <- x >= b & x <= c
  if (b > a) y[up] <- (x[up] - a) / (b - a)
  if (c > b) y[dn] <- pmax(y[dn], (c - x[dn]) / (c - b))
  y[x == b] <- 1
  height * y
}

mf_trapezoidal <- function(x, a, b, c, d, height = 1) {
  if (!(a <= b && b <= c && c <= d && a < d))
    stop("need a <= b <= c <= d with a < d")
  y <- numeric(length(x))
  y[x >= b & x <= c] <- 1
  up <- x >= a & x < b
  dn <- x > c & x <= d
  if (b > a) y[up] <- (x[up] - a) / (b - a)
  if (d > c) y[dn] <- (d - x[dn]) / (d - c)
  height * y
}

mf_gaussian <- function(x, mean, sd) {
  if (!(sd > 0)) stop("'sd' must be positive")
  exp(-0.5 * ((x - mean) / sd)^2)
}

mf_bell <- function(x, a, b, c) {
  ## generalised bell: 1 / (1 + |(x - c) / a|^(2b))
  if (!(a > 0 && b > 0)) stop("'a' and 'b' must be positive")
  1 / (1 + abs((x - c) / a)^(2 * b))
}

mf_sigmoid <- function(x, a, c) {
  1 / (1 + exp(-a * (x - c)))
}
