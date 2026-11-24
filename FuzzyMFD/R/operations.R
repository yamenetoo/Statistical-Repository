## Set operations on piecewise-linear fuzzy sets -----------------------------
##
## The union (max) and intersection (min) of two piecewise-linear membership
## functions are again piecewise linear, with knots at the knots of both sets
## and at the points where the two functions cross. The operations below add
## these points, so the results are exact (no grid approximation).

fuzzy_union <- function(..., name = NULL) {
  sets <- .as_set_list(...)
  out <- Reduce(function(A, B) .combine(A, B, pmax), sets)
  if (!is.null(name)) out$name <- name
  out
}

fuzzy_intersection <- function(..., name = NULL) {
  sets <- .as_set_list(...)
  out <- Reduce(function(A, B) .combine(A, B, pmin), sets)
  if (!is.null(name)) out$name <- name
  out
}

fuzzy_complement <- function(A, universe = range(A$x), name = NULL) {
  .check_fs(A)
  if (length(universe) != 2L || !(universe[1] < universe[2]))
    stop("'universe' must be c(lower, upper) with lower < upper")
  x <- sort(unique(c(universe, A$x[A$x >= universe[1] & A$x <= universe[2]])))
  out <- fuzzy_set(x, 1 - membership(A, x),
                   name = if (is.null(name)) paste0("not ", A$name) else name)
  out
}

agreement <- function(A, reading) {
  ## degree to which a fuzzy set and a reading agree (Figure 1.7):
  ## the membership of a crisp reading, or the height of the intersection
  ## of the set with a fuzzy reading
  .check_fs(A)
  if (inherits(reading, "fuzzyset"))
    return(fuzzy_height(fuzzy_intersection(A, reading)))
  if (is.numeric(reading)) return(membership(A, reading))
  stop("'reading' must be numeric or a 'fuzzyset'")
}

## internal -----------------------------------------------------------------

.combine <- function(A, B, op) {
  .check_fs(A)
  .check_fs(B)
  span <- diff(range(c(A$x, B$x)))
  eps <- if (span > 0) 1e-9 * span else 1e-9
  ## a set whose end knot has positive membership jumps to zero just outside
  ## it; add a point there so the jump is kept
  extra <- function(S) {
    n <- length(S$x)
    c(if (S$mu[1L] > 0) S$x[1L] - eps, if (S$mu[n] > 0) S$x[n] + eps)
  }
  x <- sort(unique(c(A$x, B$x, extra(A), extra(B))))
  a <- membership(A, x)
  b <- membership(B, x)
  d <- a - b
  cross <- numeric(0)
  if (length(x) >= 2L) {
    for (i in seq_len(length(x) - 1L)) {
      if (d[i] * d[i + 1L] < 0)
        cross <- c(cross, x[i] + (x[i + 1L] - x[i]) * d[i] / (d[i] - d[i + 1L]))
    }
  }
  if (length(cross)) {
    x <- sort(unique(c(x, cross)))
    a <- membership(A, x)
    b <- membership(B, x)
  }
  mu <- op(a, b)
  .simplify(fuzzy_set(x, mu, name = paste(A$name, B$name, sep = "|")))
}

## drop zero knots that are not needed and collinear interior knots
.simplify <- function(S, tol = 1e-12) {
  x <- S$x
  mu <- S$mu
  n <- length(x)
  if (n <= 2L) return(S)
  keep <- rep(TRUE, n)
  for (i in 2:(n - 1L)) {
    ## collinear with neighbours
    s1 <- (mu[i] - mu[i - 1L]) / (x[i] - x[i - 1L])
    s2 <- (mu[i + 1L] - mu[i]) / (x[i + 1L] - x[i])
    if (abs(s1 - s2) <= 1e-10) keep[i] <- FALSE
  }
  ## leading / trailing zeros: keep only the last leading and first trailing
  nz <- which(mu > tol)
  if (length(nz)) {
    first <- nz[1L]
    last <- nz[length(nz)]
    if (first > 2L) keep[seq_len(first - 2L)] <- FALSE
    if (last < n - 1L) keep[(last + 2L):n] <- FALSE
  }
  if (sum(keep) < 1L) keep[1L] <- TRUE
  S$x <- x[keep]
  S$mu <- mu[keep]
  S
}
