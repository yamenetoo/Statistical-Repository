## Features of the membership function (Ross 2000, Section 2.1) -----------

fuzzy_height <- function(A) {
  .check_fs(A)
  max(A$mu)
}

is_normal <- function(A, tol = 1e-9) {
  fuzzy_height(A) >= 1 - tol
}

is_subnormal <- function(A, tol = 1e-9) {
  !is_normal(A, tol = tol)
}

is_convex <- function(A, tol = 1e-9) {
  ## Equation (3): mu(y) >= min(mu(x), mu(z)) for all x < y < z.
  ## For a piecewise-linear membership function it suffices to check the
  ## knots: no knot may lie below both the highest knot to its left and the
  ## highest knot to its right.
  .check_fs(A)
  mu <- A$mu
  n <- length(mu)
  if (n < 3L) return(TRUE)
  left <- cummax(mu)
  right <- rev(cummax(rev(mu)))
  all(mu >= pmin(left, right) - tol)
}

is_fuzzy_number <- function(A, tol = 1e-9) {
  is_normal(A, tol = tol) && is_convex(A, tol = tol)
}

fuzzy_core <- function(A, tol = 1e-9) {
  ## region of full membership, mu(x) = 1
  .check_fs(A)
  .level_intervals(A, function(m) m >= 1 - tol, both = TRUE)
}

fuzzy_support <- function(A, tol = 1e-9) {
  ## (closure of the) region of non-zero membership, mu(x) > 0
  .check_fs(A)
  .level_intervals(A, function(m) m > tol, both = FALSE)
}

fuzzy_boundary <- function(A, tol = 1e-9) {
  ## region of partial membership, 0 < mu(x) < 1
  .check_fs(A)
  x <- A$x
  mu <- A$mu
  n <- length(x)
  if (n < 2L) return(.intervals_matrix(numeric(0), numeric(0)))
  lo <- hi <- numeric(0)
  for (i in seq_len(n - 1L)) {
    a <- mu[i]
    b <- mu[i + 1L]
    zero <- a <= tol && b <= tol
    full <- a >= 1 - tol && b >= 1 - tol
    if (!zero && !full) {
      lo <- c(lo, x[i])
      hi <- c(hi, x[i + 1L])
    }
  }
  .merge_intervals(lo, hi)
}

crossover_points <- function(A, level = 0.5, tol = 1e-9) {
  ## elements with mu(x) = level (0.5 by default)
  .check_fs(A)
  x <- A$x
  d <- A$mu - level
  out <- x[abs(d) <= tol]
  n <- length(x)
  if (n >= 2L) {
    for (i in seq_len(n - 1L)) {
      if (d[i] * d[i + 1L] < 0) {
        out <- c(out, x[i] + (x[i + 1L] - x[i]) * d[i] / (d[i] - d[i + 1L]))
      }
    }
  }
  sort(unique(round(out, 12)))
}

fuzzy_prototype <- function(A, tol = 1e-9) {
  ## the single element with full membership, norm(A), if it exists
  core <- fuzzy_core(A, tol = tol)
  if (nrow(core) == 1L && abs(core[1, 2] - core[1, 1]) <= tol) core[1, 1] else NA_real_
}

## internal helpers ---------------------------------------------------------

.intervals_matrix <- function(lo, hi) {
  m <- cbind(lower = lo, upper = hi)
  if (!length(lo)) m <- matrix(numeric(0), 0, 2, dimnames = list(NULL, c("lower", "upper")))
  m
}

.merge_intervals <- function(lo, hi) {
  if (!length(lo)) return(.intervals_matrix(numeric(0), numeric(0)))
  o <- order(lo)
  lo <- lo[o]
  hi <- hi[o]
  L <- lo[1L]
  H <- hi[1L]
  outL <- outH <- numeric(0)
  if (length(lo) > 1L) {
    for (i in 2:length(lo)) {
      if (lo[i] <= H + 1e-12) {
        H <- max(H, hi[i])
      } else {
        outL <- c(outL, L)
        outH <- c(outH, H)
        L <- lo[i]
        H <- hi[i]
      }
    }
  }
  .intervals_matrix(c(outL, L), c(outH, H))
}

## intervals where cond(mu) holds; `both = TRUE` requires the condition at
## both ends of a segment (used for the core), `FALSE` at either end
.level_intervals <- function(A, cond, both) {
  x <- A$x
  ok <- cond(A$mu)
  n <- length(x)
  lo <- hi <- numeric(0)
  covered <- rep(FALSE, n)
  if (n >= 2L) {
    for (i in seq_len(n - 1L)) {
      seg <- if (both) ok[i] && ok[i + 1L] else ok[i] || ok[i + 1L]
      if (seg) {
        lo <- c(lo, x[i])
        hi <- c(hi, x[i + 1L])
        covered[c(i, i + 1L)] <- TRUE
      }
    }
  }
  single <- which(ok & !covered)
  lo <- c(lo, x[single])
  hi <- c(hi, x[single])
  .merge_intervals(lo, hi)
}
