## Fuzzy sets on the real line ----------------------------------------------
##
## A fuzzy set is stored as a piecewise-linear membership function: a
## strictly increasing vector of knots `x` and the membership values `mu` at
## the knots. Between knots the membership is interpolated linearly and
## outside [min(x), max(x)] it is zero. Triangles, trapezoids, crisp intervals
## and every shape used by Ross (2000) are represented exactly; smooth shapes
## (Gaussian, bell, sigmoid) are represented on a fine grid.

.tol <- 1e-9

fuzzy_set <- function(x, mu, name = NULL) {
  x <- as.numeric(x)
  mu <- as.numeric(mu)
  if (length(x) != length(mu))
    stop("'x' and 'mu' must have the same length")
  if (length(x) < 1L)
    stop("a fuzzy set needs at least one knot")
  if (anyNA(x) || anyNA(mu) || any(!is.finite(x)))
    stop("'x' and 'mu' must be finite and free of NA")
  if (any(mu < -.tol | mu > 1 + .tol))
    stop("membership values must lie in [0, 1]")
  mu <- pmin(pmax(mu, 0), 1)
  o <- order(x)
  x <- x[o]
  mu <- mu[o]
  if (any(diff(x) <= 0))
    stop("knots 'x' must be distinct")
  structure(list(x = x, mu = mu, name = if (is.null(name)) "" else as.character(name)),
            class = "fuzzyset")
}

fuzzy_set_fun <- function(fun, lower, upper, n = 1001L, name = NULL) {
  if (!is.function(fun)) stop("'fun' must be a function")
  if (!(is.finite(lower) && is.finite(upper) && lower < upper))
    stop("need finite 'lower' < 'upper'")
  x <- seq(lower, upper, length.out = max(2L, as.integer(n)))
  mu <- fun(x)
  if (length(mu) != length(x)) stop("'fun' must be vectorised")
  fuzzy_set(x, mu, name = name)
}

fs_triangular <- function(a, b, c, height = 1, name = NULL) {
  if (!(a <= b && b <= c && a < c))
    stop("need a <= b <= c with a < c")
  .check_height(height)
  x <- c(a, b, c)
  mu <- c(0, height, 0)
  keep <- c(a < b, TRUE, b < c)
  fuzzy_set(x[keep], mu[keep], name = name)
}

fs_trapezoidal <- function(a, b, c, d, height = 1, name = NULL) {
  if (!(a <= b && b <= c && c <= d && a < d))
    stop("need a <= b <= c <= d with a < d")
  .check_height(height)
  x <- c(a, b, c, d)
  mu <- c(0, height, height, 0)
  keep <- c(a < b, TRUE, b < c, c < d)
  fuzzy_set(x[keep], mu[keep], name = name)
}

fs_rectangular <- function(a, b, name = NULL) {
  if (!(a < b)) stop("need a < b")
  fuzzy_set(c(a, b), c(1, 1), name = name)
}

fs_singleton <- function(x0, height = 1, name = NULL) {
  .check_height(height)
  fuzzy_set(x0, height, name = name)
}

fs_gaussian <- function(mean, sd, lower = mean - 4 * sd, upper = mean + 4 * sd,
                        n = 401L, name = NULL) {
  if (!(sd > 0)) stop("'sd' must be positive")
  fuzzy_set_fun(function(z) mf_gaussian(z, mean, sd), lower, upper, n = n,
                name = name)
}

membership <- function(A, z) {
  .check_fs(A)
  z <- as.numeric(z)
  if (length(A$x) == 1L)
    return(ifelse(abs(z - A$x) <= .tol, A$mu, 0))
  stats::approx(A$x, A$mu, xout = z, yleft = 0, yright = 0,
                ties = "ordered")$y
}

## internal helpers ---------------------------------------------------------

.check_fs <- function(A) {
  if (!inherits(A, "fuzzyset")) stop("expected an object of class 'fuzzyset'")
  invisible(TRUE)
}

.check_height <- function(h) {
  if (!(is.numeric(h) && length(h) == 1L && h > 0 && h <= 1))
    stop("'height' must be a single number in (0, 1]")
  invisible(TRUE)
}

.as_set_list <- function(...) {
  sets <- list(...)
  if (length(sets) == 1L && is.list(sets[[1L]]) && !inherits(sets[[1L]], "fuzzyset"))
    sets <- sets[[1L]]
  ok <- vapply(sets, inherits, logical(1), what = "fuzzyset")
  if (!length(sets) || !all(ok)) stop("all arguments must be 'fuzzyset' objects")
  sets
}

## S3 methods ---------------------------------------------------------------

print.fuzzyset <- function(x, digits = 4, ...) {
  nm <- if (nzchar(x$name)) paste0(" '", x$name, "'") else ""
  cat("Fuzzy set", nm, " (piecewise linear, ", length(x$x), " knots)\n", sep = "")
  n <- length(x$x)
  show <- if (n > 12L) c(1:6, (n - 5L):n) else seq_len(n)
  tab <- cbind(x = format(signif(x$x[show], digits)),
               mu = format(signif(x$mu[show], digits)))
  rownames(tab) <- show
  if (n > 12L) {
    tab <- rbind(tab[1:6, , drop = FALSE], c("...", "..."), tab[7:12, , drop = FALSE])
    rownames(tab) <- c(1:6, "", (n - 5L):n)
  }
  print(tab, quote = FALSE, right = TRUE)
  cat("height =", signif(fuzzy_height(x), digits),
      if (is_normal(x)) "(normal)" else "(sub-normal)",
      "| convex =", is_convex(x), "\n")
  invisible(x)
}

summary.fuzzyset <- function(object, ...) {
  out <- list(name = object$name,
              height = fuzzy_height(object),
              normal = is_normal(object),
              convex = is_convex(object),
              fuzzy_number = is_fuzzy_number(object),
              core = fuzzy_core(object),
              support = fuzzy_support(object),
              boundary = fuzzy_boundary(object),
              crossover = crossover_points(object),
              prototype = fuzzy_prototype(object),
              centroid = defuzz_centroid(object))
  class(out) <- "summary.fuzzyset"
  out
}

print.summary.fuzzyset <- function(x, digits = 4, ...) {
  iv <- function(m) {
    if (!nrow(m)) return("(empty)")
    paste(sprintf("[%s, %s]", signif(m[, 1], digits), signif(m[, 2], digits)),
          collapse = " U ")
  }
  cat("Fuzzy set", if (nzchar(x$name)) paste0("'", x$name, "'") else "", "\n")
  cat("  height       :", signif(x$height, digits), "\n")
  cat("  normal       :", x$normal, "\n")
  cat("  convex       :", x$convex, "\n")
  cat("  fuzzy number :", x$fuzzy_number, "\n")
  cat("  core         :", iv(x$core), "\n")
  cat("  support      :", iv(x$support), "\n")
  cat("  boundaries   :", iv(x$boundary), "\n")
  cat("  crossover    :", if (length(x$crossover)) signif(x$crossover, digits) else "(none)", "\n")
  cat("  prototype    :", if (is.na(x$prototype)) "(none)" else signif(x$prototype, digits), "\n")
  cat("  centroid     :", signif(x$centroid, digits), "\n")
  invisible(x)
}
