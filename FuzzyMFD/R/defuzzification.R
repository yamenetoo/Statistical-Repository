## Defuzzification (Ross 2000, Section 5) ------------------------------------

defuzzify <- function(sets, method = c("centroid", "max", "meanmax", "weighted"),
                      ...) {
  ## `sets` is one fuzzy set or a list of fuzzy sets (the parts of a fuzzy
  ## output, Equation 14). The centroid, max-membership and mean-max methods
  ## work on the union of the parts; the weighted-average method on the parts.
  method <- match.arg(method)
  if (inherits(sets, "fuzzyset")) sets <- list(sets)
  sets <- .as_set_list(sets)
  if (method == "weighted") return(defuzz_weighted(sets, ...))
  C <- if (length(sets) == 1L) sets[[1L]] else fuzzy_union(sets)
  switch(method,
         centroid = defuzz_centroid(C),
         max = defuzz_max(C, ...),
         meanmax = defuzz_meanmax(C))
}

defuzz_max <- function(A, which = c("first", "last", "all")) {
  ## Equation (15), the height method: z* with mu(z*) >= mu(z) for all z
  which <- match.arg(which)
  .check_fs(A)
  M <- .maximisers(A)
  switch(which,
         first = M[1, 1],
         last = M[nrow(M), 2],
         all = M)
}

defuzz_centroid <- function(A) {
  ## Equation (16), centre of area: int mu(z) z dz / int mu(z) dz, computed
  ## exactly for the piecewise-linear membership function
  .check_fs(A)
  x <- A$x
  m <- A$mu
  n <- length(x)
  if (n == 1L) return(x)
  h <- diff(x)
  m0 <- m[-n]
  m1 <- m[-1L]
  x0 <- x[-n]
  x1 <- x[-1L]
  area <- sum(h * (m0 + m1) / 2)
  if (area <= 0) {
    warning("the membership function has zero area; returning the mean-max value")
    return(defuzz_meanmax(A))
  }
  moment <- sum(h / 6 * (m0 * (2 * x0 + x1) + m1 * (x0 + 2 * x1)))
  moment / area
}

defuzz_weighted <- function(sets, centers = NULL, weights = NULL) {
  ## Equation (17): z* = sum mu(zbar) zbar / sum mu(zbar), each part weighted
  ## by its maximum membership; zbar is the centre of each (symmetric) part,
  ## taken by default as the middle of its maximum-membership plateau
  sets <- .as_set_list(sets)
  if (is.null(centers)) centers <- vapply(sets, defuzz_meanmax, numeric(1))
  if (is.null(weights)) weights <- vapply(sets, fuzzy_height, numeric(1))
  if (length(centers) != length(sets) || length(weights) != length(sets))
    stop("'centers' and 'weights' must have one value per set")
  sym <- vapply(sets, .is_symmetric, logical(1))
  if (!all(sym))
    warning("the weighted-average method is meant for symmetric membership functions")
  sum(weights * centers) / sum(weights)
}

defuzz_meanmax <- function(A) {
  ## Equation (18), middle of maxima: (a + b) / 2 where [a, b] spans the
  ## locations of maximum membership
  .check_fs(A)
  M <- .maximisers(A)
  (M[1, 1] + M[nrow(M), 2]) / 2
}

## internal -----------------------------------------------------------------

.maximisers <- function(A, tol = 1e-9) {
  h <- max(A$mu)
  .level_intervals(A, function(m) m >= h - tol, both = TRUE)
}

.is_symmetric <- function(A, tol = 1e-6) {
  c0 <- defuzz_meanmax(A)
  r <- max(abs(A$x - c0))
  if (r == 0) return(TRUE)
  z <- c0 + seq(0, r, length.out = 51)
  all(abs(membership(A, z) - membership(A, 2 * c0 - z)) <= tol)
}
