## Internal numerical helpers ------------------------------------------------

## log(1 - exp(-x)) for x >= 0, accurate over the whole range (Maechler, 2012)
.log1mexp <- function(x) {
  out <- x
  i <- !is.na(x) & x <= log(2)
  out[i] <- log(-expm1(-x[i]))
  out[!i] <- log1p(-exp(-x[!i]))
  out
}

## log(1 + exp(x)) without overflow
.log1pexp <- function(x) {
  out <- x
  i <- !is.na(x) & x <= 35
  out[i] <- log1p(exp(x[i]))
  out
}

## k * v with the convention 0 * (+-Inf) = 0 (k is a scalar)
.mul <- function(k, v) {
  if (isTRUE(k == 0)) rep(0, length(v)) else k * v
}

## Turn a log-survival into the requested probability scale
.from_logS <- function(lS, lower.tail = TRUE, log.p = FALSE) {
  if (!lower.tail) return(if (log.p) lS else exp(lS))
  lF <- .log1mexp(-lS)
  if (log.p) lF else exp(lF)
}

## Turn a log-cdf into the requested probability scale
.from_logF <- function(lF, lower.tail = TRUE, log.p = FALSE) {
  if (lower.tail) return(if (log.p) lF else exp(lF))
  lS <- .log1mexp(-lF)
  if (log.p) lS else exp(lS)
}

## Shape parameters ------------------------------------------------------------

.shape_names <- c("a", "b", "alpha", "beta")

.shape_ok <- function(a, b, alpha, beta) {
  v <- c(a, b, alpha, beta)
  length(v) == 4L && all(is.finite(v)) && all(v > 0)
}

.check_shape <- function(a, b, alpha, beta) {
  if (length(a) != 1L || length(b) != 1L || length(alpha) != 1L ||
      length(beta) != 1L)
    stop("'a', 'b', 'alpha' and 'beta' must be scalars", call. = FALSE)
  if (!is.numeric(c(a, b, alpha, beta)) || !.shape_ok(a, b, alpha, beta))
    stop("'a', 'b', 'alpha' and 'beta' must be positive and finite",
         call. = FALSE)
  invisible(TRUE)
}

## Baseline parameters: fill defaults, order, validate
.check_par <- function(par, B) {
  if (length(B$par) == 0L) return(numeric(0))
  if (is.null(par)) return(B$par)
  par <- unlist(par)
  if (!is.numeric(par))
    stop("'par' must be a numeric vector or a list of numbers", call. = FALSE)
  if (is.null(names(par))) {
    if (length(par) != length(B$par))
      stop("unnamed 'par' must have length ", length(B$par), " (",
           paste(names(B$par), collapse = ", "), ")", call. = FALSE)
    names(par) <- names(B$par)
  }
  bad <- setdiff(names(par), names(B$par))
  if (length(bad))
    stop("unknown baseline parameter(s) for '", B$name, "': ",
         paste(bad, collapse = ", "), "; expected: ",
         paste(names(B$par), collapse = ", "), call. = FALSE)
  out <- B$par
  out[names(par)] <- par
  if (any(!is.finite(out)) || any(out <= B$lower) || any(out >= B$upper))
    stop("baseline parameter(s) outside their admissible range", call. = FALSE)
  out
}

## Collect and validate everything a distribution function needs
.prep <- function(a, b, alpha, beta, baseline, par) {
  .check_shape(a, b, alpha, beta)
  B <- .get_baseline(baseline)
  list(a = a, b = b, alpha = alpha, beta = beta, B = B,
       par = .check_par(par, B))
}

## Unconstrained reparametrisation used by the optimiser -----------------------

.to_z <- function(p, lo, up) {
  z <- p
  for (k in seq_along(p)) {
    if (is.finite(lo[k]) && is.finite(up[k])) {
      z[k] <- stats::qlogis((p[k] - lo[k]) / (up[k] - lo[k]))
    } else if (is.finite(lo[k])) {
      z[k] <- log(p[k] - lo[k])
    } else if (is.finite(up[k])) {
      z[k] <- log(up[k] - p[k])
    }
  }
  z
}

.from_z <- function(z, lo, up) {
  p <- z
  for (k in seq_along(z)) {
    if (is.finite(lo[k]) && is.finite(up[k])) {
      p[k] <- lo[k] + (up[k] - lo[k]) * stats::plogis(z[k])
    } else if (is.finite(lo[k])) {
      p[k] <- lo[k] + exp(z[k])
    } else if (is.finite(up[k])) {
      p[k] <- up[k] - exp(z[k])
    }
  }
  p
}

## Numerical integration with a stricter tolerance and a safe fallback
.integrate <- function(f, lower, upper, ...) {
  args <- list(...)
  if (is.null(args$subdivisions)) args$subdivisions <- 2000L
  strict <- args
  if (is.null(strict$rel.tol)) strict$rel.tol <- 1e-9
  res <- tryCatch(do.call(stats::integrate,
                          c(list(f = f, lower = lower, upper = upper), strict)),
                  error = function(e) e)
  if (inherits(res, "error")) {
    res <- tryCatch(do.call(stats::integrate,
                            c(list(f = f, lower = lower, upper = upper), args)),
                    error = function(e) e)
  }
  if (inherits(res, "error")) {
    warning("numerical integration failed: ", conditionMessage(res),
            call. = FALSE)
    return(NA_real_)
  }
  res$value
}
