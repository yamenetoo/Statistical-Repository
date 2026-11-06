#' Define a distribution for stress-strength work
#'
#' Builds a distribution object from a family name or from user-supplied
#' functions. When only \code{name} is given, the functions
#' \code{d<name>}, \code{p<name>}, \code{q<name>} and \code{r<name>} are looked
#' up (for example \code{"weibull"} finds \code{dweibull}, \code{pweibull}, ...)
#' and the parameter names are read from the arguments of the density, so any
#' family that follows R's \code{d/p/q/r} convention works without extra code.
#' Only \code{d} and \code{p} are required; missing quantile and random
#' functions are replaced by numerical inversion of \code{p}.
#'
#' @param name Family name, e.g. \code{"weibull"}, \code{"gamma"},
#'   \code{"toppleone"}, \code{"genrayleigh"}.
#' @param d,p,q,r Optional density, distribution, quantile and random
#'   functions. They override those found from \code{name}. The first argument
#'   of each is the point (or n), the remaining named arguments are parameters.
#' @param par Optional named numeric vector of parameter values. Its names
#'   select which arguments of \code{d} are treated as parameters.
#' @param constraints Optional named list giving the range of each parameter:
#'   \code{"real"}, \code{"positive"}, \code{"unit"} or a numeric interval
#'   \code{c(a, b)}. Unlisted parameters are guessed from their names
#'   (location-type names are real, everything else positive).
#' @param lower,upper Support of the distribution. Guessed from the density
#'   when omitted; give them explicitly for bounded supports other than (0, 1).
#' @param envir Environment in which to look up the functions.
#' @return An object of class \code{"ss_dist"}.
#' @examples
#' X <- ss_dist("toppleone")
#' Y <- ss_dist("weibull", par = c(shape = 2, scale = 1))
#' Z <- ss_dist("gamma")            # parameters: shape, rate
#' @export
ss_dist <- function(name = NULL, d = NULL, p = NULL, q = NULL, r = NULL,
                    par = NULL, constraints = NULL, lower = NULL, upper = NULL,
                    envir = parent.frame()) {
  bi <- if (!is.null(name)) .ss_builtin()[[name]] else NULL
  if (!is.null(bi)) {
    if (is.null(d)) d <- bi$d
    if (is.null(p)) p <- bi$p
    if (is.null(q)) q <- bi$q
    if (is.null(r)) r <- bi$r
    if (is.null(lower)) lower <- bi$lower
    if (is.null(upper)) upper <- bi$upper
    constraints <- modifyList(bi$constraints, as.list(constraints))
  } else if (!is.null(name)) {
    if (is.null(d)) d <- .find_fun("d", name, envir)
    if (is.null(p)) p <- .find_fun("p", name, envir)
    if (is.null(q)) q <- .find_fun("q", name, envir)
    if (is.null(r)) r <- .find_fun("r", name, envir)
  }
  if (!is.function(d) || !is.function(p))
    stop("a density 'd' and a distribution function 'p' are required",
         if (!is.null(name)) paste0("; none found for '", name, "'"))

  fm <- formals(d)
  cand <- setdiff(names(fm), .skip_args)
  if (!is.null(par)) {
    if (is.null(names(par))) {
      if (length(par) > length(cand)) stop("too many parameter values")
      names(par) <- cand[seq_along(par)]
    }
    pn <- names(par)
    if (!all(pn %in% names(fm)))
      stop("unknown parameter(s): ", paste(setdiff(pn, names(fm)), collapse = ", "))
  } else {
    pn <- cand
    if (all(c("rate", "scale") %in% pn)) pn <- setdiff(pn, "scale")
  }
  if (length(pn) == 0) stop("no parameters found in the density's arguments")

  cons <- lapply(pn, function(nm) {
    if (!is.null(constraints[[nm]])) .bounds(constraints[[nm]]) else .guess_bounds(nm)
  })
  names(cons) <- pn

  is_num <- vapply(fm, is.numeric, logical(1))
  fdef <- as.list(fm[is_num])
  default <- vapply(pn, function(nm) .default_value(cons[[nm]], fdef[[nm]]), numeric(1))
  names(default) <- pn

  if (!is.null(par)) {
    par <- par[pn]
    .check_inside(par, cons)
  }

  if (is.null(lower) || is.null(upper)) {
    s <- .guess_support(d, if (is.null(par)) default else par)
    if (is.null(lower)) lower <- s[1]
    if (is.null(upper)) upper <- s[2]
    message("support of '", if (is.null(name)) "user distribution" else name,
            "' taken as (", lower, ", ", upper,
            "); pass 'lower'/'upper' if this is wrong")
  }

  structure(list(name = name, d = d, p = p, q = q, r = r, par = par,
                 par_names = pn, constraints = cons, default = default,
                 lower = lower, upper = upper,
                 has_log = "log" %in% names(fm)),
            class = "ss_dist")
}

#' @export
print.ss_dist <- function(x, ...) {
  cat("Distribution:", if (is.null(x$name)) "<user-supplied>" else x$name, "\n")
  cat("Support: (", x$lower, ", ", x$upper, ")\n", sep = "")
  cat("Parameters:\n")
  for (nm in x$par_names) {
    b <- x$constraints[[nm]]
    cat("  ", nm, " in (", b[1], ", ", b[2], ")",
        if (!is.null(x$par)) paste0(" = ", format(x$par[[nm]])), "\n", sep = "")
  }
  invisible(x)
}

# ---------------------------------------------------------------------------
# internal helpers
# ---------------------------------------------------------------------------

.skip_args <- c("x", "q", "p", "n", "log", "log.p", "lower.tail", "...")

.find_fun <- function(prefix, name, envir) {
  fname <- paste0(prefix, name)
  if (exists(fname, envir = envir, mode = "function"))
    get(fname, envir = envir, mode = "function") else NULL
}

.bounds <- function(con) {
  if (is.character(con)) {
    con <- match.arg(con, c("real", "positive", "unit"))
    return(switch(con, real = c(-Inf, Inf), positive = c(0, Inf), unit = c(0, 1)))
  }
  if (!is.numeric(con) || length(con) != 2 || !(con[1] < con[2]))
    stop("a constraint must be 'real', 'positive', 'unit' or c(a, b) with a < b")
  as.numeric(con)
}

.guess_bounds <- function(nm) {
  nm <- tolower(nm)
  if (nm %in% c("mean", "mu", "location", "loc", "m", "ncp", "xi",
                "median", "meanlog")) return(c(-Inf, Inf))
  if (nm %in% c("prob", "pi")) return(c(0, 1))
  c(0, Inf)
}

.default_value <- function(b, fdefault = NULL) {
  if (is.numeric(fdefault) && length(fdefault) == 1 && is.finite(fdefault) &&
      fdefault > b[1] && fdefault < b[2]) return(fdefault)
  if (is.finite(b[1]) && is.finite(b[2])) return(mean(b))
  if (is.finite(b[1])) return(b[1] + 1)
  if (is.finite(b[2])) return(b[2] - 1)
  0
}

.guess_support <- function(d, par) {
  dd <- function(x) tryCatch(suppressWarnings(do.call(d, c(list(x), as.list(par)))),
                             error = function(e) NA_real_)
  lower <- if (isTRUE(any(dd(c(-0.5, -2, -10)) > 0))) -Inf else 0
  upper <- if (isTRUE(any(dd(c(1.5, 3, 10)) > 0))) Inf else 1
  c(lower, upper)
}

.check_inside <- function(par, cons) {
  for (nm in names(par)) {
    b <- cons[[nm]]
    if (!is.finite(par[[nm]]) || par[[nm]] <= b[1] || par[[nm]] >= b[2])
      stop("parameter '", nm, "' = ", par[[nm]], " is outside (", b[1], ", ", b[2], ")")
  }
  invisible(TRUE)
}

.name_par <- function(par, dist) {
  if (is.null(par)) return(NULL)
  par <- unlist(par)
  if (is.null(names(par))) {
    if (length(par) != length(dist$par_names))
      stop("expected ", length(dist$par_names), " parameter value(s)")
    names(par) <- dist$par_names
  }
  if (!all(dist$par_names %in% names(par)))
    stop("missing parameter(s): ",
         paste(setdiff(dist$par_names, names(par)), collapse = ", "))
  par[dist$par_names]
}

.need_par <- function(dist, par) {
  if (is.null(par)) par <- dist$par
  if (is.null(par))
    stop("parameters of '", if (is.null(dist$name)) "distribution" else dist$name,
         "' are not set; pass them explicitly or via ss_dist(par = )")
  .name_par(par, dist)
}

.inside <- function(theta, dist) {
  all(vapply(seq_along(theta), function(j) {
    b <- dist$constraints[[j]]
    is.finite(theta[[j]]) && theta[[j]] > b[1] && theta[[j]] < b[2]
  }, logical(1)))
}

.dens <- function(dist, x, par, log = FALSE) {
  if (log && dist$has_log)
    return(do.call(dist$d, c(list(x), as.list(par), list(log = TRUE))))
  v <- do.call(dist$d, c(list(x), as.list(par)))
  if (log) log(v) else v
}

.cdf <- function(dist, x, par) do.call(dist$p, c(list(x), as.list(par)))

.quant <- function(dist, u, par) {
  if (!is.null(dist$q)) return(do.call(dist$q, c(list(u), as.list(par))))
  vapply(u, function(ui) .num_quantile(dist, ui, par), numeric(1))
}

.rand <- function(dist, n, par) {
  if (!is.null(dist$r)) return(do.call(dist$r, c(list(n), as.list(par))))
  .quant(dist, runif(n), par)
}

.num_quantile <- function(dist, u, par) {
  if (u <= 0) return(dist$lower)
  if (u >= 1) return(dist$upper)
  f <- function(x) .cdf(dist, x, par) - u
  lo <- if (is.finite(dist$lower)) dist$lower else -1
  hi <- if (is.finite(dist$upper)) dist$upper else 1
  if (hi <= lo) hi <- lo + 1
  i <- 0
  while (!is.finite(dist$lower) && f(lo) > 0 && i < 200) { lo <- 2 * lo - 1; i <- i + 1 }
  i <- 0
  while (!is.finite(dist$upper) && f(hi) < 0 && i < 200) { hi <- 2 * hi + 1; i <- i + 1 }
  uniroot(f, c(lo, hi), tol = 1e-12)$root
}

# --- parameter transformations (constrained <-> free scale) ---------------

.tf1 <- function(t, b) {
  a <- b[1]; z <- b[2]
  if (is.finite(a) && is.finite(z)) return(qlogis((t - a) / (z - a)))
  if (is.finite(a)) return(log(t - a))
  if (is.finite(z)) return(log(z - t))
  t
}

.ft1 <- function(e, b) {
  a <- b[1]; z <- b[2]
  if (is.finite(a) && is.finite(z)) return(a + (z - a) * plogis(e))
  if (is.finite(a)) return(a + exp(e))
  if (is.finite(z)) return(z - exp(e))
  e
}

.dft1 <- function(e, b) {           # d theta / d eta (signed)
  a <- b[1]; z <- b[2]
  if (is.finite(a) && is.finite(z)) return((z - a) * plogis(e) * plogis(-e))
  if (is.finite(a)) return(exp(e))
  if (is.finite(z)) return(-exp(e))
  1
}

.ldft1 <- function(e, b) {          # log |d theta / d eta|
  a <- b[1]; z <- b[2]
  if (is.finite(a) && is.finite(z))
    return(log(z - a) + plogis(e, log.p = TRUE) + plogis(-e, log.p = TRUE))
  if (is.finite(a) || is.finite(z)) return(e)
  0
}

.to_free <- function(theta, dist) {
  theta <- .name_par(theta, dist)
  vapply(seq_along(theta), function(j) .tf1(theta[[j]], dist$constraints[[j]]),
         numeric(1))
}

.from_free <- function(eta, dist) {
  setNames(vapply(seq_along(eta), function(j) .ft1(eta[[j]], dist$constraints[[j]]),
                  numeric(1)), dist$par_names)
}

.jac <- function(eta, dist) {
  vapply(seq_along(eta), function(j) .dft1(eta[[j]], dist$constraints[[j]]), numeric(1))
}

.logjac <- function(eta, dist) {
  sum(vapply(seq_along(eta), function(j) .ldft1(eta[[j]], dist$constraints[[j]]),
             numeric(1)))
}

.row <- function(m, i) setNames(as.numeric(m[i, ]), colnames(m))
