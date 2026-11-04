# Helpers to call baseline d/p/q/r functions generically

.location_names <- c("mean", "mu", "location", "m", "xi", "meanlog", "center", "mode")

# Locate d/p/q/r<dist>: global/search path first, then extraDistr, then stats.
# Returns NULL when nothing is found.
.find_fun <- function(prefix, dist) {
  fn <- paste0(prefix, dist)
  f <- tryCatch(get(fn, mode = "function"), error = function(e) NULL)
  if (is.null(f)) {
    if (requireNamespace("extraDistr", quietly = TRUE) &&
        exists(fn, envir = asNamespace("extraDistr"), inherits = FALSE))
      f <- get(fn, envir = asNamespace("extraDistr"))
    else if (exists(fn, envir = asNamespace("stats"), inherits = FALSE))
      f <- get(fn, envir = asNamespace("stats"))
  }
  f
}

.get_fun <- function(prefix, dist) {
  f <- .find_fun(prefix, dist)
  if (is.null(f)) stop("Cannot find function '", prefix, dist, "'.", call. = FALSE)
  f
}

.params_from <- function(f) {
  nm <- names(formals(f))[-1]
  nm <- setdiff(nm, c("log", "lower.tail", "log.p", "ncp", "..."))
  # distributions such as gamma accept rate OR scale: keep rate only
  if (all(c("rate", "scale") %in% nm)) nm <- setdiff(nm, "scale")
  nm
}

# A resolved baseline: the d/p/q functions are looked up once. Every exported
# function accepts either a name or such an object, so that the optimiser and
# the bootstrap loops never repeat the lookup.
.as_base <- function(dist) {
  if (inherits(dist, "unitbase")) return(dist)
  if (!is.character(dist) || length(dist) != 1L || is.na(dist))
    stop("'dist' must be the name of a baseline distribution.", call. = FALSE)
  d <- .get_fun("d", dist)
  structure(list(name = dist, d = d, p = .find_fun("p", dist),
                 q = .find_fun("q", dist), params = .params_from(d)),
            class = "unitbase")
}

.accepts <- function(f, arg) arg %in% names(formals(f))

.invoke <- function(f, x, par, extra = list()) {
  a <- c(list(x), as.list(par), extra)
  names(a)[1] <- names(formals(f))[1]
  do.call(f, a)
}

.need <- function(bl, what) {
  f <- bl[[what]]
  if (is.null(f))
    stop("Cannot find function '", what, bl$name, "' (needed for this operation).",
         call. = FALSE)
  f
}

# log-density of the baseline (falls back to log(d) when `log` is unsupported)
.call_d <- function(bl, y, par) {
  f <- bl$d
  if (.accepts(f, "log")) .invoke(f, y, par, list(log = TRUE))
  else log(.invoke(f, y, par))
}

# baseline cdf; handles lower.tail / log.p even if the function lacks them
.call_p <- function(bl, y, par, lower.tail = TRUE, log.p = FALSE) {
  f <- .need(bl, "p")
  has_lt <- .accepts(f, "lower.tail"); has_lg <- .accepts(f, "log.p")
  extra <- list()
  if (has_lt) extra$lower.tail <- lower.tail
  if (has_lg) extra$log.p <- log.p
  val <- .invoke(f, y, par, extra)
  if (!has_lt && !lower.tail)
    val <- if (has_lg && log.p) log1p(-exp(val)) else 1 - val
  if (!has_lg && log.p) val <- log(val)
  val
}

# baseline quantile function; handles lower.tail / log.p even if unsupported
.call_q <- function(bl, p, par, lower.tail = TRUE, log.p = FALSE) {
  f <- .need(bl, "q")
  has_lt <- .accepts(f, "lower.tail"); has_lg <- .accepts(f, "log.p")
  if (!has_lg && log.p) { p <- exp(p); log.p <- FALSE }
  if (!has_lt && !lower.tail) {
    p <- if (log.p) log1p(-exp(p)) else 1 - p
    lower.tail <- TRUE
  }
  extra <- list()
  if (has_lt) extra$lower.tail <- lower.tail
  if (has_lg) extra$log.p <- log.p
  .invoke(f, p, par, extra)
}

#' Parameter names of a baseline distribution
#'
#' Reads the parameter names from the formals of the `d<dist>` function,
#' dropping the first argument and arguments such as `log`.
#' @param dist Name of the baseline distribution, e.g. `"gamma"`.
#' @return Character vector of parameter names.
#' @examples
#' baseline_params("gamma")
#' @export
baseline_params <- function(dist) .as_base(dist)$params

.name_par <- function(par, dist) {
  if (is.list(par)) par <- unlist(par)
  if (is.null(names(par))) {
    nm <- .as_base(dist)$params
    if (length(par) > length(nm))
      stop("Too many parameters for baseline '", .as_base(dist)$name, "'.",
           call. = FALSE)
    names(par) <- nm[seq_along(par)]
  }
  par
}

#' Check that a baseline/transformation pair is valid
#'
#' Verifies that the support of Y is compatible with the transformation. When
#' the baseline cdf is available the probability mass that Y places outside
#' the admissible domain of the transformation is computed exactly
#' (for example `P(Y < 0)` for `"rational"`); otherwise a grid of quantiles is
#' inspected. The range of V actually reached is also reported.
#' @param dist Baseline distribution name.
#' @param par Baseline parameters (named vector/list or unnamed in formals order).
#' @param trans Transformation name or object.
#' @param tol Largest tolerated probability mass outside the admissible domain.
#' @return A list with `valid` (logical), `range` (approximate support of V),
#'   `mass_outside` and `message`.
#' @examples
#' check_domain("norm", c(0, 1), "logit")
#' check_domain("norm", c(0, 1), "rational")
#' @export
check_domain <- function(dist, par, trans, tol = 1e-8) {
  bl <- .as_base(dist)
  tr <- get_transformation(trans)
  par <- .name_par(par, bl)
  bad <- function(msg, mass = NA_real_)
    list(valid = FALSE, range = c(NA_real_, NA_real_), mass_outside = mass,
         message = msg)

  mass <- 0
  if (tr$domain != "real") {
    if (!is.null(bl$p)) {
      mass <- tryCatch(as.numeric(.call_p(bl, 0, par)), error = function(e) NA_real_)
      if (length(mass) != 1L || is.na(mass))
        return(bad("cannot evaluate the baseline cdf at 0"))
    }
  }
  p <- c(1e-6, seq(0.001, 0.999, length.out = 200), 1 - 1e-6)
  y <- if (!is.null(bl$q))
    suppressWarnings(tryCatch(.call_q(bl, p, par), error = function(e) NULL)) else NULL
  if (is.null(y) && is.null(bl$p))
    return(bad("baseline has neither a cdf nor a quantile function"))
  if (!is.null(y)) {
    yf <- y[is.finite(y)]
    if (length(yf) < 3L) return(bad("baseline quantiles not finite"))
    if (tr$domain != "real" && is.null(bl$p)) {
      mass <- if (all(yf >= 0)) 0 else NA_real_
      if (is.na(mass)) return(bad(paste0("baseline support incompatible with '",
                                         tr$name, "' (needs ", tr$domain, " support)")))
    }
  }
  if (mass > tol)
    return(bad(sprintf("baseline puts mass %.3g outside the domain of '%s' (needs %s support)",
                       mass, tr$name, tr$domain), mass))
  rng <- if (!is.null(y)) range(tr$g(yf)) else c(NA_real_, NA_real_)
  msg <- if (!anyNA(rng) && diff(rng) < 0.5)
    "valid, but V covers only a small part of (0,1)" else "valid"
  list(valid = TRUE, range = rng, mass_outside = mass, message = msg)
}

#' Detect which transformations suit a baseline distribution
#'
#' @param dist Baseline distribution name.
#' @param par Baseline parameters.
#' @return A data frame listing each transformation, whether it is valid and
#'   the approximate range of V.
#' @examples
#' detect_transformations("gamma", c(shape = 3, rate = 1))
#' @export
detect_transformations <- function(dist, par) {
  bl <- .as_base(dist)
  nm <- names(.registry())
  res <- lapply(nm, function(n) {
    d <- check_domain(bl, par, n)
    data.frame(transformation = n, valid = d$valid,
               v_min = d$range[1], v_max = d$range[2], note = d$message,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, res)
}
