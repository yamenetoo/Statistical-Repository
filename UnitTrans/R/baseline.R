# Helpers to call baseline d/p/q/r functions generically

.location_names <- c("mean", "mu", "location", "m", "xi", "meanlog", "center", "mode")

.get_fun <- function(prefix, dist) {
  fn <- paste0(prefix, dist)
  f <- tryCatch(get(fn, mode = "function"),
                error = function(e) NULL)
  if (is.null(f)) {
    if (requireNamespace("extraDistr", quietly = TRUE) &&
        exists(fn, envir = asNamespace("extraDistr"), inherits = FALSE))
      f <- get(fn, envir = asNamespace("extraDistr"))
    else if (exists(fn, envir = asNamespace("stats"), inherits = FALSE))
      f <- get(fn, envir = asNamespace("stats"))
  }
  if (is.null(f)) stop("Cannot find function '", fn, "'.", call. = FALSE)
  f
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
baseline_params <- function(dist) {
  f <- .get_fun("d", dist)
  nm <- names(formals(f))[-1]
  nm <- setdiff(nm, c("log", "lower.tail", "log.p", "ncp", "..."))
  # distributions such as gamma accept rate OR scale: keep rate only
  if (all(c("rate", "scale") %in% nm)) nm <- setdiff(nm, "scale")
  nm
}

.name_par <- function(par, dist) {
  if (is.list(par)) par <- unlist(par)
  if (is.null(names(par))) {
    nm <- baseline_params(dist)
    if (length(par) > length(nm))
      stop("Too many parameters for baseline '", dist, "'.", call. = FALSE)
    names(par) <- nm[seq_along(par)]
  }
  par
}

.call_base <- function(prefix, dist, x, par, ...) {
  f <- .get_fun(prefix, dist)
  par <- .name_par(par, dist)
  a <- c(list(x), as.list(par), list(...))
  names(a)[1] <- names(formals(f))[1]
  do.call(f, a)
}

#' Check that a baseline/transformation pair is valid
#'
#' Evaluates the baseline quantile function on a fine grid of probabilities to
#' verify that the support of Y is compatible with the transformation, and
#' reports the range of V actually reached.
#' @param dist Baseline distribution name.
#' @param par Baseline parameters (named vector/list or unnamed in formals order).
#' @param trans Transformation name or object.
#' @return A list with `valid` (logical), `range` (approximate support of V)
#'   and `message`.
#' @examples
#' check_domain("norm", c(0, 1), "logit")
#' check_domain("norm", c(0, 1), "rational")
#' @export
check_domain <- function(dist, par, trans) {
  tr <- get_transformation(trans)
  p <- c(1e-6, seq(0.001, 0.999, length.out = 200), 1 - 1e-6)
  y <- suppressWarnings(.call_base("q", dist, p, par))
  if (any(!is.finite(y)) && !all(is.finite(y[-c(1, length(y))])))
    return(list(valid = FALSE, range = c(NA, NA), message = "baseline quantiles not finite"))
  ok <- switch(tr$domain,
               real = TRUE,
               nonneg = all(y[is.finite(y)] >= 0),
               positive = all(y[is.finite(y)] > 0))
  if (!ok)
    return(list(valid = FALSE, range = c(NA, NA),
                message = paste0("baseline support incompatible with '", tr$name,
                                 "' (needs ", tr$domain, " support)")))
  v <- tr$g(y[is.finite(y)])
  rng <- range(v)
  msg <- if (diff(rng) < 0.5)
    "valid, but V covers only a small part of (0,1)" else "valid"
  list(valid = TRUE, range = rng, message = msg)
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
  nm <- names(.registry())
  res <- lapply(nm, function(n) {
    d <- check_domain(dist, par, n)
    data.frame(transformation = n, valid = d$valid,
               v_min = d$range[1], v_max = d$range[2], note = d$message,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, res)
}
