# Internal helpers: resolving a distribution, calling its functions, fitting it.

#' Resolve a distribution specification
#'
#' @param distr Either a single character string such as `"weibull"` (the
#'   functions `dweibull`, `qweibull` and `rweibull` are then looked up), or a
#'   named list with functions `d`, `q` and `r`.
#' @param env Environment in which to look up named functions.
#' @return A list with elements `name`, `d`, `q`, `r`.
#' @noRd
.resolve_distr <- function(distr, env = parent.frame()) {
  if (is.character(distr) && length(distr) == 1L) {
    find <- function(prefix) {
      nm <- paste0(prefix, distr)
      f <- get0(nm, envir = env, mode = "function", inherits = TRUE)
      if (is.null(f)) {
        stop(sprintf("Could not find function '%s' for distribution '%s'.",
                     nm, distr), call. = FALSE)
      }
      f
    }
    return(list(name = distr, d = find("d"), q = find("q"), r = find("r")))
  }
  if (is.list(distr) && all(c("d", "q", "r") %in% names(distr)) &&
      all(vapply(distr[c("d", "q", "r")], is.function, logical(1)))) {
    nm <- if (!is.null(distr$name)) as.character(distr$name) else "custom"
    return(list(name = nm, d = distr$d, q = distr$q, r = distr$r))
  }
  stop("'distr' must be a distribution name (e.g. \"weibull\") or a list ",
       "with functions d, q and r.", call. = FALSE)
}

# Call f(x, <named parameters>)
.call_fun <- function(f, x, par) do.call(f, c(list(x), par))

# Check that params is a named list of finite numbers
.check_params <- function(params, arg = "params") {
  params <- as.list(params)
  if (length(params) == 0L || is.null(names(params)) || any(names(params) == "")) {
    stop(sprintf("'%s' must be a fully named list, e.g. list(shape = 2, scale = 1).",
                 arg), call. = FALSE)
  }
  if (!all(vapply(params, function(p) is.numeric(p) && length(p) == 1L && is.finite(p),
                  logical(1)))) {
    stop(sprintf("All entries of '%s' must be single finite numbers.", arg),
         call. = FALSE)
  }
  params
}

#' Build a maximum likelihood fitter
#'
#' Returns a function of a sample `x` giving a named list of estimates, or
#' `NULL` when the optimisation fails.
#' @noRd
.make_mle_fitter <- function(d, start, method, lower, upper) {
  nms <- names(start)
  start_vec <- unlist(start)
  use_bounds <- method %in% c("L-BFGS-B", "Brent")
  function(x) {
    nll <- function(par) {
      pl <- stats::setNames(as.list(par), nms)
      dens <- tryCatch(suppressWarnings(.call_fun(d, x, pl)),
                       error = function(e) NULL)
      if (is.null(dens) || length(dens) != length(x) ||
          any(!is.finite(dens)) || any(dens <= 0)) {
        return(1e10)
      }
      -sum(log(dens))
    }
    args <- list(par = start_vec, fn = nll, method = method)
    if (use_bounds) {
      args$lower <- lower
      args$upper <- upper
    }
    fit <- tryCatch(suppressWarnings(do.call(stats::optim, args)),
                    error = function(e) NULL)
    if (is.null(fit) || fit$convergence != 0 || !is.finite(fit$value) ||
        fit$value >= 1e10) {
      return(NULL)
    }
    stats::setNames(as.list(fit$par), nms)
  }
}
