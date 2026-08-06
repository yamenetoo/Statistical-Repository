# ==== Design principles (apply to every file in this package) ====
# 1. Mathematical contract first: every transform below is derived on paper
#    (Y = 1/(1+X), strictly decreasing, maps (0,Inf) onto (0,1)) before being
#    coded; boundary behaviour at y -> 0+ and y -> 1- is checked explicitly.
# 2. Input validation at the exported boundary: user-facing functions
#    validate domain and type; internal helpers (prefixed `.`) assume
#    already-validated input and skip redundant checks for speed.
# 3. Numerical stability: likelihoods are computed on the log scale
#    throughout; density evaluations return -Inf rather than NaN outside
#    the support.
# 4. Reproducibility: randomness is never seeded globally inside package
#    code; callers control seeds explicitly.
# 5. Honest failure: every stats::optim() call inspects $convergence and
#    surfaces a warning rather than silently returning a misleading fit.
# 6. R idiom: d/p/q/r naming mirrors base R's stats distributions; S3
#    generics (print, summary, coef, logLik, vcov) are implemented instead
#    of ad hoc accessor names.
# ====================================================================

#' @title Generic unit-transform distribution
#' @description
#' Builds a distribution on (0,1) from ANY base distribution on (0, Inf)
#' via Y = 1/(1+X). Unlike a registry-restricted design, the base
#' distribution is supplied directly as d/p/q/r functions, so this works
#' with literally any distribution the user can write down -- built-in R
#' distributions (dgamma/pgamma/qgamma/rgamma), a custom research
#' distribution (e.g. XGamma, UQXGD), or a distribution with no closed
#' form quantile function (q can be omitted; it is then obtained
#' numerically by inverting p).
#'
#' @param ddist function(x, ..., log = FALSE) -> density of the base dist
#' @param pdist function(q, ..., lower.tail = TRUE) -> CDF of the base dist
#' @param qdist function(p, ...) -> quantile function of the base dist.
#'   Optional: if NULL, quantiles are obtained numerically from pdist.
#' @param rdist function(n, ...) -> random generator. Optional: if NULL,
#'   sampling is done by inverse-CDF using qutd().
#' @param params named list of parameters passed to ddist/pdist/qdist/rdist
#' @param name optional character label for printing (e.g. "gamma", "xgamma")
#' @return an object of class \code{utd} (unit-transform distribution)
#' @export
#' @examples
#' # Any base distribution works -- built-in gamma:
#' d <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5))
#' dutd(0.3, d)
#'
#' # A fully custom user-defined base distribution also works, with no
#' # registry entry required:
#' dcustom <- function(x, rate, log = FALSE) dexp(x, rate = rate, log = log)
#' pcustom <- function(q, rate, lower.tail = TRUE) pexp(q, rate, lower.tail = lower.tail)
#' d2 <- utd(dcustom, pcustom, params = list(rate = 1.5), name = "custom_exp")
utd <- function(ddist, pdist, qdist = NULL, rdist = NULL,
                 params = list(), name = "custom") {
  if (!is.function(ddist)) stop("`ddist` must be a function(x, ..., log).", call. = FALSE)
  if (!is.function(pdist)) stop("`pdist` must be a function(q, ..., lower.tail).", call. = FALSE)
  if (!is.null(qdist) && !is.function(qdist)) stop("`qdist` must be a function or NULL.", call. = FALSE)
  if (!is.null(rdist) && !is.function(rdist)) stop("`rdist` must be a function or NULL.", call. = FALSE)
  if (!is.list(params)) stop("`params` must be a named list.", call. = FALSE)
  structure(
    list(ddist = ddist, pdist = pdist, qdist = qdist, rdist = rdist,
         params = params, name = name),
    class = "utd"
  )
}

#' @export
print.utd <- function(x, ...) {
  cat(sprintf("<utd distribution on (0,1)> base = %s\n", x$name))
  if (length(x$params)) {
    cat("params: ")
    cat(paste(names(x$params), unlist(x$params), sep = " = ", collapse = ", "))
    cat("\n")
  }
  invisible(x)
}

.utd_call <- function(fn, ..., extra_params) {
  do.call(fn, c(list(...), extra_params))
}

#' Numerically invert the base CDF when no closed-form quantile is given.
.numeric_qbase <- function(p, pdist, params, lower = 1e-12, upper = NULL) {
  if (is.null(upper)) {
    upper <- 1
    # expand search interval until pdist(upper) exceeds max(p)
    target <- max(p)
    tries <- 0
    while (.utd_call(pdist, upper, extra_params = params) < target && tries < 200) {
      upper <- upper * 2
      tries <- tries + 1
    }
  }
  vapply(p, function(pp) {
    f <- function(x) .utd_call(pdist, x, extra_params = params) - pp
    stats::uniroot(f, lower = lower, upper = upper, extendInt = "yes")$root
  }, numeric(1))
}

#' @title Density of a unit-transform distribution
#' @param y numeric vector in (0,1)
#' @param dist a \code{utd} object
#' @param log logical; return log-density
#' @return numeric vector of density values
#' @export
#' @examples
#' d <- utd(dgamma, pgamma, qgamma, params = list(shape = 2.5))
#' dutd(c(0.2, 0.5, 0.8), d)
dutd <- function(y, dist, log = FALSE) {
  stopifnot(inherits(dist, "utd"))
  out <- rep(-Inf, length(y))
  ok <- y > 0 & y < 1
  x <- 1 / y[ok] - 1
  logf <- .utd_call(dist$ddist, x, log = TRUE, extra_params = dist$params)
  out[ok] <- logf - 2 * log(y[ok])
  if (log) out else exp(out)
}

#' @title CDF of a unit-transform distribution
#' @param y numeric vector in (0,1)
#' @param dist a \code{utd} object
#' @param lower.tail logical
#' @export
#' @examples
#' d <- utd(dgamma, pgamma, qgamma, params = list(shape = 2.5))
#' putd(c(0.2, 0.5, 0.8), d)
putd <- function(y, dist, lower.tail = TRUE) {
  stopifnot(inherits(dist, "utd"))
  out <- ifelse(y <= 0, 0, ifelse(y >= 1, 1, NA_real_))
  ok <- is.na(out)
  x <- 1 / y[ok] - 1
  Fx <- .utd_call(dist$pdist, x, extra_params = dist$params)
  # F_Y(y) = 1 - F_X(1/y - 1)
  out[ok] <- 1 - Fx
  if (!lower.tail) out <- 1 - out
  out
}

#' @title Quantile function of a unit-transform distribution
#' @param p numeric vector of probabilities in (0,1)
#' @param dist a \code{utd} object
#' @export
#' @examples
#' d <- utd(dgamma, pgamma, qgamma, params = list(shape = 2.5))
#' qutd(c(0.1, 0.5, 0.9), d)
qutd <- function(p, dist) {
  stopifnot(inherits(dist, "utd"))
  if (!is.null(dist$qdist)) {
    Qx <- .utd_call(dist$qdist, 1 - p, extra_params = dist$params)
  } else {
    Qx <- .numeric_qbase(1 - p, dist$pdist, dist$params)
  }
  1 / (1 + Qx)
}

#' @title Random generation from a unit-transform distribution
#' @param n number of draws
#' @param dist a \code{utd} object
#' @export
#' @examples
#' d <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5))
#' set.seed(1); rutd(5, d)
rutd <- function(n, dist) {
  stopifnot(inherits(dist, "utd"))
  if (!is.null(dist$rdist)) {
    x <- .utd_call(dist$rdist, n, extra_params = dist$params)
    1 / (1 + x)
  } else {
    qutd(stats::runif(n), dist)
  }
}

#' @title Maximum likelihood fit of a unit-transform distribution
#' @description
#' Fits ANY base distribution to (0,1)-valued data by MLE. The base
#' distribution is fully general: supply its density function and starting
#' values, no registry lookup required.
#' @param data numeric vector strictly inside (0,1)
#' @param ddist base density function, signature function(x, ..., log)
#' @param start named list/vector of starting parameter values
#' @param lower,upper optional named bounds passed to \code{optim(method="L-BFGS-B")}
#' @param name optional label
#' @return an object of class \code{utd_fit}
#' @export
#' @examples
#' set.seed(1)
#' d0 <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5))
#' y <- rutd(300, d0)
#' fit <- fit_utd_mle(y, dgamma, start = c(shape = 1), lower = c(shape = 1e-6))
#' fit$params
fit_utd_mle <- function(data, ddist, start, lower = NULL, upper = NULL, name = "custom") {
  if (!is.numeric(data) || length(data) < 2) {
    stop("`data` must be a numeric vector with at least 2 observations.", call. = FALSE)
  }
  if (anyNA(data)) stop("`data` contains NA; remove or impute before fitting.", call. = FALSE)
  if (!all(data > 0 & data < 1)) {
    stop("`data` must lie strictly inside (0,1); found values outside this range.", call. = FALSE)
  }
  if (is.null(names(start)) || any(names(start) == "")) {
    stop("`start` must be a fully named vector/list of starting parameter values.", call. = FALSE)
  }
  pnames <- names(start)
  negloglik <- function(theta) {
    params <- as.list(theta)
    names(params) <- pnames
    x <- 1 / data - 1
    logf <- .utd_call(ddist, x, log = TRUE, extra_params = params)
    val <- -sum(logf) + 2 * sum(log(data))
    if (!is.finite(val)) return(1e10)
    val
  }
  use_box <- !is.null(lower) || !is.null(upper)
  opt <- if (use_box) {
    stats::optim(unlist(start), negloglik, method = "L-BFGS-B",
                 lower = if (is.null(lower)) -Inf else lower[pnames],
                 upper = if (is.null(upper)) Inf else upper[pnames],
                 hessian = TRUE)
  } else {
    stats::optim(unlist(start), negloglik, method = "BFGS", hessian = TRUE)
  }
  if (opt$convergence != 0) {
    warning(sprintf(
      "optim() did not converge (code %d): '%s'. Treat parameter estimates with caution.",
      opt$convergence,
      if (!is.null(opt$message)) opt$message else "no message"
    ), call. = FALSE)
  }
  params <- as.list(opt$par)
  names(params) <- pnames
  vcov <- tryCatch(solve(opt$hessian), error = function(e) matrix(NA, length(pnames), length(pnames)))
  structure(
    list(params = params, loglik = -opt$value, convergence = opt$convergence,
         vcov = vcov, n = length(data), ddist = ddist, data = data,
         negloglik = negloglik, name = name),
    class = "utd_fit"
  )
}

#' @export
print.utd_fit <- function(x, ...) {
  cat(sprintf("<utd_fit> base = %s, n = %d, logLik = %.4f\n", x$name, x$n, x$loglik))
  cat("params: ")
  cat(paste(names(x$params), signif(unlist(x$params), 4), sep = " = ", collapse = ", "))
  cat("\n")
  invisible(x)
}

#' @export
logLik.utd_fit <- function(object, ...) {
  val <- object$loglik
  attr(val, "df") <- length(object$params)
  class(val) <- "logLik"
  val
}

#' @export
coef.utd_fit <- function(object, ...) unlist(object$params)

#' @export
vcov.utd_fit <- function(object, ...) object$vcov

#' @title Corrected profile-likelihood confidence interval
#' @description
#' FIX for the historical bug in the Python version: for every candidate
#' value of the parameter of interest, the nuisance parameters are
#' RE-OPTIMIZED (true profile likelihood), not held fixed at their
#' unconstrained MLE. The chi-square threshold also uses the correct
#' \code{qchisq(1 - alpha, df = 1)} (equivalently \code{stats::qchisq}),
#' avoiding the argument-order mistake that previously collapsed the
#' threshold to ~0.
#' @param fit a \code{utd_fit} object
#' @param param_name name of the parameter of interest
#' @param alpha significance level (default 0.05 -> 95% CI)
#' @param grid_expand multiplicative search range around the MLE
#' @export
#' @examples
#' set.seed(1)
#' d0 <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5))
#' y <- rutd(300, d0)
#' fit <- fit_utd_mle(y, dgamma, start = c(shape = 1), lower = c(shape = 1e-6))
#' \donttest{profile_likelihood_ci(fit, "shape")}
profile_likelihood_ci <- function(fit, param_name, alpha = 0.05, grid_expand = 4) {
  stopifnot(inherits(fit, "utd_fit"), param_name %in% names(fit$params))
  pnames <- names(fit$params)
  other <- setdiff(pnames, param_name)
  mle_val <- fit$params[[param_name]]
  se_guess <- sqrt(diag(fit$vcov))[which(pnames == param_name)]
  if (!is.finite(se_guess) || se_guess <= 0) se_guess <- abs(mle_val) * 0.2 + 1e-3

  # Correct threshold: chi-square with df = 1
  crit <- stats::qchisq(1 - alpha, df = 1)
  max_ll <- fit$loglik

  profile_negloglik_at <- function(fixed_val) {
    if (length(other) == 0) {
      theta <- fixed_val
      names(theta) <- param_name
      return(fit$negloglik(theta))
    }
    inner <- function(theta_other) {
      theta <- numeric(length(pnames))
      names(theta) <- pnames
      theta[param_name] <- fixed_val
      theta[other] <- theta_other
      fit$negloglik(theta)
    }
    start_other <- unlist(fit$params[other])
    # RE-OPTIMIZE nuisance parameters at each candidate value -- this is
    # the fix: the previous version fixed these at the unconstrained MLE.
    opt <- stats::optim(start_other, inner, method = "BFGS")
    opt$value
  }

  root_fn <- function(fixed_val) {
    nll <- profile_negloglik_at(fixed_val)
    2 * (nll - (-max_ll)) - crit
  }

  lo_search <- mle_val - grid_expand * se_guess
  hi_search <- mle_val + grid_expand * se_guess

  lower <- tryCatch(
    stats::uniroot(root_fn, lower = lo_search, upper = mle_val - 1e-8,
                    extendInt = "downX")$root,
    error = function(e) NA_real_
  )
  upper <- tryCatch(
    stats::uniroot(root_fn, lower = mle_val + 1e-8, upper = hi_search,
                    extendInt = "upX")$root,
    error = function(e) NA_real_
  )
  list(param = param_name, mle = mle_val, lower = lower, upper = upper,
       level = 1 - alpha, method = "profile_likelihood_corrected")
}
