# ==== Design principles: see R/core.R ====
# Each method below is a distinct, well-established bias-reduction strategy
# from the statistics literature (Cox & Snell 1968 for the analytic term;
# Quenouille/Tukey for jackknife; Efron for bootstrap; Firth 1993 for
# penalized likelihood). No method is presented as universally best --
# mc_evaluate_bias_correction() lets the user judge empirically for their
# own base distribution and sample size.

#' @title Jackknife bias-corrected estimate
#' @description Leave-one-out jackknife bias correction, applicable to any
#'   base distribution fitted via \code{fit_utd_mle}.
#' @param data numeric vector in (0,1)
#' @param ddist base density function
#' @param start,lower,upper as in \code{fit_utd_mle}
#' @return list with \code{original}, \code{jackknife_estimate},
#'   \code{bias_estimate}, \code{corrected}
#' @export
#' @examples
#' set.seed(1)
#' d0 <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5))
#' y <- rutd(80, d0)
#' \donttest{bias_correct_jackknife(y, dgamma, start = c(shape = 1), lower = c(shape = 1e-6))}
bias_correct_jackknife <- function(data, ddist, start, lower = NULL, upper = NULL) {
  n <- length(data)
  full <- fit_utd_mle(data, ddist, start, lower, upper)
  theta_hat <- unlist(full$params)
  loo <- vapply(seq_len(n), function(i) {
    fit_i <- tryCatch(
      suppressWarnings(fit_utd_mle(data[-i], ddist, start, lower, upper)),
      error = function(e) NULL
    )
    if (is.null(fit_i)) return(rep(NA_real_, length(theta_hat)))
    unlist(fit_i$params)
  }, numeric(length(theta_hat)))
  loo <- matrix(loo, ncol = n)
  jack_mean <- rowMeans(loo, na.rm = TRUE)
  bias_est <- (n - 1) * (jack_mean - theta_hat)
  corrected <- theta_hat - bias_est
  names(corrected) <- names(theta_hat)
  list(original = theta_hat, jackknife_estimate = jack_mean,
       bias_estimate = bias_est, corrected = corrected)
}

#' @title Bootstrap bias-corrected estimate
#' @param fit a \code{utd_fit} object (used for parametric resampling)
#' @param data original data (used for nonparametric resampling)
#' @param n_boot number of bootstrap replications
#' @param type "parametric" resamples from the fitted \code{utd}
#'   distribution; "nonparametric" resamples the data with replacement
#' @param sim_dist for type = "parametric": a \code{utd} object built from
#'   \code{fit$params} used to generate bootstrap samples
#' @export
#' @examples
#' set.seed(1)
#' d0 <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5))
#' y <- rutd(150, d0)
#' fit <- fit_utd_mle(y, dgamma, start = c(shape = 1), lower = c(shape = 1e-6))
#' \donttest{
#' bias_correct_bootstrap(fit, y, n_boot = 200, type = "nonparametric")
#' }
bias_correct_bootstrap <- function(fit, data, n_boot = 1000,
                                    type = c("nonparametric", "parametric"),
                                    sim_dist = NULL) {
  type <- match.arg(type)
  theta_hat <- unlist(fit$params)
  n <- length(data)
  boots <- matrix(NA_real_, nrow = n_boot, ncol = length(theta_hat))
  for (b in seq_len(n_boot)) {
    y_b <- if (type == "nonparametric") {
      sample(data, n, replace = TRUE)
    } else {
      if (is.null(sim_dist)) stop("`sim_dist` is required for type = 'parametric'.", call. = FALSE)
      rutd(n, sim_dist)
    }
    fit_b <- tryCatch(
      suppressWarnings(fit_utd_mle(y_b, fit$ddist, start = as.list(theta_hat))),
      error = function(e) NULL
    )
    if (!is.null(fit_b)) boots[b, ] <- unlist(fit_b$params)
  }
  boot_mean <- colMeans(boots, na.rm = TRUE)
  bias_est <- boot_mean - theta_hat
  corrected <- theta_hat - bias_est
  names(corrected) <- names(theta_hat)
  list(original = theta_hat, bootstrap_mean = boot_mean,
       bias_estimate = bias_est, corrected = corrected, boots = boots)
}

#' @title Median-unbiased estimate via bootstrap median
#' @description A robust alternative to mean-based bias correction: uses
#'   the bootstrap median rather than the bootstrap mean, which is less
#'   sensitive to a skewed sampling distribution of the estimator.
#' @inheritParams bias_correct_bootstrap
#' @export
#' @examples
#' set.seed(1)
#' d0 <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5))
#' y <- rutd(150, d0)
#' fit <- fit_utd_mle(y, dgamma, start = c(shape = 1), lower = c(shape = 1e-6))
#' \donttest{bias_correct_median(fit, y, n_boot = 200)}
bias_correct_median <- function(fit, data, n_boot = 1000) {
  theta_hat <- unlist(fit$params)
  n <- length(data)
  boots <- matrix(NA_real_, nrow = n_boot, ncol = length(theta_hat))
  for (b in seq_len(n_boot)) {
    y_b <- sample(data, n, replace = TRUE)
    fit_b <- tryCatch(
      suppressWarnings(fit_utd_mle(y_b, fit$ddist, start = as.list(theta_hat))),
      error = function(e) NULL
    )
    if (!is.null(fit_b)) boots[b, ] <- unlist(fit_b$params)
  }
  boot_median <- apply(boots, 2, stats::median, na.rm = TRUE)
  corrected <- 2 * theta_hat - boot_median
  names(corrected) <- names(theta_hat)
  list(original = theta_hat, bootstrap_median = boot_median, corrected = corrected)
}

#' @title Firth-type penalized-likelihood bias reduction
#' @description
#' Adds Jeffreys'-prior-style penalty \eqn{+ \tfrac12 \log|I(\theta)|} to the
#' log-likelihood before maximizing, following Firth (1993). The
#' information matrix is approximated numerically from the Hessian of the
#' unpenalized log-likelihood at each trial point, so this works for any
#' base distribution without requiring an analytic expected-information
#' formula.
#' @param data numeric vector in (0,1)
#' @param ddist base density function
#' @param start,lower,upper as in \code{fit_utd_mle}
#' @export
#' @examples
#' set.seed(1)
#' d0 <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5))
#' y <- rutd(60, d0)
#' \donttest{bias_correct_firth(y, dgamma, start = c(shape = 1), lower = c(shape = 1e-6))}
bias_correct_firth <- function(data, ddist, start, lower = NULL, upper = NULL) {
  pnames <- names(start)
  base_negloglik <- function(theta) {
    params <- as.list(theta); names(params) <- pnames
    x <- 1 / data - 1
    val <- -sum(.utd_call(ddist, x, log = TRUE, extra_params = params)) + 2 * sum(log(data))
    if (!is.finite(val)) return(1e10)
    val
  }
  penalized_negloglik <- function(theta) {
    h <- tryCatch(numDeriv_hessian(base_negloglik, theta), error = function(e) NULL)
    penalty <- if (is.null(h)) 0 else {
      logdet <- tryCatch(as.numeric(determinant(h, logarithm = TRUE)$modulus),
                          error = function(e) 0)
      -0.5 * logdet
    }
    base_negloglik(theta) + penalty
  }
  use_box <- !is.null(lower) || !is.null(upper)
  opt <- if (use_box) {
    stats::optim(unlist(start), penalized_negloglik, method = "L-BFGS-B",
                 lower = if (is.null(lower)) -Inf else lower[pnames],
                 upper = if (is.null(upper)) Inf else upper[pnames])
  } else {
    stats::optim(unlist(start), penalized_negloglik, method = "BFGS")
  }
  params <- as.list(opt$par); names(params) <- pnames
  list(corrected = unlist(params), convergence = opt$convergence)
}

#' Finite-difference Hessian (dependency-free stand-in for numDeriv::hessian)
.numDeriv_hessian_step <- 1e-4
numDeriv_hessian <- function(f, x, h = .numDeriv_hessian_step) {
  n <- length(x)
  H <- matrix(0, n, n)
  f0 <- f(x)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      xi <- x; xj <- x
      xi[i] <- xi[i] + h; xj[j] <- xj[j] + h
      xij <- x; xij[i] <- xij[i] + h; xij[j] <- xij[j] + h
      H[i, j] <- (f(xij) - f(xi) - f(xj) + f0) / (h * h)
    }
  }
  (H + t(H)) / 2
}

#' @title Double bootstrap bias correction
#' @description
#' Corrects for second-order bias by bootstrapping the bootstrap bias
#' estimate itself (a "bootstrap-after-bootstrap"). Computationally
#' expensive (\code{n_outer * n_inner} fits); intended for small samples
#' where first-order bootstrap correction is suspected to be insufficient.
#' @param data numeric vector in (0,1)
#' @param ddist base density function
#' @param start,lower,upper as in \code{fit_utd_mle}
#' @param n_outer,n_inner outer/inner bootstrap replication counts
#' @export
#' @examples
#' set.seed(1)
#' d0 <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5))
#' y <- rutd(60, d0)
#' \donttest{
#' bias_correct_double_bootstrap(y, dgamma, start = c(shape = 1),
#'                                lower = c(shape = 1e-6), n_outer = 20, n_inner = 20)
#' }
bias_correct_double_bootstrap <- function(data, ddist, start, lower = NULL, upper = NULL,
                                           n_outer = 200, n_inner = 100) {
  n <- length(data)
  fit0 <- fit_utd_mle(data, ddist, start, lower, upper)
  theta_hat <- unlist(fit0$params)
  outer_bias <- matrix(NA_real_, nrow = n_outer, ncol = length(theta_hat))
  for (o in seq_len(n_outer)) {
    y_o <- sample(data, n, replace = TRUE)
    fit_o <- tryCatch(suppressWarnings(fit_utd_mle(y_o, ddist, as.list(theta_hat), lower, upper)),
                       error = function(e) NULL)
    if (is.null(fit_o)) next
    theta_o <- unlist(fit_o$params)
    inner_ests <- matrix(NA_real_, nrow = n_inner, ncol = length(theta_hat))
    for (k in seq_len(n_inner)) {
      y_ok <- sample(y_o, n, replace = TRUE)
      fit_ok <- tryCatch(suppressWarnings(fit_utd_mle(y_ok, ddist, as.list(theta_o), lower, upper)),
                          error = function(e) NULL)
      if (!is.null(fit_ok)) inner_ests[k, ] <- unlist(fit_ok$params)
    }
    outer_bias[o, ] <- colMeans(inner_ests, na.rm = TRUE) - theta_o
  }
  second_order_bias <- colMeans(outer_bias, na.rm = TRUE)
  corrected <- theta_hat - second_order_bias
  names(corrected) <- names(theta_hat)
  list(original = theta_hat, second_order_bias = second_order_bias, corrected = corrected)
}

#' @title Unified bias-correction dispatcher
#' @param fit a \code{utd_fit} object
#' @param data original data
#' @param method one of "jackknife", "bootstrap", "median", "firth", "double_bootstrap"
#' @param ... passed through to the underlying method
#' @export
#' @examples
#' set.seed(1)
#' d0 <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5))
#' y <- rutd(100, d0)
#' fit <- fit_utd_mle(y, dgamma, start = c(shape = 1), lower = c(shape = 1e-6))
#' \donttest{bias_correct(fit, y, method = "jackknife")}
bias_correct <- function(fit, data, method = c("jackknife", "bootstrap", "median",
                                                 "firth", "double_bootstrap"), ...) {
  method <- match.arg(method)
  switch(method,
    jackknife        = bias_correct_jackknife(data, fit$ddist, start = fit$params, ...),
    bootstrap        = bias_correct_bootstrap(fit, data, ...),
    median           = bias_correct_median(fit, data, ...),
    firth            = bias_correct_firth(data, fit$ddist, start = fit$params, ...),
    double_bootstrap = bias_correct_double_bootstrap(data, fit$ddist, start = fit$params, ...)
  )
}
