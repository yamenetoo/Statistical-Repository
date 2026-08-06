# ==== Design principles: see R/core.R ====
# Every function here accepts a `true_dist` (a utd() object) plus an
# estimator specification -- never a hard-coded distribution name. Random
# draws are generated with an explicitly supplied seed (never set globally),
# and every replicate's fit is checked for convergence before being counted.

#' @title Run a Monte Carlo study of an estimator's sampling behaviour
#' @description
#' Repeatedly simulates data from `true_dist`, refits the same base
#' distribution family, and records bias, variance, and RMSE of the
#' estimator across replications. Works for any base distribution family,
#' not only ones present in a registry.
#' @param true_dist a \code{utd} object with an \code{rdist} (or, failing
#'   that, sampling falls back to inverse-CDF via \code{qutd})
#' @param true_params named numeric vector of the true parameter values
#'   used to generate data (must match names used in \code{start})
#' @param ddist base density function used for refitting
#' @param start named starting values passed to \code{fit_utd_mle}
#' @param n sample size per replication
#' @param n_sim number of Monte Carlo replications
#' @param lower,upper optional box constraints passed to \code{fit_utd_mle}
#' @param seed optional integer seed for reproducibility
#' @return an object of class \code{utd_mc_sim}: a list with
#'   \code{estimates} (n_sim x p matrix), \code{true_params}, \code{bias},
#'   \code{variance}, \code{mse}, \code{rmse}, \code{n_converged}, \code{n_sim}
#' @export
#' @examples
#' set.seed(1)
#' true_dist <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5))
#' \donttest{
#' mc <- run_monte_carlo(true_dist, true_params = c(shape = 2.5), ddist = dgamma,
#'                        start = c(shape = 1), lower = c(shape = 1e-6),
#'                        n = 100, n_sim = 200, seed = 42)
#' summary(mc)
#' }
run_monte_carlo <- function(true_dist, true_params, ddist, start, n, n_sim = 1000,
                             lower = NULL, upper = NULL, seed = NULL) {
  if (!inherits(true_dist, "utd")) stop("`true_dist` must be a utd object.", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  pnames <- names(true_params)
  if (is.null(pnames) || any(pnames == "")) {
    stop("`true_params` must be a fully named numeric vector.", call. = FALSE)
  }
  estimates <- matrix(NA_real_, nrow = n_sim, ncol = length(pnames),
                       dimnames = list(NULL, pnames))
  n_converged <- 0L
  for (i in seq_len(n_sim)) {
    y <- rutd(n, true_dist)
    fit <- tryCatch(
      suppressWarnings(fit_utd_mle(y, ddist, start = start, lower = lower, upper = upper)),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      estimates[i, ] <- unlist(fit$params)[pnames]
      n_converged <- n_converged + (fit$convergence == 0)
    }
  }
  bias <- colMeans(estimates, na.rm = TRUE) - true_params
  variance <- apply(estimates, 2, stats::var, na.rm = TRUE)
  mse <- bias^2 + variance
  structure(
    list(estimates = estimates, true_params = true_params, bias = bias,
         variance = variance, mse = mse, rmse = sqrt(mse),
         n_converged = n_converged, n_sim = n_sim, n = n),
    class = "utd_mc_sim"
  )
}

#' @export
summary.utd_mc_sim <- function(object, ...) {
  df <- data.frame(
    true = object$true_params,
    mean_estimate = colMeans(object$estimates, na.rm = TRUE),
    bias = object$bias,
    sd = sqrt(object$variance),
    rmse = object$rmse
  )
  structure(list(table = df, n = object$n, n_sim = object$n_sim,
                 n_converged = object$n_converged),
            class = "summary.utd_mc_sim")
}

#' @export
print.summary.utd_mc_sim <- function(x, ...) {
  cat(sprintf("Monte Carlo study: n = %d, replications = %d (%d converged)\n",
              x$n, x$n_sim, x$n_converged))
  print(round(x$table, 4))
  invisible(x)
}

#' @title Compare several estimators on the same simulated data
#' @param true_dist a \code{utd} object to simulate from
#' @param true_params named numeric vector of true parameter values
#' @param estimator_specs a named list of estimator specifications, each a
#'   list with elements \code{ddist}, \code{start}, and optional
#'   \code{lower}/\code{upper} -- e.g.
#'   \code{list(mle = list(ddist = dgamma, start = c(shape = 1)))}
#' @param n,n_sim,seed as in \code{run_monte_carlo}
#' @export
#' @examples
#' true_dist <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5))
#' specs <- list(mle = list(ddist = dgamma, start = c(shape = 1), lower = c(shape = 1e-6)))
#' \donttest{compare_estimators_mc(true_dist, c(shape = 2.5), specs, n = 100, n_sim = 100, seed = 1)}
compare_estimators_mc <- function(true_dist, true_params, estimator_specs, n, n_sim = 1000, seed = NULL) {
  if (!is.list(estimator_specs) || is.null(names(estimator_specs))) {
    stop("`estimator_specs` must be a named list of specifications.", call. = FALSE)
  }
  results <- lapply(names(estimator_specs), function(nm) {
    spec <- estimator_specs[[nm]]
    mc <- run_monte_carlo(true_dist, true_params, spec$ddist, spec$start, n, n_sim,
                           lower = spec$lower, upper = spec$upper, seed = seed)
    data.frame(estimator = nm, param = names(true_params),
               bias = mc$bias, rmse = mc$rmse, n_converged = mc$n_converged)
  })
  do.call(rbind, results)
}

#' @title Study estimator bias/MSE as a function of sample size
#' @param true_dist a \code{utd} object
#' @param true_params named numeric vector of true parameter values
#' @param ddist,start,lower,upper as in \code{fit_utd_mle}
#' @param n_grid numeric vector of sample sizes to evaluate
#' @param n_sim replications per sample size
#' @param seed optional seed
#' @export
#' @examples
#' true_dist <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5))
#' \donttest{
#' mc_sample_size_study(true_dist, c(shape = 2.5), dgamma, c(shape = 1),
#'                       lower = c(shape = 1e-6), n_grid = c(30, 100), n_sim = 50, seed = 1)
#' }
mc_sample_size_study <- function(true_dist, true_params, ddist, start, lower = NULL,
                                  upper = NULL, n_grid, n_sim = 500, seed = NULL) {
  out <- lapply(n_grid, function(n) {
    mc <- run_monte_carlo(true_dist, true_params, ddist, start, n, n_sim, lower, upper, seed)
    data.frame(n = n, param = names(true_params), bias = mc$bias,
               variance = mc$variance, mse = mc$mse, rmse = mc$rmse)
  })
  do.call(rbind, out)
}

#' @export
plot.utd_mc_sim <- function(x, param = NULL, type = c("density", "bias"), ...) {
  type <- match.arg(type)
  if (is.null(param)) param <- colnames(x$estimates)[1]
  vals <- x$estimates[, param]
  vals <- vals[!is.na(vals)]
  if (type == "density") {
    plot(stats::density(vals), main = sprintf("Sampling distribution: %s", param),
         xlab = param)
    graphics::abline(v = x$true_params[[param]], col = "red", lty = 2)
  } else {
    graphics::barplot(x$bias, main = "Bias per parameter", ylab = "bias")
    graphics::abline(h = 0)
  }
  invisible(x)
}
