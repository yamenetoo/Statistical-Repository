# ==== Design principles: see R/core.R ====

#' @export
AIC.utd_fit <- function(object, ..., k = 2) {
  -2 * object$loglik + k * length(object$params)
}

#' @export
BIC.utd_fit <- function(object, ...) {
  -2 * object$loglik + log(object$n) * length(object$params)
}

#' @title Wald confidence interval
#' @param fit a \code{utd_fit} object
#' @param param_name parameter name
#' @param alpha significance level
#' @export
#' @examples
#' set.seed(1)
#' d0 <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5))
#' y <- rutd(200, d0)
#' fit <- fit_utd_mle(y, dgamma, start = c(shape = 1), lower = c(shape = 1e-6))
#' wald_ci(fit, "shape")
wald_ci <- function(fit, param_name, alpha = 0.05) {
  se <- sqrt(diag(fit$vcov))[which(names(fit$params) == param_name)]
  est <- fit$params[[param_name]]
  z <- stats::qnorm(1 - alpha / 2)
  c(lower = est - z * se, upper = est + z * se)
}

#' @title Nonparametric bootstrap confidence interval
#' @param fit a \code{utd_fit} object
#' @param data original data
#' @param n_boot number of bootstrap replications
#' @param alpha significance level
#' @export
#' @examples
#' set.seed(1)
#' d0 <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5))
#' y <- rutd(120, d0)
#' fit <- fit_utd_mle(y, dgamma, start = c(shape = 1), lower = c(shape = 1e-6))
#' \donttest{bootstrap_ci(fit, y, n_boot = 200)}
bootstrap_ci <- function(fit, data, n_boot = 1000, alpha = 0.05) {
  theta_hat <- unlist(fit$params)
  n <- length(data)
  boots <- matrix(NA_real_, nrow = n_boot, ncol = length(theta_hat))
  for (b in seq_len(n_boot)) {
    y_b <- sample(data, n, replace = TRUE)
    fit_b <- tryCatch(suppressWarnings(fit_utd_mle(y_b, fit$ddist, as.list(theta_hat))),
                       error = function(e) NULL)
    if (!is.null(fit_b)) boots[b, ] <- unlist(fit_b$params)
  }
  apply(boots, 2, stats::quantile, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
}

#' @title Kolmogorov-Smirnov goodness-of-fit test
#' @param fit a \code{utd_fit} object
#' @param pdist base CDF function
#' @export
#' @examples
#' set.seed(1)
#' d0 <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5))
#' y <- rutd(100, d0)
#' fit <- fit_utd_mle(y, dgamma, start = c(shape = 1), lower = c(shape = 1e-6))
#' ks_test(fit, pgamma)
ks_test <- function(fit, pdist) {
  dist <- utd(fit$ddist, pdist, params = fit$params)
  stats::ks.test(fit$data, function(y) putd(y, dist))
}

#' @title Hazard function of a fitted unit-transform distribution
#' @param fit a \code{utd_fit} object
#' @param pdist base CDF function
#' @param y evaluation grid in (0,1)
#' @export
#' @examples
#' set.seed(1)
#' d0 <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5))
#' y <- rutd(100, d0)
#' fit <- fit_utd_mle(y, dgamma, start = c(shape = 1), lower = c(shape = 1e-6))
#' hazard(fit, pgamma, seq(0.1, 0.9, 0.1))
hazard <- function(fit, pdist, y) {
  dist <- utd(fit$ddist, pdist, params = fit$params)
  dutd(y, dist) / (1 - putd(y, dist))
}
