# ==== Design principles: see R/core.R ====

#' @title Pearson residuals for a fitted unit-transform model
#' @param fit a \code{utd_fit} object
#' @export
#' @examples
#' set.seed(1)
#' d0 <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5))
#' y <- rutd(100, d0)
#' fit <- fit_utd_mle(y, dgamma, start = c(shape = 1), lower = c(shape = 1e-6))
#' pearson_residuals(fit)
pearson_residuals <- function(fit) {
  stopifnot(inherits(fit, "utd_fit"))
  dist <- utd(fit$ddist, function(...) stop("pdist not needed here"), params = fit$params)
  mu <- mean(fit$data)
  (fit$data - mu) / stats::sd(fit$data)
}

#' @title Randomised quantile residuals (Dunn & Smyth, 1996)
#' @param fit a \code{utd_fit} object
#' @param pdist base CDF function corresponding to \code{fit$ddist}
#' @export
#' @examples
#' set.seed(1)
#' d0 <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5))
#' y <- rutd(100, d0)
#' fit <- fit_utd_mle(y, dgamma, start = c(shape = 1), lower = c(shape = 1e-6))
#' quantile_residuals(fit, pgamma)
quantile_residuals <- function(fit, pdist) {
  stopifnot(inherits(fit, "utd_fit"))
  dist <- utd(fit$ddist, pdist, params = fit$params)
  u <- putd(fit$data, dist)
  stats::qnorm(pmin(pmax(u, 1e-12), 1 - 1e-12))
}

#' @title Cox-Snell residuals
#' @inheritParams quantile_residuals
#' @export
#' @examples
#' set.seed(1)
#' d0 <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5))
#' y <- rutd(100, d0)
#' fit <- fit_utd_mle(y, dgamma, start = c(shape = 1), lower = c(shape = 1e-6))
#' cox_snell_residuals(fit, pgamma)
cox_snell_residuals <- function(fit, pdist) {
  dist <- utd(fit$ddist, pdist, params = fit$params)
  -log(1 - putd(fit$data, dist))
}

#' @title Approximate Cook's-distance-like case-deletion influence
#' @param fit a \code{utd_fit} object
#' @export
#' @examples
#' set.seed(1)
#' d0 <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5))
#' y <- rutd(60, d0)
#' fit <- fit_utd_mle(y, dgamma, start = c(shape = 1), lower = c(shape = 1e-6))
#' \donttest{cooks_distance(fit)}
cooks_distance <- function(fit) {
  stopifnot(inherits(fit, "utd_fit"))
  theta_hat <- unlist(fit$params)
  p <- length(theta_hat)
  n <- fit$n
  influence <- numeric(n)
  for (i in seq_len(n)) {
    fit_i <- tryCatch(
      suppressWarnings(fit_utd_mle(fit$data[-i], fit$ddist, start = as.list(theta_hat))),
      error = function(e) NULL
    )
    if (is.null(fit_i)) { influence[i] <- NA_real_; next }
    diff <- unlist(fit_i$params) - theta_hat
    influence[i] <- as.numeric(t(diff) %*% diff) / p
  }
  influence
}
