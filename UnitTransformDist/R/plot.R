# ==== Design principles: see R/core.R ====
# Base graphics only -- ggplot2 stays in Suggests, never Imports, so the
# package's plotting functions work even if ggplot2 is not installed.

#' @title Plot the density of a unit-transform distribution
#' @param dist a \code{utd} object
#' @param n_points number of grid points
#' @export
#' @examples
#' d <- utd(dgamma, pgamma, qgamma, params = list(shape = 2.5))
#' plot_pdf(d)
plot_pdf <- function(dist, n_points = 500) {
  y <- seq(1e-4, 1 - 1e-4, length.out = n_points)
  plot(y, dutd(y, dist), type = "l", xlab = "y", ylab = "f(y)",
       main = sprintf("Density: %s", dist$name))
}

#' @title Plot the CDF of a unit-transform distribution
#' @inheritParams plot_pdf
#' @export
#' @examples
#' d <- utd(dgamma, pgamma, qgamma, params = list(shape = 2.5))
#' plot_cdf(d)
plot_cdf <- function(dist, n_points = 500) {
  y <- seq(1e-4, 1 - 1e-4, length.out = n_points)
  plot(y, putd(y, dist), type = "l", xlab = "y", ylab = "F(y)",
       main = sprintf("CDF: %s", dist$name))
}

#' @title Quantile-quantile plot against a fitted unit-transform model
#' @param fit a \code{utd_fit} object
#' @param pdist base CDF function
#' @export
#' @examples
#' set.seed(1)
#' d0 <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5))
#' y <- rutd(100, d0)
#' fit <- fit_utd_mle(y, dgamma, start = c(shape = 1), lower = c(shape = 1e-6))
#' qq_plot(fit, pgamma)
qq_plot <- function(fit, pdist) {
  dist <- utd(fit$ddist, pdist, params = fit$params)
  n <- fit$n
  emp <- sort(fit$data)
  theo <- qutd((seq_len(n) - 0.5) / n, dist)
  plot(theo, emp, xlab = "theoretical quantiles", ylab = "empirical quantiles",
       main = "Q-Q plot")
  graphics::abline(0, 1, col = "red", lty = 2)
}
