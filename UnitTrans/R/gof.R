.gof_stats <- function(u) {
  u <- sort(pmin(pmax(u, 1e-12), 1 - 1e-12)); n <- length(u); i <- seq_len(n)
  ks <- max(pmax(i / n - u, u - (i - 1) / n))
  cvm <- 1 / (12 * n) + sum((u - (2 * i - 1) / (2 * n))^2)
  ad <- -n - mean((2 * i - 1) * (log(u) + log(1 - rev(u))))
  c(KS = ks, CvM = cvm, AD = ad)
}

#' Goodness-of-fit tests with parametric bootstrap
#'
#' Computes Kolmogorov-Smirnov, Cramer-von Mises and Anderson-Darling
#' statistics for a fitted model and obtains p-values by a parametric bootstrap
#' that re-estimates the parameters in every replication (so the p-values are
#' valid with estimated parameters).
#'
#' @param fit A `transfit` object.
#' @param B Number of bootstrap replications.
#' @param refit Re-estimate the parameters in each replication.
#' @param seed Optional seed.
#' @return A data frame with statistic and bootstrap p-value for each test.
#' @examples
#' set.seed(3)
#' fit <- fit_trans(rtrans(100, "norm", c(0, 1), "logit"), "norm", "logit")
#' gof_trans(fit, B = 50)
#' @export
gof_trans <- function(fit, B = 200, refit = TRUE, seed = NULL) {
  stopifnot(inherits(fit, "transfit"))
  if (!is.null(seed)) set.seed(seed)
  obs <- .gof_stats(ptrans(fit$data, fit$dist, fit$coef, fit$trans))
  boot <- matrix(NA_real_, B, 3, dimnames = list(NULL, names(obs)))
  for (b in seq_len(B)) {
    xs <- rtrans(fit$n, fit$dist, fit$coef, fit$trans)
    xs <- pmin(pmax(xs, 1e-10), 1 - 1e-10)
    par_b <- if (refit)
      tryCatch(fit_trans(xs, fit$dist, fit$trans, start = fit$coef, n_start = 1)$coef,
               error = function(e) NULL) else fit$coef
    if (is.null(par_b)) next
    boot[b, ] <- .gof_stats(ptrans(xs, fit$dist, par_b, fit$trans))
  }
  pv <- vapply(names(obs), function(s) mean(boot[, s] >= obs[s], na.rm = TRUE), 0)
  data.frame(test = names(obs), statistic = unname(obs), p.value = unname(pv),
             valid_boot = sum(stats::complete.cases(boot)), row.names = NULL)
}
