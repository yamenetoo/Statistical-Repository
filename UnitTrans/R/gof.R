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
#' valid with estimated parameters). The p-value is
#' `(1 + number of T* >= T) / (B + 1)`, so it is never exactly zero.
#'
#' @param fit A `transfit` object.
#' @param B Number of bootstrap replications.
#' @param refit Re-estimate the parameters in each replication.
#' @param seed Optional seed (the caller's random number stream is restored).
#' @return A data frame with statistic and bootstrap p-value for each test, and
#'   `valid_boot`, the number of replications that could be fitted. A warning
#'   is issued when some replications fail.
#' @examples
#' set.seed(3)
#' fit <- fit_trans(rtrans(100, "norm", c(0, 1), "logit"), "norm", "logit")
#' gof_trans(fit, B = 50)
#' @export
gof_trans <- function(fit, B = 200, refit = TRUE, seed = NULL) {
  stopifnot(inherits(fit, "transfit"))
  bl <- .as_base(fit$dist); tr <- get_transformation(fit$trans)
  obs <- .gof_stats(ptrans(fit$data, bl, fit$coef, tr))
  boot <- .with_seed(seed, {
    m <- matrix(NA_real_, B, 3, dimnames = list(NULL, names(obs)))
    for (b in seq_len(B)) {
      xs <- .clip01(rtrans(fit$n, bl, fit$coef, tr))
      par_b <- if (refit)
        tryCatch(.fit_core(xs, bl, tr, start = fit$coef, lower = fit$lower,
                           upper = fit$upper, n_start = 1, hessian = FALSE,
                           check_support = FALSE)$coef,
                 error = function(e) NULL) else fit$coef
      if (is.null(par_b)) next
      m[b, ] <- .gof_stats(ptrans(xs, bl, par_b, tr))
    }
    m
  })
  ok <- stats::complete.cases(boot)
  nv <- sum(ok)
  if (nv == 0L) stop("All bootstrap replications failed.", call. = FALSE)
  if (nv < B)
    warning(B - nv, " of ", B, " bootstrap replications failed and were dropped.",
            call. = FALSE)
  pv <- vapply(names(obs), function(s) (1 + sum(boot[ok, s] >= obs[[s]])) / (nv + 1), 0)
  data.frame(test = names(obs), statistic = unname(obs), p.value = unname(pv),
             valid_boot = nv, row.names = NULL)
}
