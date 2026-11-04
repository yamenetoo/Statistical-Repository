#' Monte-Carlo simulation study of the MLE
#'
#' Repeatedly simulates from a transformed distribution, refits it and
#' summarizes bias, root mean squared error and the coverage of Wald intervals.
#' Coverage is computed over the replications with available standard errors.
#'
#' @param dist,par,trans True model.
#' @param n Sample size(s).
#' @param R Number of replications.
#' @param seed Optional seed (the caller's random number stream is restored).
#' @return A data frame with one row per sample size and parameter, including
#'   the number of `failed` replications.
#' @examples
#' simulation_study("norm", c(mean = 0, sd = 1), "logit", n = c(50, 100), R = 20, seed = 1)
#' @export
simulation_study <- function(dist, par, trans, n = c(50, 100, 200), R = 200, seed = NULL) {
  bl <- .as_base(dist); tr <- get_transformation(trans)
  par <- .name_par(par, bl)
  zcrit <- stats::qnorm(0.975)
  .with_seed(seed, {
    out <- list()
    for (nn in n) {
      est <- matrix(NA_real_, R, length(par)); cov <- est
      for (i in seq_len(R)) {
        x <- .clip01(rtrans(nn, bl, par, tr))
        f <- tryCatch(.fit_core(x, bl, tr, start = par, n_start = 1,
                                hessian = TRUE, check_support = FALSE),
                      error = function(e) NULL)
        if (is.null(f)) next
        est[i, ] <- f$coef
        cov[i, ] <- abs(f$coef - par) <= zcrit * f$se
      }
      out[[length(out) + 1]] <- data.frame(
        n = nn, parameter = names(par), true = unname(par),
        mean_est = colMeans(est, na.rm = TRUE),
        bias = colMeans(est, na.rm = TRUE) - par,
        rmse = sqrt(colMeans((sweep(est, 2, par))^2, na.rm = TRUE)),
        coverage = colMeans(cov, na.rm = TRUE),
        failed = sum(!stats::complete.cases(est)), row.names = NULL)
    }
    do.call(rbind, out)
  })
}
