#' Phase II monitoring of new samples
#'
#' Fits the distribution to each new sample, computes the monitored statistic
#' and compares it with the control limits of a chart built by [bootqcc()].
#'
#' @param chart An object returned by [bootqcc()].
#' @param data A matrix or data frame with one sample per row, or a list of
#'   numeric vectors.
#' @return A data frame of class `"bootqcc_monitor"` with columns `sample`,
#'   `statistic`, `LCL`, `UCL` and `status` (`"in control"`, `"out of control"`
#'   or `"fit failed"`).
#' @references
#' Shewhart, W. A. (1931). \emph{Economic Control of Quality of Manufactured
#' Product}. Van Nostrand, New York.
#'
#' Montgomery, D. C. (2019). \emph{Introduction to Statistical Quality Control},
#' 8th edition. Wiley, Hoboken, NJ.
#'
#' @examples
#' ch <- bootqcc("weibull", list(shape = 2, scale = 1), u = 0.1, n = 5,
#'               alpha = 0.05, B = 200, M = 3, arl_rep = 5, seed = 1)
#' set.seed(2)
#' new <- matrix(rweibull(50, 2, 1.6), ncol = 5)
#' bootqcc_monitor(ch, new)
#' @export
bootqcc_monitor <- function(chart, data) {
  stopifnot(inherits(chart, "bootqcc"))
  samples <- if (is.list(data) && !is.data.frame(data)) data
             else lapply(seq_len(nrow(data)), function(i) as.numeric(unlist(data[i, ])))
  lim <- chart$limits
  stat <- vapply(samples, function(x) {
    est <- chart$fitter(x)
    if (is.null(est)) return(NA_real_)
    s <- tryCatch(suppressWarnings(chart$statistic(est)), error = function(e) NA_real_)
    if (length(s) != 1L || !is.finite(s)) NA_real_ else s
  }, numeric(1))
  status <- ifelse(is.na(stat), "fit failed",
                   ifelse(stat < lim[["LCL"]] | stat > lim[["UCL"]],
                          "out of control", "in control"))
  structure(data.frame(sample = seq_along(stat), statistic = stat,
                       LCL = lim[["LCL"]], UCL = lim[["UCL"]],
                       status = status, stringsAsFactors = FALSE),
            class = c("bootqcc_monitor", "data.frame"))
}

#' Run a simulation study over a grid of settings
#'
#' Runs [bootqcc()] for every combination of sample size, quantile level and
#' parameter vector, and collects the limits and ARL summaries in one table.
#' This reproduces the Monte Carlo study usually reported in control chart
#' papers.
#'
#' @param distr,... Passed to [bootqcc()].
#' @param param_grid A data frame with one row per parameter vector and one
#'   column per parameter, named as the distribution parameters.
#' @param u_values Vector of quantile levels.
#' @param n_values Vector of subgroup sizes.
#' @param verbose Print the setting being run.
#' @param env Environment used to look up distribution functions by name.
#' @return A data frame with columns `n`, `u`, the parameters, `MLCL`, `MUCL`,
#'   `SDLCL`, `SDUCL`, `MARL`, `SDARL` and `censored`.
#' @references
#' Jones, L. A. and Woodall, W. H. (1998). The performance of bootstrap control
#' charts. \emph{Journal of Quality Technology}, 30(4), 362--375.
#'
#' Al-Lami, M. A. H., Jabbari Khamnei, H. and Heydari, A. A. (2025). A parametric
#' bootstrap control chart for Lindley Geometric percentiles.
#' \emph{PLOS ONE}, 20(2), e0316449.
#' \doi{10.1371/journal.pone.0316449}
#'
#' @examples
#' g <- expand.grid(shape = c(1.5, 2), scale = 1)
#' bootqcc_grid("weibull", g, u_values = c(0.1, 0.9), n_values = 5,
#'              alpha = 0.1, B = 100, M = 2, arl_rep = 5, seed = 1)
#' @export
bootqcc_grid <- function(distr, param_grid, u_values = 0.5, n_values = 5, ...,
                         verbose = FALSE, env = parent.frame()) {
  param_grid <- as.data.frame(param_grid)
  rows <- list()
  for (n in n_values) {
    for (u in u_values) {
      for (i in seq_len(nrow(param_grid))) {
        pars <- as.list(param_grid[i, , drop = FALSE])
        if (verbose) message(sprintf("n = %s, u = %s, %s", n, u,
                                     paste(names(pars), unlist(pars), sep = "=",
                                           collapse = ", ")))
        ch <- bootqcc(distr, pars, u = u, n = n, ..., env = env)
        rows[[length(rows) + 1L]] <- data.frame(
          n = n, u = u, param_grid[i, , drop = FALSE],
          MLCL = ch$limits[["LCL"]], MUCL = ch$limits[["UCL"]],
          SDLCL = ch$sd_limits[["LCL"]], SDUCL = ch$sd_limits[["UCL"]],
          MARL = ch$arl$mean, SDARL = ch$arl$sd, censored = ch$arl$censored,
          row.names = NULL)
      }
    }
  }
  do.call(rbind, rows)
}
