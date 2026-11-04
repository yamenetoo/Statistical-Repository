#' @export
print.bootqcc <- function(x, digits = 4, ...) {
  what <- if (is.na(x$u)) "custom statistic" else sprintf("%g-quantile", x$u)
  cat("Parametric bootstrap control chart\n")
  cat(sprintf("Distribution : %s (%s)\n", x$distr_name,
              paste(names(x$params), unlist(x$params), sep = " = ", collapse = ", ")))
  cat(sprintf("Statistic    : %s, subgroup size n = %d\n", what, x$n))
  cat(sprintf("Settings     : alpha = %g, B = %d, M = %d\n", x$alpha, x$B, x$M))
  cat(sprintf("LCL = %s (SD %s)\n", format(x$limits[["LCL"]], digits = digits),
              format(x$sd_limits[["LCL"]], digits = digits)))
  cat(sprintf("UCL = %s (SD %s)\n", format(x$limits[["UCL"]], digits = digits),
              format(x$sd_limits[["UCL"]], digits = digits)))
  if (!is.null(x$arl)) {
    cat(sprintf("In-control ARL = %s (SD %s) from %d runs",
                format(x$arl$mean, digits = digits),
                format(x$arl$sd, digits = digits), length(x$arl$runs)))
    if (x$arl$censored > 0) cat(sprintf(", %d censored", x$arl$censored))
    cat("\n")
  }
  invisible(x)
}

#' Plot a bootstrap control chart
#'
#' Without `data`, shows the distribution of the bootstrap statistic with the
#' control limits. With `data` (a matrix, data frame or list of samples), draws
#' the Phase II control chart of those samples.
#'
#' @param x An object returned by [bootqcc()].
#' @param data Optional new samples, see [bootqcc_monitor()].
#' @param ... Further graphical parameters.
#' @return `x`, invisibly.
#' @export
plot.bootqcc <- function(x, data = NULL, ...) {
  lim <- x$limits
  if (is.null(data)) {
    graphics::hist(x$boot_stats, breaks = "FD", col = "grey85", border = "white",
                   main = "Bootstrap distribution of the statistic",
                   xlab = "Statistic", ...)
    graphics::abline(v = lim, col = "red", lty = 2, lwd = 2)
  } else {
    plot(bootqcc_monitor(x, data), ...)
  }
  invisible(x)
}

#' @export
plot.bootqcc_monitor <- function(x, ...) {
  rng <- range(c(x$statistic, x$LCL[1], x$UCL[1]), na.rm = TRUE)
  out <- !is.na(x$statistic) & x$status == "out of control"
  graphics::plot(x$sample, x$statistic, type = "b", pch = 19,
                 col = ifelse(out, "red", "black"), ylim = rng,
                 xlab = "Sample", ylab = "Statistic",
                 main = "Bootstrap control chart", ...)
  graphics::abline(h = c(x$LCL[1], x$UCL[1]), col = "red", lty = 2, lwd = 2)
  invisible(x)
}
