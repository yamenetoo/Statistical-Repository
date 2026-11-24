## Plotting -------------------------------------------------------------------

.polyline <- function(A) {
  ## knots plus the drops to zero outside the set
  x <- A$x
  m <- A$mu
  n <- length(x)
  xs <- c(x[1L], x, x[n])
  ms <- c(0, m, 0)
  list(x = xs, y = ms)
}

plot.fuzzyset <- function(x, xlim = NULL, ylim = c(0, 1.05),
                          xlab = "x", ylab = expression(mu(x)),
                          main = if (nzchar(x$name)) x$name else "",
                          col = "black", lwd = 2, ...) {
  p <- .polyline(x)
  if (is.null(xlim)) {
    r <- range(x$x)
    pad <- if (diff(r) > 0) 0.05 * diff(r) else 1
    xlim <- r + c(-pad, pad)
  }
  graphics::plot(p$x, p$y, type = "l", xlim = xlim, ylim = ylim, xlab = xlab,
                 ylab = ylab, main = main, col = col, lwd = lwd, ...)
  invisible(x)
}

lines.fuzzyset <- function(x, col = "black", lwd = 2, ...) {
  p <- .polyline(x)
  graphics::lines(p$x, p$y, col = col, lwd = lwd, ...)
  invisible(x)
}

plot_fuzzy_sets <- function(sets, col = NULL, lwd = 2, lty = 1, xlim = NULL,
                            ylim = c(0, 1.05), xlab = "x",
                            ylab = expression(mu(x)), main = "",
                            legend = TRUE, legend_pos = "topright", ...) {
  sets <- .as_set_list(sets)
  k <- length(sets)
  if (is.null(col)) col <- grDevices::hcl.colors(max(k, 2L), "Dark 3")[seq_len(k)]
  col <- rep_len(col, k)
  lty <- rep_len(lty, k)
  if (is.null(xlim)) {
    r <- range(unlist(lapply(sets, `[[`, "x")))
    pad <- if (diff(r) > 0) 0.05 * diff(r) else 1
    xlim <- r + c(-pad, pad)
  }
  graphics::plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
                 main = main, ...)
  for (i in seq_len(k)) lines.fuzzyset(sets[[i]], col = col[i], lwd = lwd, lty = lty[i])
  nm <- vapply(sets, `[[`, character(1), "name")
  if (isTRUE(legend) && any(nzchar(nm)))
    graphics::legend(legend_pos, legend = nm, col = col, lty = lty, lwd = lwd,
                     bty = "n")
  invisible(sets)
}

## Example sets from Ross (2000) ---------------------------------------------

ross2000_sets <- function(which = c("map", "output", "height", "voltage")) {
  which <- match.arg(which)
  switch(which,
         ## Figures 1.20-1.22: three measurements of mean arterial pressure (mm Hg)
         map = list(B1 = fs_trapezoidal(0, 1, 4, 5, height = 0.3, name = "B1"),
                    B2 = fs_trapezoidal(3, 4, 6, 7, height = 0.5, name = "B2"),
                    B3 = fs_trapezoidal(5, 6, 7, 8, height = 1, name = "B3")),
         ## Figure 1.15: the two parts of a typical fuzzy output
         output = list(C1 = fs_trapezoidal(0, 2, 6, 8, height = 0.5, name = "C1"),
                       C2 = fs_triangular(6, 8, 10, name = "C2")),
         ## Figure 1.1: crisp set A = [5, 7] and fuzzy set H = "near 6 feet"
         height = list(A = fs_rectangular(5, 7, name = "A (crisp)"),
                       H = fs_triangular(4.5, 6, 7.5, name = "H (near 6 ft)")),
         ## Figure 1.7: "low voltage" and "medium voltage" sets (illustrative)
         voltage = list(low = fs_triangular(0, 3, 7, name = "low voltage"),
                        medium = fs_trapezoidal(5, 7, 9, 12, name = "medium voltage")))
}
