## Methods for "bkwgfit" objects -------------------------------------------------

.model_label <- function(fit) {
  sm <- bkwg_submodels()
  paste0(sm$name[sm$submodel == fit$submodel], " (baseline: ",
         fit$baseline$name, ")")
}

print.bkwgfit <- function(x, digits = 4, ...) {
  cat("\n", .model_label(x), "\n", sep = "")
  cat("Method: ", .method_label(x$method), ";  n = ", x$n, "\n\n", sep = "")
  if (x$k) {
    cat("Estimates:\n")
    print(signif(x$coefficients, digits))
  } else {
    cat("No free parameters.\n")
  }
  if (length(x$fixed)) {
    cat("Fixed: ", paste(names(x$fixed), "=", format(x$fixed), collapse = ", "),
        "\n", sep = "")
  }
  cat("\nlogLik = ", format(x$loglik, digits = digits + 2),
      ",  AIC = ", format(x$AIC, digits = digits + 2),
      ",  BIC = ", format(x$BIC, digits = digits + 2), "\n", sep = "")
  invisible(x)
}

summary.bkwgfit <- function(object, ...) {
  est <- object$coefficients
  se <- object$se
  z <- est / se
  tab <- cbind(Estimate = est, `Std. Error` = se, `z value` = z,
               `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  ic <- c(logLik = object$loglik, AIC = object$AIC, CAIC = object$CAIC,
          BIC = object$BIC, HQIC = object$HQIC)
  structure(list(fit = object, coefficients = tab, ic = ic),
            class = "summary.bkwgfit")
}

print.summary.bkwgfit <- function(x, digits = 4, ...) {
  f <- x$fit
  cat("\n", .model_label(f), "\n", sep = "")
  cat("Method: ", .method_label(f$method), ";  n = ", f$n, "\n\n", sep = "")
  if (nrow(x$coefficients)) {
    cat("Coefficients:\n")
    stats::printCoefmat(x$coefficients, digits = digits, has.Pvalue = TRUE,
                        P.values = TRUE, na.print = "NA")
  }
  if (length(f$fixed))
    cat("Fixed: ", paste(names(f$fixed), "=", format(f$fixed), collapse = ", "),
        "\n", sep = "")
  cat("\n")
  print(signif(x$ic, digits + 2))
  if (f$method %in% c("mle", "mps") && f$k && all(is.na(f$se)))
    cat("(standard errors not available)\n")
  if (!f$method %in% c("mle", "mps"))
    cat("(standard errors are only computed for method = 'mle' or 'mps')\n")
  invisible(x)
}

coef.bkwgfit <- function(object, complete = FALSE, ...) {
  if (complete) object$estimate else object$coefficients
}

vcov.bkwgfit <- function(object, ...) object$vcov

logLik.bkwgfit <- function(object, ...) {
  structure(object$loglik, df = object$k, nobs = object$n, class = "logLik")
}

nobs.bkwgfit <- function(object, ...) object$n

confint.bkwgfit <- function(object, parm, level = 0.95,
                            type = c("wald", "log"), ...) {
  type <- match.arg(type)
  est <- object$coefficients
  se <- object$se
  if (missing(parm)) {
    parm <- names(est)
  } else if (is.numeric(parm)) {
    parm <- names(est)[parm]
  }
  z <- stats::qnorm(1 - (1 - level) / 2)
  ci <- cbind(est[parm] - z * se[parm], est[parm] + z * se[parm])
  if (type == "log") {
    lo0 <- c(a = 0, b = 0, alpha = 0, beta = 0, object$baseline$lower)[parm]
    pos <- !is.na(lo0) & lo0 == 0
    if (any(pos)) {
      e <- est[parm][pos]
      w <- z * se[parm][pos] / e
      ci[pos, ] <- cbind(e * exp(-w), e * exp(w))
    }
  }
  pr <- c((1 - level) / 2, 1 - (1 - level) / 2)
  dimnames(ci) <- list(parm, paste(format(100 * pr, trim = TRUE,
                                          scientific = FALSE, digits = 3), "%"))
  ci
}

residuals.bkwgfit <- function(object, type = c("quantile", "pit", "cox-snell"),
                              ...) {
  type <- match.arg(type)
  u <- do.call(pbkwg, c(list(q = object$x), .fit_args(object)))
  switch(type,
         pit = u,
         quantile = stats::qnorm(u),
         `cox-snell` = -log1p(-u))
}

quantile.bkwgfit <- function(x, probs = c(0.25, 0.5, 0.75), ...) {
  q <- do.call(qbkwg, c(list(p = probs), .fit_args(x)))
  names(q) <- paste0(format(100 * probs, trim = TRUE), "%")
  q
}

plot.bkwgfit <- function(x, y, which = 1:4, ...) {
  args <- .fit_args(x)
  d <- function(t) do.call(dbkwg, c(list(x = t), args))
  P <- function(t) do.call(pbkwg, c(list(q = t), args))
  Q <- function(u) do.call(qbkwg, c(list(p = u), args))
  xs <- sort(x$x)
  n <- length(xs)
  pp <- (seq_len(n) - 0.5) / n
  grid <- seq(min(xs), max(xs), length.out = 400)
  which <- intersect(which, 1:4)
  if (length(which) > 1L) {
    op <- graphics::par(mfrow = c(ceiling(length(which) / 2), 2))
    on.exit(graphics::par(op))
  }
  if (1 %in% which) {
    h <- graphics::hist(xs, plot = FALSE)
    dg <- d(grid)
    ym <- max(c(h$density, dg[is.finite(dg)]))
    graphics::hist(xs, freq = FALSE, ylim = c(0, ym), col = "grey90",
                   border = "grey60", main = "Histogram and fitted density",
                   xlab = "x")
    graphics::lines(grid, dg, lwd = 2, col = "steelblue")
  }
  if (2 %in% which) {
    graphics::plot(stats::ecdf(xs), main = "Empirical and fitted CDF",
                   xlab = "x", ylab = "F(x)", col = "grey40")
    graphics::lines(grid, P(grid), lwd = 2, col = "steelblue")
  }
  if (3 %in% which) {
    graphics::plot(pp, P(xs), xlab = "Empirical probabilities",
                   ylab = "Fitted probabilities", main = "P-P plot",
                   pch = 19, cex = 0.6)
    graphics::abline(0, 1, col = "red")
  }
  if (4 %in% which) {
    graphics::plot(Q(pp), xs, xlab = "Fitted quantiles",
                   ylab = "Sample quantiles", main = "Q-Q plot",
                   pch = 19, cex = 0.6)
    graphics::abline(0, 1, col = "red")
  }
  invisible(x)
}
