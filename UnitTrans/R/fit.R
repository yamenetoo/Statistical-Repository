.default_start <- function(dist) {
  nm <- baseline_params(dist)
  s <- ifelse(nm %in% .location_names, 0, 1)
  names(s) <- nm
  s
}

.default_bounds <- function(nm) {
  lo <- ifelse(nm %in% .location_names, -Inf, 1e-8)
  names(lo) <- nm
  lo
}

#' Log-likelihood of the transformed distribution
#'
#' @param par Baseline parameters.
#' @param data Observations in (0,1).
#' @param dist Baseline distribution name.
#' @param trans Transformation name or object.
#' @return Numeric scalar.
#' @examples
#' loglik_trans(c(0, 1), c(.2, .4, .7), "norm", "logit")
#' @export
loglik_trans <- function(par, data, dist, trans) {
  ll <- dtrans(data, dist, par, trans, log = TRUE)
  if (anyNA(ll) || any(!is.finite(ll))) return(-Inf)
  sum(ll)
}

#' Maximum likelihood fit of a transformed distribution
#'
#' Fits the distribution of \eqn{V = g(Y)} to data in (0,1) by maximising the
#' log-likelihood over the baseline parameters.
#'
#' @param data Numeric vector with values strictly inside (0,1).
#' @param dist Baseline distribution name.
#' @param trans Transformation name or object.
#' @param start Optional named starting values. Defaults to 0 for location-type
#'   parameters and 1 otherwise.
#' @param lower,upper Optional bounds for the parameters.
#' @param n_start Number of starting points (the default plus random
#'   perturbations); the best optimum is kept.
#' @param seed Optional seed for the random perturbations.
#' @return An object of class `transfit` with `coef`, `se`, `vcov`, `logLik`,
#'   `AIC`, `BIC`, `AICc`, `n`, `convergence`.
#' @examples
#' set.seed(1)
#' x <- rtrans(200, "norm", c(0.5, 1), "logit")
#' fit <- fit_trans(x, "norm", "logit")
#' fit
#' @export
fit_trans <- function(data, dist, trans, start = NULL, lower = NULL,
                      upper = NULL, n_start = 3, seed = NULL) {
  tr <- get_transformation(trans)
  data <- data[!is.na(data)]
  if (any(data <= 0 | data >= 1))
    stop("All data must lie strictly inside (0,1).", call. = FALSE)
  if (is.null(start)) start <- .default_start(dist) else start <- .name_par(start, dist)
  nm <- names(start)
  k <- length(start)
  if (is.null(lower)) lower <- .default_bounds(nm)[nm] else lower <- rep_len(lower, k)
  if (is.null(upper)) upper <- rep(Inf, k) else upper <- rep_len(upper, k)
  lower <- unname(lower); upper <- unname(upper)
  nll <- function(th) {
    names(th) <- nm
    v <- -loglik_trans(th, data, dist, tr)
    if (!is.finite(v)) 1e10 else v
  }
  if (!is.null(seed)) set.seed(seed)
  starts <- list(unname(start))
  if (n_start > 1)
    for (i in seq_len(n_start - 1)) {
      s <- unname(start) * exp(stats::rnorm(k, 0, 0.5)) + stats::rnorm(k, 0, 0.2)
      starts[[length(starts) + 1]] <- pmin(pmax(s, lower + 1e-6), upper)
    }
  best <- NULL
  for (s in starts) {
    s <- pmin(pmax(s, ifelse(is.finite(lower), lower + 1e-6, -Inf)), upper)
    o <- tryCatch(stats::optim(s, nll, method = "L-BFGS-B", lower = lower,
                               upper = upper, control = list(maxit = 1000)),
                  error = function(e) NULL)
    if (!is.null(o) && (is.null(best) || o$value < best$value)) best <- o
  }
  if (is.null(best) || best$value >= 1e10)
    stop("Optimisation failed for ", dist, " + ", tr$name, call. = FALSE)
  est <- best$par; names(est) <- nm
  H <- tryCatch(stats::optimHess(est, nll), error = function(e) NULL)
  vc <- tryCatch(solve(H), error = function(e) NULL)
  if (is.null(vc)) { vc <- matrix(NA_real_, k, k) }
  dimnames(vc) <- list(nm, nm)
  se <- suppressWarnings(sqrt(diag(vc)))
  ll <- -best$value
  n <- length(data)
  aic <- 2 * k - 2 * ll
  structure(list(coef = est, se = se, vcov = vc, logLik = ll, k = k, n = n,
                 AIC = aic, BIC = k * log(n) - 2 * ll,
                 AICc = if (n - k - 1 > 0) aic + 2 * k * (k + 1) / (n - k - 1) else NA_real_,
                 dist = dist, trans = tr$name, data = data,
                 convergence = best$convergence),
            class = "transfit")
}

#' @export
print.transfit <- function(x, ...) {
  cat("Transformed distribution fit\n")
  cat("  baseline:", x$dist, "| transformation:", x$trans, "| n =", x$n, "\n\n")
  print(round(cbind(Estimate = x$coef, `Std.Err` = x$se), 5))
  cat(sprintf("\nlogLik = %.3f  AIC = %.3f  BIC = %.3f  AICc = %.3f\n",
              x$logLik, x$AIC, x$BIC, x$AICc))
  if (x$convergence != 0) cat("Warning: optimiser reported convergence code", x$convergence, "\n")
  invisible(x)
}

#' @export
coef.transfit <- function(object, ...) object$coef
#' @export
vcov.transfit <- function(object, ...) object$vcov
#' @export
logLik.transfit <- function(object, ...)
  structure(object$logLik, df = object$k, nobs = object$n, class = "logLik")

#' @export
summary.transfit <- function(object, ...) {
  z <- object$coef / object$se
  tab <- cbind(Estimate = object$coef, `Std.Err` = object$se, z = z,
               `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  structure(list(table = tab, fit = object), class = "summary.transfit")
}

#' @export
print.summary.transfit <- function(x, ...) {
  print(x$fit)
  cat("\nWald tests:\n")
  print(round(x$table, 5))
  invisible(x)
}

#' Confidence intervals for a fitted model
#'
#' @param object A `transfit` object.
#' @param parm Ignored (all parameters).
#' @param level Confidence level.
#' @param method `"wald"` (from the Hessian) or `"bootstrap"` (parametric
#'   percentile bootstrap).
#' @param B Bootstrap replications.
#' @param ... Unused.
#' @return Matrix with lower and upper limits.
#' @examples
#' set.seed(1)
#' fit <- fit_trans(rtrans(150, "norm", c(0, 1), "logit"), "norm", "logit")
#' confint(fit)
#' @export
confint.transfit <- function(object, parm, level = 0.95, method = c("wald", "bootstrap"),
                             B = 200, ...) {
  method <- match.arg(method)
  a <- (1 - level) / 2
  if (method == "wald") {
    z <- stats::qnorm(1 - a)
    ci <- cbind(object$coef - z * object$se, object$coef + z * object$se)
  } else {
    bs <- t(replicate(B, {
      xs <- rtrans(object$n, object$dist, object$coef, object$trans)
      tryCatch(fit_trans(xs, object$dist, object$trans, start = object$coef,
                         n_start = 1)$coef,
               error = function(e) rep(NA_real_, length(object$coef)))
    }))
    ci <- t(apply(bs, 2, stats::quantile, probs = c(a, 1 - a), na.rm = TRUE))
  }
  dimnames(ci) <- list(names(object$coef), paste0(format(100 * c(a, 1 - a), trim = TRUE), " %"))
  ci
}

#' Plot a fitted model
#'
#' Histogram with fitted density, plus fitted vs empirical CDF, Q-Q and P-P.
#' @param x A `transfit` object.
#' @param ... Unused.
#' @return Invisibly `x`.
#' @examples
#' set.seed(1)
#' fit <- fit_trans(rtrans(150, "norm", c(0, 1), "logit"), "norm", "logit")
#' plot(fit)
#' @export
plot.transfit <- function(x, ...) {
  op <- graphics::par(mfrow = c(2, 2)); on.exit(graphics::par(op))
  d <- x$data; n <- length(d); xs <- sort(d)
  grid <- seq(min(d) / 2, min(1 - 1e-4, max(d) + (1 - max(d)) / 2), length.out = 400)
  graphics::hist(d, probability = TRUE, col = "grey90", main = "Density",
                 xlab = "v", breaks = "FD")
  graphics::lines(grid, dtrans(grid, x$dist, x$coef, x$trans), lwd = 2, col = "darkgreen")
  plot(stats::ecdf(d), main = "CDF", xlab = "v", ylab = "F(v)")
  graphics::lines(grid, ptrans(grid, x$dist, x$coef, x$trans), lwd = 2, col = "darkgreen")
  pp <- (seq_len(n) - 0.5) / n
  plot(qtrans(pp, x$dist, x$coef, x$trans), xs, main = "Q-Q",
       xlab = "Theoretical", ylab = "Empirical"); graphics::abline(0, 1, col = "red")
  plot(pp, ptrans(xs, x$dist, x$coef, x$trans), main = "P-P",
       xlab = "Empirical", ylab = "Theoretical"); graphics::abline(0, 1, col = "red")
  invisible(x)
}

#' Fit many baseline/transformation combinations and rank them
#'
#' Fits every valid combination of baseline distributions and transformations
#' and ranks them by an information criterion.
#'
#' @param data Numeric vector in (0,1).
#' @param dists Character vector of baseline distribution names.
#' @param trans Character vector of transformation names (default: all).
#' @param criterion `"AIC"`, `"BIC"` or `"AICc"`.
#' @param verbose Print skipped combinations.
#' @return An object of class `transfit_all`: a ranking `table` and the list of
#'   `fits`.
#' @examples
#' set.seed(2)
#' x <- rtrans(200, "gamma", c(2, 1), "rational")
#' res <- fit_all(x, c("gamma", "lnorm", "weibull"), c("rational", "exponential"))
#' res$table
#' best_model(res)
#' @export
fit_all <- function(data, dists, trans = names(.registry()),
                    criterion = c("AIC", "BIC", "AICc"), verbose = FALSE) {
  criterion <- match.arg(criterion)
  fits <- list(); rows <- list()
  for (d in dists) for (t in trans) {
    key <- paste(d, t, sep = " + ")
    f <- tryCatch({
      st <- .default_start(d)
      dom <- check_domain(d, st, t)
      if (!dom$valid) stop(dom$message)
      fit_trans(data, d, t)
    }, error = function(e) { if (verbose) message("skip ", key, ": ", conditionMessage(e)); NULL })
    if (is.null(f)) next
    fits[[key]] <- f
    rows[[key]] <- data.frame(model = key, baseline = d, transformation = t,
                              k = f$k, logLik = f$logLik, AIC = f$AIC,
                              BIC = f$BIC, AICc = f$AICc, stringsAsFactors = FALSE)
  }
  if (!length(rows)) stop("No valid combination could be fitted.", call. = FALSE)
  tab <- do.call(rbind, rows)
  tab <- tab[order(tab[[criterion]]), ]
  tab$delta <- tab[[criterion]] - min(tab[[criterion]])
  rownames(tab) <- NULL
  structure(list(table = tab, fits = fits, criterion = criterion),
            class = "transfit_all")
}

#' @export
print.transfit_all <- function(x, ...) {
  cat("Ranking by", x$criterion, "\n")
  print(x$table[, c("model", "k", "logLik", "AIC", "BIC", "AICc", "delta")],
        digits = 5, row.names = FALSE)
  invisible(x)
}

#' Extract the best model from a [fit_all()] result
#'
#' @param object A `transfit_all` object.
#' @param criterion Criterion to use; defaults to the one used for ranking.
#' @return The best `transfit` object.
#' @export
best_model <- function(object, criterion = object$criterion) {
  tab <- object$table
  object$fits[[tab$model[which.min(tab[[criterion]])]]]
}
