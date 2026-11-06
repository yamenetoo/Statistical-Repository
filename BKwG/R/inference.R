## Goodness of fit, hypothesis tests and model comparison ------------------------

## W* and A* (Chen and Balakrishnan, 1995) and the Kolmogorov-Smirnov statistic
.gof_stats <- function(x, Fn) {
  n <- length(x)
  xs <- sort(x)
  i <- seq_len(n)
  eps <- 1e-12
  u <- pmin(pmax(Fn(xs), eps), 1 - eps)
  D <- max(i / n - u, u - (i - 1) / n)
  y <- stats::qnorm(u)
  v <- stats::pnorm((y - mean(y)) / stats::sd(y))
  v <- pmin(pmax(v, eps), 1 - eps)
  W2 <- sum((v - (2 * i - 1) / (2 * n))^2) + 1 / (12 * n)
  A2 <- -n - mean((2 * i - 1) * (log(v) + log(1 - rev(v))))
  c(`W*` = W2 * (1 + 0.5 / n),
    `A*` = A2 * (1 + 0.75 / n + 2.25 / n^2),
    KS = D)
}

.fit_cdf <- function(fit) {
  args <- .fit_args(fit)
  function(q) do.call(pbkwg, c(list(q = q), args))
}

gof_bkwg <- function(fit, B = 0, seed = NULL) {
  if (!inherits(fit, "bkwgfit")) stop("'fit' must be a 'bkwgfit' object")
  Fn <- .fit_cdf(fit)
  obs <- .gof_stats(fit$x, Fn)
  ks_p <- suppressWarnings(stats::ks.test(fit$x, Fn)$p.value)
  pval <- c(`W*` = NA_real_, `A*` = NA_real_, KS = ks_p)
  Bdone <- 0L
  if (B > 0) {
    if (!is.null(seed)) set.seed(seed)
    args <- .fit_args(fit)
    n <- fit$n
    boot <- matrix(NA_real_, B, 3L)
    for (r in seq_len(B)) {
      xb <- do.call(rbkwg, c(list(n = n), args))
      fb <- tryCatch(suppressWarnings(
        fit_bkwg(xb, baseline = fit$baseline, submodel = fit$submodel,
                 fixed = fit$fixed_user, start = fit$coefficients,
                 method = fit$method, nstart = 1, hessian = FALSE)),
        error = function(e) NULL)
      if (!is.null(fb)) boot[r, ] <- .gof_stats(xb, .fit_cdf(fb))
    }
    ok <- stats::complete.cases(boot)
    Bdone <- sum(ok)
    if (Bdone) {
      pval <- vapply(1:3, function(j)
        (1 + sum(boot[ok, j] >= obs[j])) / (Bdone + 1), numeric(1))
      names(pval) <- names(obs)
    }
  }
  tab <- data.frame(statistic = obs, p.value = pval,
                    row.names = c("Cramer-von Mises (W*)",
                                  "Anderson-Darling (A*)",
                                  "Kolmogorov-Smirnov (KS)"))
  structure(list(table = tab, B = Bdone, fit = fit), class = "bkwggof")
}

print.bkwggof <- function(x, digits = 4, ...) {
  cat("\nGoodness of fit: ", .model_label(x$fit), "\n\n", sep = "")
  print(signif(x$table, digits))
  if (x$B > 0) {
    cat("\np-values from a parametric bootstrap with", x$B, "replicates.\n")
  } else {
    cat("\nKS p-value is asymptotic and ignores parameter estimation;",
        "use B > 0 for parametric bootstrap p-values.\n")
  }
  invisible(x)
}

## Likelihood ratio test of a nested sub-model
lrt_bkwg <- function(fit, restricted = "KwG") {
  if (!inherits(fit, "bkwgfit")) stop("'fit' must be a 'bkwgfit' object")
  if (fit$method != "mle")
    stop("the likelihood ratio test requires maximum likelihood fits")
  if (is.character(restricted)) {
    fit0 <- fit_bkwg(fit$x, baseline = fit$baseline, submodel = restricted,
                     fixed = fit$fixed_user, method = "mle", hessian = FALSE)
  } else if (inherits(restricted, "bkwgfit")) {
    fit0 <- restricted
    if (fit0$method != "mle")
      stop("the likelihood ratio test requires maximum likelihood fits")
    if (fit0$n != fit$n || any(fit0$x != fit$x))
      stop("both models must be fitted to the same data")
  } else {
    stop("'restricted' must be a sub-model name or a 'bkwgfit' object")
  }
  if (fit0$baseline$name != fit$baseline$name)
    stop("both models must use the same baseline")
  if (!all(names(fit$fixed) %in% names(fit0$fixed)) ||
      any(fit0$fixed[names(fit$fixed)] != fit$fixed))
    stop("the restricted model is not nested in the full model")
  df <- fit$k - fit0$k
  if (df <= 0) stop("the restricted model must have fewer free parameters")
  lr <- 2 * (fit$loglik - fit0$loglik)
  if (lr < -1e-6)
    warning("restricted log-likelihood exceeds the full one; the full fit ",
            "may not be at the global maximum (refit with more 'nstart')",
            call. = FALSE)
  lr <- max(lr, 0)
  sm <- bkwg_submodels()
  structure(list(
    statistic = c(LR = lr), parameter = c(df = df),
    p.value = stats::pchisq(lr, df, lower.tail = FALSE),
    method = "Likelihood ratio test for nested BKw-G models",
    data.name = paste0(fit$data.name, ";  H0: ",
                       sm$name[sm$submodel == fit0$submodel], "  vs  H1: ",
                       sm$name[sm$submodel == fit$submodel],
                       "  (baseline: ", fit$baseline$name, ")"),
    estimate = c(`logLik H1` = fit$loglik, `logLik H0` = fit0$loglik)),
    class = "htest")
}

## Wald test of H0: theta_j = theta_j0 for a set of free parameters
wald_bkwg <- function(fit, H0) {
  if (!inherits(fit, "bkwgfit")) stop("'fit' must be a 'bkwgfit' object")
  H0 <- unlist(H0)
  if (is.null(names(H0)) || any(!names(H0) %in% fit$free))
    stop("'H0' must be named with free parameters among: ",
         paste(fit$free, collapse = ", "))
  nm <- names(H0)
  V <- fit$vcov[nm, nm, drop = FALSE]
  if (any(!is.finite(V))) stop("the covariance matrix is not available")
  d <- fit$coefficients[nm] - H0
  W <- drop(t(d) %*% solve(V, d))
  structure(list(
    statistic = c(W = W), parameter = c(df = length(nm)),
    p.value = stats::pchisq(W, length(nm), lower.tail = FALSE),
    method = "Wald test",
    data.name = paste0(fit$data.name, ";  H0: ",
                       paste(nm, "=", H0, collapse = ", ")),
    estimate = fit$coefficients[nm]),
    class = "htest")
}

## Fit several sub-models (and baselines) and tabulate the usual criteria
compare_bkwg <- function(x, baseline = "weibull",
                         submodels = c("BKwG", "KwG", "BetaG", "EKwG", "BEG",
                                       "EG", "G"),
                         method = "mle", ...) {
  bl <- if (inherits(baseline, "bkwg_baseline")) list(baseline) else
    as.list(baseline)
  fits <- list()
  rows <- list()
  for (bb in bl) {
    B <- .get_baseline(bb)
    for (sm in submodels) {
      lab <- if (sm == "G") B$name else paste0(sm, "-", B$name)
      f <- tryCatch(fit_bkwg(x, baseline = B, submodel = sm, method = method,
                             ...),
                    error = function(e) {
                      warning("fit of ", lab, " failed: ", conditionMessage(e),
                              call. = FALSE)
                      NULL
                    })
      if (is.null(f)) next
      fits[[lab]] <- f
      g <- gof_bkwg(f)$table
      rows[[lab]] <- data.frame(
        model = lab, k = f$k, logLik = f$loglik, AIC = f$AIC, CAIC = f$CAIC,
        BIC = f$BIC, HQIC = f$HQIC, W = g$statistic[1], A = g$statistic[2],
        KS = g$statistic[3], KS.p = g$p.value[3], stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) stop("no model could be fitted")
  out <- do.call(rbind, rows)
  out <- out[order(out$AIC), , drop = FALSE]
  rownames(out) <- NULL
  names(out)[names(out) == "W"] <- "W*"
  names(out)[names(out) == "A"] <- "A*"
  attr(out, "fits") <- fits
  out
}
