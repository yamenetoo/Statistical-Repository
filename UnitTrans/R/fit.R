# ---- internal helpers ------------------------------------------------------

# Run `expr` with a given seed and restore the caller's RNG state afterwards.
.with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  genv <- globalenv()
  had <- exists(".Random.seed", envir = genv, inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = genv, inherits = FALSE)
  on.exit({
    if (had) assign(".Random.seed", old, envir = genv)
    else if (exists(".Random.seed", envir = genv, inherits = FALSE))
      rm(".Random.seed", envir = genv)
  }, add = TRUE)
  set.seed(seed)
  expr
}

.default_start <- function(dist) {
  nm <- .as_base(dist)$params
  s <- ifelse(nm %in% .location_names, 0, 1)
  names(s) <- nm
  s
}

.default_bounds <- function(nm) {
  lo <- ifelse(nm %in% .location_names, -Inf, 1e-8)
  names(lo) <- nm
  lo
}

# Map parameters to an unconstrained scale and back, according to the bounds.
.to_u <- function(x, lo, up) {
  u <- as.numeric(x)
  for (i in seq_along(u)) {
    u[i] <- if (is.finite(lo[i]) && is.finite(up[i]))
      stats::qlogis((x[i] - lo[i]) / (up[i] - lo[i]))
    else if (is.finite(lo[i])) log(x[i] - lo[i])
    else if (is.finite(up[i])) -log(up[i] - x[i])
    else x[i]
  }
  u
}

.from_u <- function(u, lo, up) {
  x <- as.numeric(u)
  for (i in seq_along(x)) {
    x[i] <- if (is.finite(lo[i]) && is.finite(up[i]))
      lo[i] + (up[i] - lo[i]) * stats::plogis(u[i])
    else if (is.finite(lo[i])) lo[i] + exp(u[i])
    else if (is.finite(up[i])) up[i] - exp(-u[i])
    else u[i]
  }
  x
}

.clamp_inside <- function(s, lo, up) {
  eps <- 1e-6
  s <- ifelse(is.finite(lo), pmax(s, lo + eps), s)
  ifelse(is.finite(up), pmin(s, up - eps), s)
}

# Starting values by matching sample and model quantiles.
.qm_start <- function(data, bl, tr, nm, base, lower, upper) {
  pr <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  eq <- stats::quantile(data, pr, names = FALSE)
  k <- length(nm)
  f <- function(u) {
    th <- .from_u(u, lower, upper); names(th) <- nm
    qv <- suppressWarnings(tryCatch(qtrans(pr, bl, th, tr), error = function(e) NULL))
    if (is.null(qv) || any(!is.finite(qv))) return(1e10)
    sum((qv - eq)^2)
  }
  u0 <- .to_u(.clamp_inside(base, lower, upper), lower, upper)
  o <- suppressWarnings(tryCatch(
    stats::optim(u0, f, method = if (k > 1) "Nelder-Mead" else "BFGS",
                 control = list(maxit = 500)),
    error = function(e) NULL))
  if (is.null(o) || !is.finite(o$value) || o$value >= 1e10) return(NULL)
  x <- .from_u(o$par, lower, upper)
  names(x) <- nm
  x
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

# Core optimiser. `data` is already validated; `dist` and `trans` are resolved.
.fit_core <- function(data, bl, tr, start = NULL, lower = NULL, upper = NULL,
                      n_start = 3, hessian = TRUE, check_support = TRUE,
                      adjusted = FALSE) {
  user_start <- !is.null(start)
  start <- if (user_start) .name_par(start, bl) else .default_start(bl)
  nm <- names(start)
  k <- length(start)
  lower <- if (is.null(lower)) unname(.default_bounds(nm)[nm]) else rep_len(unname(lower), k)
  upper <- if (is.null(upper)) rep(Inf, k) else rep_len(unname(upper), k)
  base <- .clamp_inside(unname(start), lower, upper)

  if (check_support) {
    dom <- check_domain(bl, stats::setNames(base, nm), tr)
    if (!dom$valid)
      stop("Incompatible baseline/transformation pair: ", dom$message, call. = FALSE)
  }

  nll <- function(th) {
    names(th) <- nm
    v <- -suppressWarnings(loglik_trans(th, data, bl, tr))
    if (!is.finite(v)) 1e10 else v
  }
  nllu <- function(u) nll(.from_u(u, lower, upper))

  starts <- list(base)
  if (!user_start && n_start >= 2) {
    qs <- .qm_start(data, bl, tr, nm, base, lower, upper)
    if (!is.null(qs)) starts[[2]] <- .clamp_inside(unname(qs), lower, upper)
  }
  while (length(starts) < n_start) {
    s <- base * exp(stats::rnorm(k, 0, 0.5)) + stats::rnorm(k, 0, 0.2)
    starts[[length(starts) + 1]] <- .clamp_inside(s, lower, upper)
  }

  best <- NULL
  for (s in starts) {
    u0 <- .to_u(s, lower, upper)
    o <- tryCatch(stats::optim(u0, nllu, method = "BFGS",
                               control = list(maxit = 500, reltol = 1e-12,
                                              ndeps = rep(1e-5, k))),
                  error = function(e) NULL)
    if (is.null(o) || o$value >= 1e10 || o$convergence != 0) {
      o2 <- suppressWarnings(tryCatch(
        stats::optim(if (is.null(o) || o$value >= 1e10) u0 else o$par, nllu,
                     method = "Nelder-Mead",
                     control = list(maxit = 2000, reltol = 1e-12)),
        error = function(e) NULL))
      if (!is.null(o2) && (is.null(o) || o2$value < o$value)) o <- o2
    }
    if (!is.null(o) && (is.null(best) || o$value < best$value)) best <- o
  }
  if (is.null(best) || best$value >= 1e10)
    stop("Optimisation failed for ", bl$name, " + ", tr$name, call. = FALSE)

  est <- .from_u(best$par, lower, upper)
  names(est) <- nm
  vc <- matrix(NA_real_, k, k)
  hessian_ok <- FALSE
  if (hessian) {
    h <- 1e-4 * pmax(abs(est), 1e-2)
    h <- pmin(h, (est - lower) / 4, (upper - est) / 4)
    H <- tryCatch(stats::optimHess(est, nll, control = list(ndeps = h)),
                  error = function(e) NULL)
    if (!is.null(H) && all(is.finite(H))) {
      H <- (H + t(H)) / 2
      ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
      if (min(ev) > 1e-10 * max(abs(ev))) {
        vc <- solve(H)
        hessian_ok <- TRUE
      }
    }
  }
  dimnames(vc) <- list(nm, nm)
  se <- suppressWarnings(sqrt(diag(vc)))
  ll <- -best$value
  n <- length(data)
  aic <- 2 * k - 2 * ll
  dom_end <- check_domain(bl, est, tr)
  structure(list(coef = est, se = se, vcov = vc, logLik = ll, k = k, n = n,
                 AIC = aic, BIC = k * log(n) - 2 * ll,
                 AICc = if (n - k - 1 > 0) aic + 2 * k * (k + 1) / (n - k - 1) else NA_real_,
                 dist = bl$name, trans = tr$name, data = data,
                 convergence = best$convergence, hessian_ok = hessian_ok,
                 domain_ok = dom_end$valid, adjusted = adjusted,
                 lower = lower, upper = upper),
            class = "transfit")
}

#' Maximum likelihood fit of a transformed distribution
#'
#' Fits the distribution of \eqn{V = g(Y)} to data in (0,1) by maximising the
#' log-likelihood over the baseline parameters. The optimisation is carried
#' out on an unconstrained scale (logarithmic for positive parameters), starting
#' from the default values, from a quantile-matching solution and from random
#' perturbations; the best optimum is kept.
#'
#' @param data Numeric vector with values strictly inside (0,1).
#' @param dist Baseline distribution name.
#' @param trans Transformation name or object.
#' @param start Optional named starting values. Defaults to 0 for location-type
#'   parameters and 1 otherwise.
#' @param lower,upper Optional bounds for the parameters.
#' @param n_start Number of starting points (the default or user start, a
#'   quantile-matching start, plus random perturbations); the best optimum is
#'   kept.
#' @param seed Optional seed for the random perturbations. The caller's random
#'   number stream is restored afterwards.
#' @param adjust If `TRUE`, data are mapped into the open interval with the
#'   Smithson-Verkuilen correction `(y (n - 1) + 0.5) / n`, which allows
#'   observations equal to 0 or 1.
#' @param check_support If `TRUE` (default), stop when the baseline puts mass
#'   outside the domain of the transformation (see [check_domain()]).
#' @return An object of class `transfit` with `coef`, `se`, `vcov`, `logLik`,
#'   `AIC`, `BIC`, `AICc`, `n`, `convergence` and `hessian_ok` (whether the
#'   Hessian was positive definite, so that standard errors are available).
#' @examples
#' set.seed(1)
#' x <- rtrans(200, "norm", c(0.5, 1), "logit")
#' fit <- fit_trans(x, "norm", "logit")
#' fit
#' @export
fit_trans <- function(data, dist, trans, start = NULL, lower = NULL,
                      upper = NULL, n_start = 3, seed = NULL, adjust = FALSE,
                      check_support = TRUE) {
  tr <- get_transformation(trans)
  bl <- .as_base(dist)
  data <- data[!is.na(data)]
  if (!is.numeric(data) || length(data) < 3L)
    stop("'data' must be a numeric vector with at least 3 non-missing values.",
         call. = FALSE)
  if (adjust) {
    n <- length(data)
    if (any(data < 0 | data > 1))
      stop("With adjust = TRUE the data must lie in [0,1].", call. = FALSE)
    data <- (data * (n - 1) + 0.5) / n
  }
  if (any(data <= 0 | data >= 1))
    stop("All data must lie strictly inside (0,1). ",
         "Use adjust = TRUE to apply the Smithson-Verkuilen correction.",
         call. = FALSE)
  .with_seed(seed, .fit_core(data, bl, tr, start, lower, upper, n_start,
                             hessian = TRUE, check_support = check_support,
                             adjusted = adjust))
}

#' @export
print.transfit <- function(x, ...) {
  cat("Transformed distribution fit\n")
  cat("  baseline:", x$dist, "| transformation:", x$trans, "| n =", x$n, "\n\n")
  print(round(cbind(Estimate = x$coef, `Std.Err` = x$se), 5))
  cat(sprintf("\nlogLik = %.3f  AIC = %.3f  BIC = %.3f  AICc = %.3f\n",
              x$logLik, x$AIC, x$BIC, x$AICc))
  if (x$convergence != 0) cat("Warning: optimiser reported convergence code", x$convergence, "\n")
  if (!isTRUE(x$hessian_ok))
    cat("Note: the Hessian is not positive definite; standard errors are unavailable.\n")
  if (!isTRUE(x$domain_ok))
    cat("Warning: at the estimate the baseline puts mass outside the transformation domain.\n")
  if (isTRUE(x$adjusted)) cat("Note: data were adjusted with the Smithson-Verkuilen correction.\n")
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
nobs.transfit <- function(object, ...) object$n

#' Residuals and predictions of a fitted model
#'
#' `residuals(type = "quantile")` returns normalised quantile residuals
#' \eqn{\Phi^{-1}(F(v_i))}, which are approximately standard normal when the
#' model is adequate; `type = "pit"` returns the probability integral
#' transform \eqn{F(v_i)}, approximately uniform.
#' @param object A `transfit` object.
#' @param type `"quantile"` or `"pit"` for residuals.
#' @param newdata Numeric vector of values in (0,1) (or probabilities when
#'   `type = "quantile"`). Defaults to the data used for the fit.
#' @param ... Unused.
#' @return Numeric vector.
#' @examples
#' set.seed(1)
#' fit <- fit_trans(rtrans(150, "norm", c(0, 1), "logit"), "norm", "logit")
#' head(residuals(fit))
#' predict(fit, c(.2, .5), type = "cdf")
#' @name transfit_methods
#' @export
residuals.transfit <- function(object, type = c("quantile", "pit"), ...) {
  type <- match.arg(type)
  u <- ptrans(object$data, object$dist, object$coef, object$trans)
  if (type == "pit") u else stats::qnorm(u)
}

#' @rdname transfit_methods
#' @param type For `predict`: `"density"`, `"cdf"` or `"quantile"`.
#' @export
predict.transfit <- function(object, newdata = NULL,
                             type = c("density", "cdf", "quantile"), ...) {
  type <- match.arg(type)
  if (is.null(newdata)) newdata <- if (type == "quantile") stats::ppoints(object$n) else object$data
  switch(type,
         density = dtrans(newdata, object$dist, object$coef, object$trans),
         cdf = ptrans(newdata, object$dist, object$coef, object$trans),
         quantile = qtrans(newdata, object$dist, object$coef, object$trans))
}

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

# Profile-likelihood limits for parameter j.
.profile_ci <- function(object, j, level) {
  bl <- .as_base(object$dist); tr <- get_transformation(object$trans)
  data <- object$data; est <- object$coef; nm <- names(est); k <- length(est)
  lower <- object$lower; upper <- object$upper
  crit <- stats::qchisq(level, 1) / 2
  ll_hat <- object$logLik
  prof <- function(val) {
    if (k == 1L)
      return(suppressWarnings(loglik_trans(stats::setNames(val, nm), data, bl, tr)))
    g <- function(u) {
      th <- numeric(k); th[j] <- val
      th[-j] <- .from_u(u, lower[-j], upper[-j]); names(th) <- nm
      v <- -suppressWarnings(loglik_trans(th, data, bl, tr))
      if (is.finite(v)) v else 1e10
    }
    u0 <- .to_u(.clamp_inside(est[-j], lower[-j], upper[-j]), lower[-j], upper[-j])
    o <- suppressWarnings(stats::optim(u0, g,
           method = if (k - 1L == 1L) "BFGS" else "Nelder-Mead",
           control = list(maxit = 500, reltol = 1e-10)))
    -o$value
  }
  f <- function(val) ll_hat - prof(val) - crit
  se_j <- object$se[j]
  step0 <- if (is.finite(se_j) && se_j > 0) se_j else 0.1 * abs(est[j]) + 0.1
  find <- function(dir) {
    x0 <- est[[j]]; step <- step0
    bound <- if (dir < 0) lower[j] else upper[j]
    for (i in seq_len(60)) {
      x1 <- x0 + dir * step
      if (is.finite(bound) && dir * (x1 - bound) >= 0) x1 <- x0 + dir * abs(bound - x0) / 2
      f1 <- f(x1)
      if (is.finite(f1) && f1 >= 0)
        return(stats::uniroot(f, sort(c(x0, x1)), tol = 1e-8 * max(1, abs(x0)))$root)
      if (abs(x1 - x0) < 1e-12) break
      x0 <- x1; step <- step * 2
    }
    NA_real_
  }
  c(find(-1), find(1))
}

#' Confidence intervals for a fitted model
#'
#' @param object A `transfit` object.
#' @param parm Parameter names or indices (default: all).
#' @param level Confidence level.
#' @param method `"wald"` (from the Hessian), `"bootstrap"` (parametric
#'   percentile bootstrap) or `"profile"` (profile likelihood, usually
#'   better for positive or skewed parameters). A limit returned as `NA` means
#'   that the profile did not cross the cut-off in the searched range.
#' @param B Bootstrap replications.
#' @param seed Optional seed for the bootstrap (the caller's random number
#'   stream is restored afterwards).
#' @param ... Unused.
#' @return Matrix with lower and upper limits. For the bootstrap the number of
#'   failed replications is stored in the attribute `"failed"`.
#' @examples
#' set.seed(1)
#' fit <- fit_trans(rtrans(150, "norm", c(0, 1), "logit"), "norm", "logit")
#' confint(fit)
#' confint(fit, method = "profile")
#' @export
confint.transfit <- function(object, parm, level = 0.95,
                             method = c("wald", "bootstrap", "profile"),
                             B = 200, seed = NULL, ...) {
  method <- match.arg(method)
  nm <- names(object$coef)
  idx <- if (missing(parm)) seq_along(nm) else if (is.character(parm)) match(parm, nm) else parm
  if (anyNA(idx)) stop("Unknown parameter in 'parm'.", call. = FALSE)
  a <- (1 - level) / 2
  failed <- NULL
  if (method == "wald") {
    z <- stats::qnorm(1 - a)
    ci <- cbind(object$coef - z * object$se, object$coef + z * object$se)[idx, , drop = FALSE]
  } else if (method == "profile") {
    ci <- t(vapply(idx, function(j) .profile_ci(object, j, level), numeric(2)))
  } else {
    bl <- .as_base(object$dist); tr <- get_transformation(object$trans)
    bs <- .with_seed(seed, {
      m <- matrix(NA_real_, B, length(nm))
      for (b in seq_len(B)) {
        xs <- .clip01(rtrans(object$n, bl, object$coef, tr))
        f <- tryCatch(.fit_core(xs, bl, tr, start = object$coef, lower = object$lower,
                                upper = object$upper, n_start = 1, hessian = FALSE,
                                check_support = FALSE)$coef,
                      error = function(e) NULL)
        if (!is.null(f)) m[b, ] <- f
      }
      m
    })
    failed <- sum(!stats::complete.cases(bs))
    if (failed > 0)
      warning(failed, " of ", B, " bootstrap replications failed and were dropped.",
              call. = FALSE)
    ci <- t(apply(bs[, idx, drop = FALSE], 2, stats::quantile,
                  probs = c(a, 1 - a), na.rm = TRUE))
  }
  dimnames(ci) <- list(nm[idx], paste0(format(100 * c(a, 1 - a), trim = TRUE), " %"))
  if (!is.null(failed)) attr(ci, "failed") <- failed
  ci
}

#' Plot a fitted model
#'
#' Histogram with fitted density, plus fitted vs empirical CDF, Q-Q and P-P.
#' @param x A `transfit` object.
#' @param which Which of the four panels to draw (subset of 1:4).
#' @param ... Unused.
#' @return Invisibly `x`.
#' @examples
#' set.seed(1)
#' fit <- fit_trans(rtrans(150, "norm", c(0, 1), "logit"), "norm", "logit")
#' plot(fit)
#' @export
plot.transfit <- function(x, which = 1:4, ...) {
  which <- intersect(1:4, which)
  if (length(which) > 1L) {
    op <- graphics::par(mfrow = c(ceiling(length(which) / 2), min(2, length(which))))
    on.exit(graphics::par(op))
  }
  d <- x$data; n <- length(d); xs <- sort(d)
  grid <- seq(min(d) / 2, min(1 - 1e-4, max(d) + (1 - max(d)) / 2), length.out = 400)
  pp <- (seq_len(n) - 0.5) / n
  if (1 %in% which) {
    graphics::hist(d, probability = TRUE, col = "grey90", main = "Density",
                   xlab = "v", breaks = "FD")
    graphics::lines(grid, dtrans(grid, x$dist, x$coef, x$trans), lwd = 2, col = "darkgreen")
  }
  if (2 %in% which) {
    plot(stats::ecdf(d), main = "CDF", xlab = "v", ylab = "F(v)")
    graphics::lines(grid, ptrans(grid, x$dist, x$coef, x$trans), lwd = 2, col = "darkgreen")
  }
  if (3 %in% which) {
    plot(qtrans(pp, x$dist, x$coef, x$trans), xs, main = "Q-Q",
         xlab = "Theoretical", ylab = "Empirical"); graphics::abline(0, 1, col = "red")
  }
  if (4 %in% which) {
    plot(pp, ptrans(xs, x$dist, x$coef, x$trans), main = "P-P",
         xlab = "Empirical", ylab = "Theoretical"); graphics::abline(0, 1, col = "red")
  }
  invisible(x)
}

#' Fit many baseline/transformation combinations and rank them
#'
#' Fits every valid combination of baseline distributions and transformations
#' and ranks them by an information criterion. Combinations whose baseline
#' support is incompatible with the transformation, whose optimiser did not
#' converge or whose fit puts mass outside the transformation domain are
#' dropped (reported when `verbose = TRUE`).
#'
#' @param data Numeric vector in (0,1).
#' @param dists Character vector of baseline distribution names.
#' @param trans Character vector of transformation names (default: all).
#' @param criterion `"AIC"`, `"BIC"` or `"AICc"`.
#' @param verbose Print skipped combinations.
#' @param drop_nonconverged Drop fits whose optimiser reported non-convergence.
#' @param seed Optional seed passed to every [fit_trans()] call (the same seed
#'   for each combination, so results do not depend on the order of the
#'   combinations or on `cores`).
#' @param cores Number of processes (Unix-alike systems only; ignored on Windows).
#' @param adjust Passed to [fit_trans()].
#' @return An object of class `transfit_all`: a ranking `table` (with
#'   `delta` and Akaike `weight` for the chosen criterion) and the list of
#'   `fits`.
#' @examples
#' set.seed(2)
#' x <- rtrans(200, "gamma", c(2, 1), "rational")
#' res <- fit_all(x, c("gamma", "lnorm", "weibull"), c("rational", "exponential"))
#' res$table
#' best_model(res)
#' @export
fit_all <- function(data, dists, trans = names(.registry()),
                    criterion = c("AIC", "BIC", "AICc"), verbose = FALSE,
                    drop_nonconverged = TRUE, seed = NULL, cores = 1L,
                    adjust = FALSE) {
  criterion <- match.arg(criterion)
  combos <- expand.grid(d = dists, t = trans, stringsAsFactors = FALSE)
  keys <- paste(combos$d, combos$t, sep = " + ")
  one <- function(i) {
    tryCatch(fit_trans(data, combos$d[i], combos$t[i], seed = seed, adjust = adjust),
             error = function(e) structure(conditionMessage(e), class = "fit_failure"))
  }
  res <- if (cores > 1L && .Platform$OS.type != "windows")
    parallel::mclapply(seq_len(nrow(combos)), one, mc.cores = cores)
  else lapply(seq_len(nrow(combos)), one)
  fits <- list(); rows <- list()
  for (i in seq_along(res)) {
    f <- res[[i]]; key <- keys[i]
    if (!inherits(f, "transfit")) {
      if (verbose) message("skip ", key, ": ",
                           if (inherits(f, "fit_failure")) unclass(f) else "failed")
      next
    }
    if (drop_nonconverged && (f$convergence != 0 || !isTRUE(f$domain_ok))) {
      if (verbose) message("skip ", key, ": not converged or invalid at the estimate")
      next
    }
    fits[[key]] <- f
    rows[[key]] <- data.frame(model = key, baseline = combos$d[i],
                              transformation = combos$t[i], k = f$k,
                              logLik = f$logLik, AIC = f$AIC, BIC = f$BIC,
                              AICc = f$AICc, converged = f$convergence == 0,
                              stringsAsFactors = FALSE)
  }
  if (!length(rows)) stop("No valid combination could be fitted.", call. = FALSE)
  tab <- do.call(rbind, rows)
  tab <- tab[order(tab[[criterion]]), ]
  if (all(is.na(tab[[criterion]]))) stop("Criterion ", criterion, " is not available.", call. = FALSE)
  tab$delta <- tab[[criterion]] - min(tab[[criterion]], na.rm = TRUE)
  w <- exp(-tab$delta / 2)
  tab$weight <- w / sum(w, na.rm = TRUE)
  rownames(tab) <- NULL
  structure(list(table = tab, fits = fits, criterion = criterion),
            class = "transfit_all")
}

#' @export
print.transfit_all <- function(x, ...) {
  cat("Ranking by", x$criterion, "\n")
  print(x$table[, c("model", "k", "logLik", "AIC", "BIC", "AICc", "delta", "weight")],
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
