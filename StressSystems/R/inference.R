#' Estimate R = P(cY < X) from two samples with a confidence interval
#'
#' Fits the strength and stress distributions with the chosen method and
#' returns the plug-in estimate of \eqn{R}. For \code{method = "bayes"} the
#' estimate is the posterior mean of \eqn{R}, computed draw by draw.
#'
#' Interval types:
#' \describe{
#'   \item{\code{"delta"}}{delta method from the MLE observed information,
#'     built on the logit scale so the interval stays inside (0, 1).
#'     Requires \code{method = "mle"}.}
#'   \item{\code{"boot"}}{parametric bootstrap percentile interval; the
#'     re-estimation uses \code{method} if it is a point method, otherwise
#'     \code{base}.}
#'   \item{\code{"bayes"}}{equal-tail posterior interval (the HPD interval is
#'     also returned). Requires \code{method = "bayes"}.}
#' }
#' \code{"auto"} chooses \code{"delta"} for MLE, \code{"bayes"} for Bayes and
#' no interval otherwise.
#'
#' @param x Strength sample.
#' @param y Stress sample.
#' @param strength,stress Objects from \code{\link{ss_dist}}.
#' @param method Estimation method (see \code{\link{ss_fit}}).
#' @param ci Interval type: \code{"auto"}, \code{"delta"}, \code{"boot"},
#'   \code{"bayes"} or \code{"none"}.
#' @param level Confidence level.
#' @param c Stress scaling factor.
#' @param B Bootstrap samples for \code{ci = "boot"}.
#' @param base Point method inside jackknife and bootstrap estimators.
#' @param ndraw Maximum number of posterior draws used to compute R.
#' @param ... Passed to \code{\link{ss_fit}} (e.g. \code{iter}, \code{burnin},
#'   \code{prior} for Bayes).
#' @return An object of class \code{"ss_R_est"}: a list with the estimate,
#'   its standard error \code{se}, the interval limits \code{lower} and
#'   \code{upper}, and the two fitted models \code{fit_x} and \code{fit_y}.
#' @examples
#' set.seed(1)
#' X <- ss_dist("toppleone")
#' Y <- ss_dist("genrayleigh")
#' x <- rtoppleone(30, shape = 0.7)
#' y <- rgenrayleigh(30, shape = 1, rate = 2)
#' ss_estimate_R(x, y, X, Y, method = "mle")
#' \donttest{
#' ss_estimate_R(x, y, X, Y, method = "mle", ci = "boot", B = 200)
#' ss_estimate_R(x, y, X, Y, method = "bayes", iter = 3000, burnin = 1000)
#' }
#' @export
ss_estimate_R <- function(x, y, strength, stress, method = "mle",
                          ci = c("auto", "delta", "boot", "bayes", "none"),
                          level = 0.95, c = 1, B = 500, base = "mle",
                          ndraw = 1000, ...) {
  ci <- match.arg(ci)
  if (ci == "auto") ci <- switch(method, mle = "delta", bayes = "bayes", "none")
  fx <- ss_fit(x, strength, method = method, base = base, ...)
  fy <- ss_fit(y, stress, method = method, base = base, ...)
  alpha <- 1 - level
  out <- list(estimate = NA_real_, lower = NA_real_, upper = NA_real_, se = NA_real_,
              method = method, ci = ci, level = level, c = c,
              fit_x = fx, fit_y = fy)

  if (method == "bayes") {
    dx <- fx$draws
    dy <- fy$draws
    m <- min(ndraw, nrow(dx), nrow(dy))
    ix <- round(seq(1, nrow(dx), length.out = m))
    iy <- round(seq(1, nrow(dy), length.out = m))
    Rd <- vapply(seq_len(m), function(j)
      ss_R(strength, stress, .row(dx, ix[j]), .row(dy, iy[j]), c = c), numeric(1))
    out$estimate <- mean(Rd)
    out$se <- sd(Rd)
    out$draws <- Rd
    out$hpd <- .hpd(Rd, level)
    if (ci == "bayes") {
      q <- quantile(Rd, c(alpha / 2, 1 - alpha / 2), names = FALSE)
      out$lower <- q[1]
      out$upper <- q[2]
    }
  } else {
    out$estimate <- ss_R(strength, stress, fx$estimate, fy$estimate, c = c)
  }

  if (ci == "delta") {
    if (method != "mle" || is.null(fx$vcov_free) || is.null(fy$vcov_free)) {
      warning("delta-method interval needs method = 'mle' with an invertible information matrix")
    } else {
      g <- .grad_R(strength, stress, fx$eta, fy$eta, c)
      v <- sum(g$x * (fx$vcov_free %*% g$x)) + sum(g$y * (fy$vcov_free %*% g$y))
      R <- out$estimate
      out$se <- sqrt(max(v, 0))
      if (R > 0 && R < 1) {
        z <- qnorm(1 - alpha / 2)
        se_l <- out$se / (R * (1 - R))
        out$lower <- plogis(qlogis(R) - z * se_l)
        out$upper <- plogis(qlogis(R) + z * se_l)
      }
    }
  } else if (ci == "boot") {
    bm <- if (method %in% .point_methods) method else base
    px <- fx$estimate
    nx <- length(.check_data(x, strength))
    ny <- length(.check_data(y, stress))
    Rb <- vapply(seq_len(B), function(b) tryCatch({
      xb <- .rand(strength, nx, px)
      yb <- .rand(stress, ny, fy$estimate)
      ex <- .fit_point(xb, strength, bm, start = px, screen = FALSE)$estimate
      ey <- .fit_point(yb, stress, bm, start = fy$estimate, screen = FALSE)$estimate
      ss_R(strength, stress, ex, ey, c = c)
    }, error = function(e) NA_real_), numeric(1))
    Rb <- Rb[is.finite(Rb)]
    q <- quantile(Rb, c(alpha / 2, 1 - alpha / 2), names = FALSE)
    out$lower <- q[1]
    out$upper <- q[2]
    out$se <- sd(Rb)
    out$boot <- Rb
  } else if (ci == "bayes" && method != "bayes") {
    warning("posterior interval needs method = 'bayes'")
  }
  structure(out, class = "ss_R_est")
}

#' @export
print.ss_R_est <- function(x, ...) {
  cat("R = P(cY < X), c =", x$c, " method:", x$method, "\n")
  cat("Estimate:", format(x$estimate, digits = 6))
  if (is.finite(x$se)) cat("  SE:", format(x$se, digits = 4))
  cat("\n")
  if (is.finite(x$lower))
    cat(100 * x$level, "% ", x$ci, " interval: (", format(x$lower, digits = 5), ", ",
        format(x$upper, digits = 5), ")\n", sep = "")
  invisible(x)
}

.grad_R <- function(strength, stress, ex, ey, c, h = 1e-4) {
  f <- function(a, b) ss_R(strength, stress, .from_free(a, strength),
                           .from_free(b, stress), c = c, rel.tol = 1e-10)
  gr <- function(e, which) vapply(seq_along(e), function(j) {
    e1 <- e2 <- e
    e1[j] <- e1[j] + h
    e2[j] <- e2[j] - h
    if (which == "x") (f(e1, ey) - f(e2, ey)) / (2 * h)
    else (f(ex, e1) - f(ex, e2)) / (2 * h)
  }, numeric(1))
  list(x = gr(ex, "x"), y = gr(ey, "y"))
}

.hpd <- function(d, level) {
  d <- sort(d[is.finite(d)])
  n <- length(d)
  m <- max(1, floor(level * n))
  w <- d[m:n] - d[1:(n - m + 1)]
  i <- which.min(w)
  c(lower = d[i], upper = d[i + m - 1])
}
