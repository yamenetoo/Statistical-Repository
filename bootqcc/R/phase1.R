#' Build a chart from Phase I data, accounting for estimation uncertainty
#'
#' [bootqcc()] simulates from *known* parameters. In practice the parameters are
#' unknown and are estimated from a finite Phase I sample, and control limits
#' built from estimates have a random false-alarm rate: for a given Phase I
#' sample the chart can alarm much more (or less) often than the nominal
#' `alpha`, and the average run length is shorter than `1 / alpha`
#' (Quesenberry 1993; Jensen et al. 2006). `bootqcc_phase1()` takes the raw
#' Phase I data, estimates the parameters by maximum likelihood, builds the
#' bootstrap limits, and then uses a second (outer) parametric bootstrap
#' to measure the conditional false-alarm probability caused by the
#' estimation and, optionally, to widen the limits so that the target holds
#' (Gandy and Kvaloy 2013; Saleh et al. 2015).
#'
#' @details
#' Let \eqn{\hat\theta} be the MLE from the Phase I sample of size \eqn{N} and
#' \eqn{T} the monitored statistic. The procedure is:
#' \enumerate{
#'   \item Draw `N_ref` subgroups of size `n` from \eqn{F(\cdot;\hat\theta)},
#'     fit and compute \eqn{T}; call this reference sample \eqn{T^{ref}}.
#'     It plays the role of the distribution of a future in-control statistic.
#'   \item For \eqn{r = 1,\dots,R}: simulate a Phase I sample of size \eqn{N} from
#'     \eqn{F(\cdot;\hat\theta)}, estimate \eqn{\tilde\theta_r}, and build the
#'     bootstrap distribution \eqn{T^*_r} (`B` values) from
#'     \eqn{F(\cdot;\tilde\theta_r)}. The limits at level \eqn{\alpha'} are the
#'     \eqn{\alpha'/2} and \eqn{1-\alpha'/2} empirical quantiles of \eqn{T^*_r}.
#'   \item The conditional false-alarm probability of replicate \eqn{r} is the
#'     fraction of \eqn{T^{ref}} outside those limits, \eqn{p_r(\alpha')}.
#'     Its mean over \eqn{r} estimates the unconditional false-alarm rate, and
#'     its upper quantile measures how bad a single unlucky Phase I sample can be.
#'   \item `adjust = "mean"` chooses the level \eqn{\alpha'} for which the mean of
#'     \eqn{p_r(\alpha')} equals `alpha`; `adjust = "quantile"` chooses the
#'     largest \eqn{\alpha'} for which \eqn{p_r(\alpha') \le} `alpha` with
#'     probability `prob`; `adjust = "none"` keeps \eqn{\alpha' =} `alpha` and only
#'     reports the diagnostics. The final limits are computed from the
#'     \eqn{\hat\theta}-bootstrap with the chosen level.
#' }
#' The cost is about `R * B + N_ref` model fits, so supply a closed-form
#' `fit_fn` for large problems.
#'
#' @param data Phase I data: a matrix or data frame with one subgroup per row
#'   (then `n` is the number of columns), or a numeric vector (then `n` must be
#'   given and the observations are pooled).
#' @param distr,u,alpha,statistic,start,fit_fn,optim_method,lower,upper,max_fit_tries,env
#'   As in [bootqcc()]. `start` is required unless `fit_fn` is given.
#' @param n Subgroup size used in Phase II. Defaults to the number of columns of
#'   `data` when it is a matrix.
#' @param B Bootstrap replicates for the limits (default
#'   `max(1000, ceiling(20 / alpha))`).
#' @param R Number of outer replicates that resample the Phase I estimation.
#' @param N_ref Size of the reference sample of in-control statistics (default
#'   `max(5000, ceiling(50 / alpha))`).
#' @param adjust One of `"none"`, `"mean"`, `"quantile"`, see Details.
#' @param prob Probability for the `"quantile"` criterion.
#' @param seed Optional seed.
#' @return An object of classes `bootqcc_phase1` and `bootqcc`, so it works with
#'   [bootqcc_monitor()], [arl()] and `plot()`. In addition to the usual
#'   fields (`limits` are the final limits) it contains `params_hat`,
#'   `limits_unadjusted`, `alpha_adj` (the level used) with `alpha_mean` and
#'   `alpha_quantile` (the levels each criterion would give), `fap_mean` and `fap_quantile` (the
#'   mean and the `prob`-quantile of the conditional false-alarm probability at
#'   the nominal level, before adjustment), `fap_adjusted` (the same at the
#'   final level) and `fap` (the conditional probabilities for each outer
#'   replicate at the nominal level).
#' @references
#' Quesenberry, C. P. (1993). The effect of sample size on estimated limits for
#' X-bar and X control charts. \emph{Journal of Quality Technology}, 25(4), 237--247.
#'
#' Jensen, W. A., Jones-Farmer, L. A., Champ, C. W. and Woodall, W. H. (2006).
#' Effects of parameter estimation on control chart properties: a literature
#' review. \emph{Journal of Quality Technology}, 38(4), 349--364.
#'
#' Gandy, A. and Kvaloy, J. T. (2013). Guaranteed conditional performance of
#' control charts via bootstrap methods. \emph{Scandinavian Journal of
#' Statistics}, 40(4), 647--668.
#'
#' Saleh, N. A., Mahmoud, M. A., Keefe, M. J. and Woodall, W. H. (2015). The
#' difficulty in designing Shewhart X-bar and X control charts with estimated
#' parameters. \emph{Computational Statistics and Data Analysis}, 87, 1--16.
#' @examples
#' \donttest{
#' set.seed(1)
#' phase1 <- matrix(rexp(20 * 5, rate = 1), ncol = 5)   # 20 subgroups of size 5
#' ch <- bootqcc_phase1(phase1, "exp", u = 0.5, alpha = 0.05,
#'                      fit_fn = function(x) list(rate = 1 / mean(x)),
#'                      B = 1000, R = 100, N_ref = 5000, adjust = "mean", seed = 2)
#' ch
#' }
#' @export
bootqcc_phase1 <- function(data, distr, n = NULL, u = 0.5, alpha = 0.0027,
                           B = max(1000, ceiling(20 / alpha)), R = 200,
                           N_ref = max(5000, ceiling(50 / alpha)),
                           adjust = c("none", "mean", "quantile"), prob = 0.9,
                           statistic = NULL, start = NULL, fit_fn = NULL,
                           optim_method = "Nelder-Mead", lower = -Inf, upper = Inf,
                           max_fit_tries = 20, seed = NULL, env = parent.frame()) {
  adjust <- match.arg(adjust)
  d <- .resolve_distr(distr, env)
  if (is.data.frame(data)) data <- as.matrix(data)
  if (is.matrix(data)) {
    if (is.null(n)) n <- ncol(data)
  } else if (is.null(n)) {
    stop("'n' (the Phase II subgroup size) is required when 'data' is a vector.",
         call. = FALSE)
  }
  x <- as.numeric(data)
  x <- x[is.finite(x)]
  N <- length(x)
  n <- as.integer(n)
  if (N < 2L * n) stop("Too few Phase I observations.", call. = FALSE)
  if (!(is.numeric(alpha) && length(alpha) == 1L && alpha > 0 && alpha < 1)) {
    stop("'alpha' must be a single number in (0, 1).", call. = FALSE)
  }
  if (is.null(fit_fn) && is.null(start)) {
    stop("Give starting values in 'start' (or a closed-form 'fit_fn').", call. = FALSE)
  }
  B <- as.integer(B); R <- as.integer(R); N_ref <- as.integer(N_ref)
  if (B * alpha / 2 < 5) {
    warning(sprintf("B * alpha / 2 = %.1f is small; use B >= %d for stable limits.",
                    B * alpha / 2, as.integer(ceiling(10 / alpha))), call. = FALSE)
  }
  if (N_ref * alpha / 2 < 5) {
    warning(sprintf("N_ref * alpha / 2 = %.1f is small; use N_ref >= %d.",
                    N_ref * alpha / 2, as.integer(ceiling(10 / alpha))), call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  stat_fn <- if (is.null(statistic)) {
    function(est) .call_fun(d$q, u, est)
  } else statistic
  fitter <- if (is.null(fit_fn)) {
    .make_mle_fitter(d$d, .check_params(start, "start"), optim_method, lower, upper)
  } else fit_fn

  theta_hat <- fitter(x)
  if (is.null(theta_hat)) stop("The model could not be fitted to the Phase I data.", call. = FALSE)

  # k statistics of subgroups of size `size` drawn from `par`
  draw_stats <- function(par, k, size = n) {
    out <- numeric(k)
    got <- 0L
    tries <- 0L
    while (got < k) {
      tries <- tries + 1L
      if (tries > max_fit_tries * k) {
        stop("Too many failed fits; check the model and the starting values.", call. = FALSE)
      }
      z <- tryCatch(suppressWarnings(.call_fun(d$r, size, par)), error = function(e) NULL)
      if (is.null(z) || length(z) != size || any(!is.finite(z))) next
      est <- fitter(z)
      if (is.null(est)) next
      s <- tryCatch(suppressWarnings(stat_fn(est)), error = function(e) NULL)
      if (is.null(s) || length(s) != 1L || !is.finite(s)) next
      got <- got + 1L
      out[got] <- s
    }
    out
  }
  qlim <- function(s, a) stats::quantile(s, c(a / 2, 1 - a / 2), names = FALSE)

  # --- unadjusted limits from the estimated parameters --------------------
  boot_stats <- draw_stats(theta_hat, B)
  lim0 <- qlim(boot_stats, alpha)

  # --- outer bootstrap: conditional false-alarm probability ----------------
  t_ref <- sort(draw_stats(theta_hat, N_ref))
  fap_at <- function(lim) {
    (findInterval(lim[1], t_ref, left.open = TRUE) +
       (N_ref - findInterval(lim[2], t_ref))) / N_ref
  }
  # Draw the outer Phase I samples; keep the sorted bootstrap statistics
  tstar <- vector("list", R)
  for (r in seq_len(R)) {
    th <- NULL
    tries <- 0L
    while (is.null(th)) {
      tries <- tries + 1L
      if (tries > max_fit_tries * R) stop("Too many failed fits in the outer bootstrap.", call. = FALSE)
      xr <- tryCatch(suppressWarnings(.call_fun(d$r, N, theta_hat)), error = function(e) NULL)
      if (!is.null(xr) && length(xr) == N && all(is.finite(xr))) th <- fitter(xr)
    }
    tstar[[r]] <- sort(draw_stats(th, B))
  }
  fap_for <- function(a) vapply(tstar, function(s) fap_at(qlim(s, a)), numeric(1))
  fap0 <- fap_for(alpha)

  # criteria increase with the level, so bisect on the log scale
  solve_level <- function(crit) {
    if (crit(alpha) <= alpha) return(alpha)   # already conservative enough
    lo <- log(alpha / 100); hi <- log(alpha)
    for (i in 1:25) {
      mid <- (lo + hi) / 2
      if (crit(exp(mid)) > alpha) hi <- mid else lo <- mid
    }
    exp(lo)
  }
  alpha_mean <- solve_level(function(a) mean(fap_for(a)))
  alpha_quantile <- solve_level(function(a) stats::quantile(fap_for(a), prob, names = FALSE))
  alpha_adj <- switch(adjust, none = alpha, mean = alpha_mean, quantile = alpha_quantile)
  limits <- if (alpha_adj == alpha) lim0 else qlim(boot_stats, alpha_adj)
  fap1 <- if (alpha_adj == alpha) fap0 else fap_for(alpha_adj)

  structure(
    list(call = match.call(), distr = d, distr_name = d$name, params = theta_hat,
         params_hat = theta_hat, u = if (is.null(statistic)) u else NA_real_,
         n = n, N = N, alpha = alpha, alpha_adj = alpha_adj,
         alpha_mean = alpha_mean, alpha_quantile = alpha_quantile, adjust = adjust, prob = prob,
         B = B, M = 1L, R = R, N_ref = N_ref,
         limits = c(LCL = limits[1], UCL = limits[2]),
         limits_unadjusted = c(LCL = lim0[1], UCL = lim0[2]),
         sd_limits = c(LCL = NA_real_, UCL = NA_real_),
         lcl_runs = limits[1], ucl_runs = limits[2], boot_stats = boot_stats,
         fap = fap0, fap_mean = mean(fap0),
         fap_quantile = stats::quantile(fap0, prob, names = FALSE),
         fap_adjusted = fap1,
         statistic = stat_fn, fitter = fitter, max_run_length = 1e5, arl = NULL),
    class = c("bootqcc_phase1", "bootqcc"))
}

#' @export
print.bootqcc_phase1 <- function(x, digits = 4, ...) {
  what <- if (is.na(x$u)) "custom statistic" else sprintf("%g-quantile", x$u)
  cat("Parametric bootstrap control chart from Phase I data\n")
  cat(sprintf("Distribution : %s, estimated from N = %d observations\n", x$distr_name, x$N))
  cat("Estimates    : ", paste(names(x$params_hat), signif(unlist(x$params_hat), digits),
                              sep = " = ", collapse = ", "), "\n", sep = "")
  cat(sprintf("Statistic    : %s, subgroup size n = %d\n", what, x$n))
  cat(sprintf("Nominal alpha = %g; adjustment = %s; level used = %s\n",
              x$alpha, x$adjust, format(x$alpha_adj, digits = 3)))
  cat(sprintf("LCL = %s, UCL = %s", format(x$limits[["LCL"]], digits = digits),
              format(x$limits[["UCL"]], digits = digits)))
  if (x$adjust != "none") {
    cat(sprintf("  (unadjusted: %s, %s)", format(x$limits_unadjusted[["LCL"]], digits = digits),
                format(x$limits_unadjusted[["UCL"]], digits = digits)))
  }
  cat("\n")
  cat(sprintf("Conditional false-alarm probability at the nominal level: mean %s, %g%% quantile %s\n",
              format(x$fap_mean, digits = 3), 100 * x$prob, format(x$fap_quantile, digits = 3)))
  cat(sprintf("                       at the level used:               mean %s, %g%% quantile %s\n",
              format(mean(x$fap_adjusted), digits = 3), 100 * x$prob,
              format(stats::quantile(x$fap_adjusted, x$prob, names = FALSE), digits = 3)))
  invisible(x)
}
