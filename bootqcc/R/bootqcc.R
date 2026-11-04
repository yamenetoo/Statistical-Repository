#' Parametric bootstrap control chart for any distribution
#'
#' Computes Phase I control limits for a quantile (or any user-defined
#' statistic) of a parametric distribution, then estimates the in-control
#' average run length (ARL) by Monte Carlo simulation.
#'
#' The procedure, repeated `M` times, is:
#' 1. draw `B` bootstrap samples of size `n` from the distribution with the
#'    given `params`;
#' 2. fit the distribution to each sample by maximum likelihood;
#' 3. evaluate the monitored statistic (the `u`-quantile by default) at each
#'    fit;
#' 4. take the `alpha/2` and `1 - alpha/2` empirical quantiles of the `B`
#'    statistics as one lower and upper limit.
#'
#' The reported limits are the means of the `M` limits. The ARL is then
#' estimated by drawing samples until the statistic leaves the limits.
#'
#' @param distr A distribution name such as `"weibull"` or `"gamma"` (the
#'   functions `d<name>`, `q<name>` and `r<name>` must exist), or a named list
#'   `list(d = , q = , r = , name = )` of user-defined functions. Each
#'   function takes the data (or `n`, or probabilities) as its first argument
#'   followed by the distribution parameters by name.
#' @param params Named list of the true (in-control) parameters, for example
#'   `list(shape = 2, scale = 1)`.
#' @param u Probability level of the monitored quantile, in (0, 1). Ignored if
#'   `statistic` is given.
#' @param n Subgroup (sample) size.
#' @param alpha Total false-alarm probability used for the limits.
#' @param B Number of bootstrap replicates per Monte Carlo iteration.
#' @param M Number of Monte Carlo iterations of the limit computation.
#' @param arl_rep Number of simulated run lengths used to estimate the ARL.
#' @param statistic Optional function taking the named list of fitted
#'   parameters and returning one number, to monitor a statistic other than a
#'   quantile (a mean, a reliability, and so on).
#' @param start Named list of starting values for the maximum likelihood
#'   optimisation. Defaults to `params`.
#' @param fit_fn Optional custom estimator: a function of a sample returning a
#'   named list of parameter estimates, or `NULL` on failure. Replaces the
#'   built-in maximum likelihood fit.
#' @param optim_method Method passed to [stats::optim()] for the built-in fit.
#' @param lower,upper Parameter bounds, used with `"L-BFGS-B"` or `"Brent"`.
#' @param max_run_length Run lengths are capped at this value; capped runs are
#'   counted in the `censored` element.
#' @param max_fit_tries Multiplier bounding failed fits: at most
#'   `max_fit_tries * B` attempts are made per Monte Carlo iteration.
#' @param seed Optional random seed.
#' @param verbose Print progress messages.
#' @param env Environment used to look up distribution functions by name.
#'
#' @return An object of class `"bootqcc"`, a list with elements `limits`
#'   (named vector `LCL`, `UCL`), `sd_limits`, `arl` (a list with `mean`, `sd`,
#'   `runs`, `censored`), the per-iteration limits `lcl_runs` and `ucl_runs`,
#'   the bootstrap statistics of the last iteration `boot_stats`, and the
#'   settings used.
#'
#' @examples
#' # Small settings so that the example runs quickly; use the defaults
#' # (B = 1000, M = 30, arl_rep = 100) for a real analysis.
#' ch <- bootqcc("weibull", list(shape = 2, scale = 1), u = 0.1, n = 5,
#'               alpha = 0.05, B = 100, M = 3, arl_rep = 10, seed = 1)
#' ch
#'
#' \donttest{
#' # A user-defined distribution (Kumaraswamy)
#' kw <- list(
#'   name = "kumaraswamy",
#'   d = function(x, a, b) a * b * x^(a - 1) * (1 - x^a)^(b - 1),
#'   q = function(p, a, b) (1 - (1 - p)^(1 / b))^(1 / a),
#'   r = function(n, a, b) (1 - (1 - runif(n))^(1 / b))^(1 / a)
#' )
#' bootqcc(kw, list(a = 2, b = 3), u = 0.5, n = 10,
#'         alpha = 0.05, B = 100, M = 3, arl_rep = 10, seed = 1)
#' }
#' @export
bootqcc <- function(distr, params, u = 0.5, n = 5, alpha = 0.0027,
                    B = 1000, M = 30, arl_rep = 100,
                    statistic = NULL, start = params, fit_fn = NULL,
                    optim_method = "Nelder-Mead", lower = -Inf, upper = Inf,
                    max_run_length = 1e5, max_fit_tries = 20,
                    seed = NULL, verbose = FALSE, env = parent.frame()) {
  d <- .resolve_distr(distr, env)
  params <- .check_params(params)
  start <- .check_params(start, "start")
  if (!is.null(statistic) && !is.function(statistic)) {
    stop("'statistic' must be a function of the fitted parameter list.", call. = FALSE)
  }
  if (is.null(statistic) && !(is.numeric(u) && length(u) == 1L && u > 0 && u < 1)) {
    stop("'u' must be a single number in (0, 1).", call. = FALSE)
  }
  if (!(is.numeric(alpha) && length(alpha) == 1L && alpha > 0 && alpha < 1)) {
    stop("'alpha' must be a single number in (0, 1).", call. = FALSE)
  }
  n <- as.integer(n); B <- as.integer(B); M <- as.integer(M)
  arl_rep <- as.integer(arl_rep)
  if (any(is.na(c(n, B, M, arl_rep))) || any(c(n, B, M, arl_rep) < 1L)) {
    stop("'n', 'B', 'M' and 'arl_rep' must be positive integers.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  stat_fn <- if (is.null(statistic)) {
    function(est) .call_fun(d$q, u, est)
  } else statistic
  fitter <- if (is.null(fit_fn)) {
    .make_mle_fitter(d$d, start, optim_method, lower, upper)
  } else fit_fn

  # Statistic of one sample, or NULL if the fit or the statistic fails
  sample_stat <- function(x) {
    est <- fitter(x)
    if (is.null(est)) return(NULL)
    s <- tryCatch(suppressWarnings(stat_fn(est)), error = function(e) NULL)
    if (is.null(s) || length(s) != 1L || !is.finite(s)) return(NULL)
    s
  }
  draw <- function() {
    tryCatch(suppressWarnings(.call_fun(d$r, n, params)), error = function(e) NULL)
  }
  one_stat <- function() {
    x <- draw()
    if (is.null(x) || length(x) != n || any(!is.finite(x))) return(NULL)
    sample_stat(x)
  }

  # ---- Phase I: bootstrap limits ----------------------------------------
  lcl <- ucl <- numeric(M)
  boot_stats <- NULL
  for (m in seq_len(M)) {
    stats_b <- numeric(B)
    got <- 0L
    tries <- 0L
    while (got < B) {
      tries <- tries + 1L
      if (tries > max_fit_tries * B) {
        stop("Too many failed fits while building the limits. Check the ",
             "parameters, the starting values or the distribution functions.",
             call. = FALSE)
      }
      s <- one_stat()
      if (!is.null(s)) {
        got <- got + 1L
        stats_b[got] <- s
      }
    }
    qs <- stats::quantile(stats_b, c(alpha / 2, 1 - alpha / 2), names = FALSE)
    lcl[m] <- qs[1]
    ucl[m] <- qs[2]
    boot_stats <- stats_b
    if (verbose) message(sprintf("Iteration %d/%d: LCL = %.6g, UCL = %.6g",
                                 m, M, lcl[m], ucl[m]))
  }
  limits <- c(LCL = mean(lcl), UCL = mean(ucl))
  sd_limits <- c(LCL = if (M > 1) stats::sd(lcl) else NA_real_,
                 UCL = if (M > 1) stats::sd(ucl) else NA_real_)

  obj <- structure(
    list(call = match.call(), distr = d, distr_name = d$name, params = params,
         u = if (is.null(statistic)) u else NA_real_, n = n, alpha = alpha,
         B = B, M = M, limits = limits, sd_limits = sd_limits,
         lcl_runs = lcl, ucl_runs = ucl, boot_stats = boot_stats,
         statistic = stat_fn, fitter = fitter,
         max_run_length = max_run_length, arl = NULL),
    class = "bootqcc")

  # ---- Phase II: in-control ARL -----------------------------------------
  obj$arl <- arl(obj, params = params, rep = arl_rep, verbose = verbose)
  obj
}

#' Estimate the average run length of a fitted chart
#'
#' Simulates run lengths until the monitored statistic leaves the control
#' limits. Use the in-control parameters to estimate ARL0, or shifted
#' parameters to estimate the out-of-control ARL1.
#'
#' @param chart An object returned by [bootqcc()].
#' @param params Named list of the parameters to simulate from. Defaults to
#'   the in-control parameters of the chart.
#' @param rep Number of simulated run lengths.
#' @param max_run_length Cap on a single run length. Defaults to the value
#'   stored in the chart.
#' @param verbose Print a short summary.
#' @return A list with `mean`, `sd`, `runs` and `censored` (the number of
#'   runs that reached the cap).
#' @examples
#' ch <- bootqcc("weibull", list(shape = 2, scale = 1), u = 0.1, n = 5,
#'               alpha = 0.05, B = 100, M = 3, arl_rep = 10, seed = 1)
#' # Out-of-control ARL after a scale shift
#' arl(ch, params = list(shape = 2, scale = 1.5), rep = 20)$mean
#' @export
arl <- function(chart, params = chart$params, rep = 100,
                max_run_length = chart$max_run_length, verbose = FALSE) {
  stopifnot(inherits(chart, "bootqcc"))
  params <- .check_params(params)
  lim <- chart$limits
  n <- chart$n
  runs <- numeric(rep)
  for (i in seq_len(rep)) {
    len <- 0
    repeat {
      x <- tryCatch(suppressWarnings(.call_fun(chart$distr$r, n, params)),
                    error = function(e) NULL)
      if (!is.null(x) && length(x) == n && all(is.finite(x))) {
        est <- chart$fitter(x)
        if (!is.null(est)) {
          s <- tryCatch(suppressWarnings(chart$statistic(est)),
                        error = function(e) NULL)
          if (!is.null(s) && length(s) == 1L && is.finite(s)) {
            len <- len + 1
            if (s < lim[["LCL"]] || s > lim[["UCL"]]) break
          }
        }
      }
      if (len >= max_run_length) break
    }
    runs[i] <- len
  }
  out <- list(mean = mean(runs),
              sd = if (rep > 1) stats::sd(runs) else NA_real_,
              runs = runs, censored = sum(runs >= max_run_length))
  if (verbose) {
    message(sprintf("ARL = %.2f (SD %.2f), %d run(s) censored",
                    out$mean, out$sd, out$censored))
  }
  out
}
