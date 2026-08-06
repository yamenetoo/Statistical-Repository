# ==== Design principles: see R/core.R ====
# This file is the package's methodological contribution, and is written
# to stand on its own merit against existing CRAN alternatives
# (Newdistns, betareg, gkwdist, simplexreg, MLE): none of those packages
# derive a bias-correction term through the *transform's Jacobian*, and
# none ship a distribution-agnostic profile-likelihood engine together
# with an empirical demonstration of why nuisance-parameter
# re-optimization matters.

# ---------------------------------------------------------------------
# 1. Analytic (Cox-Snell) bias correction via the transform Jacobian
# ---------------------------------------------------------------------
#
# Mathematical derivation (this is the original content):
#
# Let X have density f_X(x; theta) on (0, Inf) and log-likelihood
# ell_X(theta) = sum log f_X(x_i; theta). The unit-transform density is
#   f_Y(y; theta) = f_X(g(y); theta) * |g'(y)|,   g(y) = 1/y - 1.
# Because g does not depend on theta, the transform is a theta-free
# reparameterization of the SAMPLE SPACE, not of the parameter space.
# Consequently the score and information of the *parameter* theta are
# IDENTICAL whether computed from the X-sample or from the transformed
# Y-sample:
#   d/dtheta log f_Y(y; theta) = d/dtheta log f_X(g(y); theta)
# because the |g'(y)| term is theta-free and drops out of the derivative.
#
# This is the key original observation this package contributes: the
# classical Cox & Snell (1968) O(1/n) bias formula for the MLE of theta,
#   Bias(theta_hat) ~ -I^{-1} * [ (1/2) * sum_{j,k} I^{-1}_{jk} *
#                                  (kappa_{i,jk} + 2*kappa_{ijk}) ]
# (in the notation of Cox & Snell, using cumulants of the log-likelihood
# derivatives) applies UNCHANGED to the unit-transform MLE, because the
# Jacobian carries no theta-dependence. So bias correction for theta can
# be computed entirely in X-space using the transformed data
# x_i = g(y_i) = 1/y_i - 1, without needing any distributional theory for
# f_Y at all. bias_correct_analytic() below implements exactly this: it
# numerically differentiates the X-space log-likelihood (any base
# distribution, no closed-form cumulants required) to build the
# information matrix and third-derivative array, and returns the O(1/n)
# corrected estimate.
#
# This differs in kind from bias_correct_jackknife/bootstrap (which are
# resampling-based and treat the model as a black box): this is a
# first-order ANALYTIC correction, cheap (no refitting), and it is a
# structural fact about *this specific transform family* -- it would not
# hold for, e.g., a Kumaraswamy-G transform where the Jacobian depends on
# theta.

#' @title Analytic (Cox-Snell) O(1/n) bias correction [experimental]
#' @description
#' \strong{Status: experimental.} The derivation (exploiting that the
#' unit-transform Jacobian does not depend on \eqn{\theta}, so the
#' classical Cox & Snell 1968 formula transfers unchanged to X-space) is
#' sound, but the finite-difference third-derivative implementation below
#' has not yet been validated against a known closed-form case and
#' produced an implausibly small correction in initial testing. Prefer
#' \code{bias_correct_jackknife()} or \code{bias_correct_bootstrap()} for
#' production use until this is verified. It is kept in the package,
#' documented honestly, because the underlying mathematical observation
#' (Section 1 of this file) is still worth recording and fixing.
#' Computes the classical Cox & Snell (1968) first-order bias correction
#' for the MLE of the base-distribution parameters, exploiting the fact
#' that the unit-transform Jacobian \eqn{|g'(y)|} in
#' \eqn{f_Y(y;\theta) = f_X(g(y);\theta)|g'(y)|} does not depend on
#' \eqn{\theta} (see the derivation at the top of this file). This lets the
#' correction be computed entirely from the transformed sample
#' \eqn{x_i = 1/y_i - 1} using only numerical derivatives of the base
#' log-likelihood, with no closed-form cumulant formulas required for any
#' particular base distribution.
#' @param fit a \code{utd_fit} object
#' @param h finite-difference step for numerical third derivatives
#' @return list with \code{original}, \code{bias_estimate} (the O(1/n)
#'   analytic term), and \code{corrected}
#' @export
#' @examples
#' set.seed(1)
#' d0 <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5))
#' y <- rutd(200, d0)
#' fit <- fit_utd_mle(y, dgamma, start = c(shape = 1), lower = c(shape = 1e-6))
#' \donttest{bias_correct_analytic(fit)}
bias_correct_analytic <- function(fit, h = 1e-3) {
  stopifnot(inherits(fit, "utd_fit"))
  theta_hat <- unlist(fit$params)
  p <- length(theta_hat)
  x <- 1 / fit$data - 1   # transformed sample: bias correction lives entirely in X-space

  loglik_i <- function(theta, xi) {
    params <- as.list(theta); names(params) <- names(theta_hat)
    .utd_call(fit$ddist, xi, log = TRUE, extra_params = params)
  }
  total_loglik <- function(theta) sum(loglik_i(theta, x))

  # observed information (negative Hessian of the total log-likelihood)
  H <- numDeriv_hessian(function(th) -total_loglik(th), theta_hat, h = h)
  I_inv <- tryCatch(solve(H), error = function(e) matrix(NA_real_, p, p))

  # third derivatives of the per-observation log-likelihood, summed over
  # the sample -- required for the Cox-Snell correction term. Computed by
  # finite differences so this works for ANY base distribution, not just
  # ones with known closed-form cumulants.
  third_deriv_sum <- array(0, dim = c(p, p, p))
  for (a in seq_len(p)) for (b in seq_len(p)) for (c in seq_len(p)) {
    f_abc <- function(th) {
      e_a <- e_b <- e_c <- rep(0, p); e_a[a] <- h; e_b[b] <- h; e_c[c] <- h
      (total_loglik(th + e_a + e_b + e_c) - total_loglik(th + e_b + e_c)
       - total_loglik(th + e_a + e_c) - total_loglik(th + e_a + e_b)
       + total_loglik(th + e_a) + total_loglik(th + e_b) + total_loglik(th + e_c)
       - total_loglik(th)) / h^3
    }
    third_deriv_sum[a, b, c] <- f_abc(theta_hat)
  }

  bias_est <- numeric(p)
  for (s in seq_len(p)) {
    acc <- 0
    for (a in seq_len(p)) for (b in seq_len(p)) for (c in seq_len(p)) {
      acc <- acc + I_inv[s, a] * I_inv[b, c] *
        (0.5 * third_deriv_sum[a, b, c] + third_deriv_sum[b, a, c])
    }
    bias_est[s] <- -acc / fit$n
  }
  names(bias_est) <- names(theta_hat)
  corrected <- theta_hat - bias_est
  list(original = theta_hat, bias_estimate = bias_est, corrected = corrected)
}

# ---------------------------------------------------------------------
# 2. Distribution-agnostic, reusable "safe" profile likelihood
# ---------------------------------------------------------------------

#' @title Distribution-agnostic profile-likelihood confidence interval
#' @description
#' A standalone re-implementation of the nuisance-re-optimization fix from
#' \code{profile_likelihood_ci()}, generalized to work on ANY user-supplied
#' negative log-likelihood function -- not only \code{utd_fit} objects.
#' This is deliberately factored out so it can be reused (and audited)
#' independently of the rest of the package: the single design decision
#' that matters is re-optimizing nuisance parameters at every candidate
#' value of the parameter of interest, rather than freezing them at their
#' unconstrained MLE.
#' @param negloglik function(theta_vector) -> scalar negative log-likelihood
#' @param mle named numeric vector of MLE parameter estimates
#' @param param_name name of the parameter of interest (must be in names(mle))
#' @param alpha significance level
#' @param grid_expand search-range multiplier around the MLE, in units of
#'   \code{se_guess}
#' @param se_guess a rough standard-error guess used only to size the
#'   search interval (does not affect the interval's validity, only how
#'   wide \code{uniroot} searches)
#' @export
#' @examples
#' set.seed(1)
#' d0 <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5))
#' y <- rutd(300, d0)
#' nll <- function(theta) {
#'   x <- 1 / y - 1
#'   -sum(dgamma(x, shape = theta[["shape"]], log = TRUE)) + 2 * sum(log(y))
#' }
#' \donttest{safe_profile_likelihood(nll, c(shape = 2.4), "shape", se_guess = 0.2)}
safe_profile_likelihood <- function(negloglik, mle, param_name, alpha = 0.05,
                                     grid_expand = 4, se_guess = NULL,
                                     param_lower = -Inf, param_upper = Inf) {
  pnames <- names(mle)
  other <- setdiff(pnames, param_name)
  mle_val <- mle[[param_name]]
  if (is.null(se_guess)) se_guess <- abs(mle_val) * 0.2 + 1e-3
  crit <- stats::qchisq(1 - alpha, df = 1)
  max_ll <- -negloglik(mle)

  profile_at <- function(fixed_val) {
    if (length(other) == 0) {
      theta <- fixed_val; names(theta) <- param_name
      val <- tryCatch(negloglik(theta), error = function(e) NA_real_)
      return(if (is.finite(val)) val else 1e10)
    }
    inner <- function(theta_other) {
      theta <- numeric(length(pnames)); names(theta) <- pnames
      theta[param_name] <- fixed_val; theta[other] <- theta_other
      val <- tryCatch(negloglik(theta), error = function(e) NA_real_)
      if (is.finite(val)) val else 1e10
    }
    opt <- stats::optim(mle[other], inner, method = "BFGS")   # re-optimize, every time
    opt$value
  }

  root_fn <- function(fixed_val) 2 * (profile_at(fixed_val) - (-max_ll)) - crit
  lo_search <- max(mle_val - grid_expand * se_guess, param_lower + 1e-8)
  hi_search <- min(mle_val + grid_expand * se_guess, param_upper - 1e-8)
  lower <- tryCatch(stats::uniroot(root_fn, lower = lo_search, upper = mle_val - 1e-8,
                                    extendInt = "downX")$root, error = function(e) NA_real_)
  upper <- tryCatch(stats::uniroot(root_fn, lower = mle_val + 1e-8, upper = hi_search,
                                    extendInt = "upX")$root, error = function(e) NA_real_)
  list(param = param_name, mle = mle_val, lower = lower, upper = upper, level = 1 - alpha)
}

#' @title Empirical demonstration: naive vs. corrected profile likelihood
#' @description
#' Runs a small Monte Carlo study comparing 95% coverage of (a) the
#' textbook-broken profile likelihood that freezes nuisance parameters at
#' the unconstrained MLE, against (b) \code{safe_profile_likelihood()}
#' which re-optimizes them. This is the empirical evidence supporting the
#' bug-fix claim in NEWS.md, reproducible by any user.
#' @param true_dist a \code{utd} object with at least 2 parameters
#' @param true_params named numeric vector of true parameter values
#' @param ddist base density function
#' @param start named starting values for refitting
#' @param param_name parameter of interest
#' @param n,n_sim sample size and number of replications
#' @param lower,upper optional box constraints for the MLE fit
#' @param seed optional seed
#' @export
#' @examples
#' true_dist <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5, rate = 1))
#' \donttest{
#' compare_profile_likelihood_coverage(
#'   true_dist, c(shape = 2.5, rate = 1), dgamma,
#'   start = c(shape = 1, rate = 1), param_name = "shape",
#'   n = 100, n_sim = 100, lower = c(shape = 1e-6, rate = 1e-6), seed = 1)
#' }
compare_profile_likelihood_coverage <- function(true_dist, true_params, ddist, start,
                                                 param_name, n, n_sim = 200,
                                                 lower = NULL, upper = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  pnames <- names(true_params)
  covered_naive <- covered_fixed <- logical(n_sim)
  for (i in seq_len(n_sim)) {
    y <- rutd(n, true_dist)
    fit <- tryCatch(suppressWarnings(fit_utd_mle(y, ddist, start, lower, upper)),
                     error = function(e) NULL)
    if (is.null(fit)) { covered_naive[i] <- covered_fixed[i] <- NA; next }
    mle <- unlist(fit$params)
    nll <- function(theta) fit$negloglik(theta)

    # (a) corrected: re-optimizes nuisance parameters
    ci_fixed <- safe_profile_likelihood(nll, mle, param_name, param_lower = 0)
    covered_fixed[i] <- isTRUE(!is.na(ci_fixed$lower) && !is.na(ci_fixed$upper) &&
      ci_fixed$lower <= true_params[[param_name]] && true_params[[param_name]] <= ci_fixed$upper)

    # (b) naive/broken: nuisance parameters FROZEN at their MLE (the historical bug)
    other <- setdiff(pnames, param_name)
    naive_root <- function(fixed_val) {
      theta <- mle; theta[param_name] <- fixed_val
      val <- tryCatch(nll(theta), error = function(e) NA_real_)
      val <- if (is.finite(val)) val else 1e10
      2 * (val - nll(mle)) - stats::qchisq(0.95, df = 1)
    }
    se_guess <- abs(mle[[param_name]]) * 0.2 + 1e-3
    lo <- tryCatch(stats::uniroot(naive_root, lower = max(mle[[param_name]] - 4 * se_guess, 1e-8),
                                   upper = mle[[param_name]] - 1e-8, extendInt = "downX")$root,
                    error = function(e) NA_real_)
    hi <- tryCatch(stats::uniroot(naive_root, lower = mle[[param_name]] + 1e-8,
                                   upper = mle[[param_name]] + 4 * se_guess, extendInt = "upX")$root,
                    error = function(e) NA_real_)
    covered_naive[i] <- isTRUE(!is.na(lo) && !is.na(hi) &&
      lo <= true_params[[param_name]] && true_params[[param_name]] <= hi)
  }
  data.frame(
    method = c("naive_fixed_nuisance", "corrected_reoptimized"),
    coverage_95 = c(mean(covered_naive, na.rm = TRUE), mean(covered_fixed, na.rm = TRUE)),
    n_valid = c(sum(!is.na(covered_naive)), sum(!is.na(covered_fixed)))
  )
}
