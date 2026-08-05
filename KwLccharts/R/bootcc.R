# =========================================================================
# Maximum likelihood estimation and parametric bootstrap percentile
# control chart for the Kumaraswamy-Lindley (KwL) distribution.
#
# Requires: lambertW.R, KwL_distribution.R  (source before this file)
# =========================================================================

# ---- Negative log-likelihood ---------------------------------------------
#' Negative log-likelihood of the KwL distribution
#'
#' @param par Numeric vector \code{c(theta, a, b)}.
#' @param y Numeric vector of observations.
#' @return Negative log-likelihood (a large sentinel value if \code{par} is
#'   infeasible or the likelihood is non-finite).
#' @export
negloglik_KwL <- function(par, y) {
  theta <- par[1]; a <- par[2]; b <- par[3]
  if (theta <= 0 || a <= 0 || b <= 0) return(1e10)
  ll <- sum(dKwL(y, theta, a, b, log = TRUE))
  if (!is.finite(ll)) return(1e10)
  -ll
}

# ---- MLE fit ---------------------------------------------------------------
#' Maximum-likelihood fit of the KwL distribution
#'
#' Maximizes the KwL log-likelihood via \code{\link[stats]{optim}}
#' (\code{L-BFGS-B}), restarting from five dispersed starting values and
#' keeping the best (highest-likelihood) converged solution, to reduce
#' sensitivity to local optima -- important for small subgroup sizes.
#'
#' @param y Numeric vector of observations (\eqn{y_i>0}).
#' @param start Named numeric vector, initial values for \code{c(theta,a,b)}
#'   used as the first starting point.
#' @return A list with elements \code{theta}, \code{a}, \code{b} (the MLEs)
#'   and \code{value} (the minimized negative log-likelihood), or \code{NULL}
#'   if no starting value converged.
#' @examples
#' set.seed(1)
#' y <- rKwL(200, theta = 1, a = 1.5, b = 2)
#' fit_KwL(y)
#' @export
fit_KwL <- function(y, start = c(theta = 1, a = 1, b = 1)) {
  starts <- list(
    start,
    c(theta = 0.5, a = 1, b = 1),
    c(theta = 2,   a = 1, b = 1),
    c(theta = 1,   a = 2, b = 0.5),
    c(theta = 1,   a = 0.5, b = 2)
  )
  best <- NULL
  for (s in starts) {
    opt <- try(
      stats::optim(
        par = s, fn = negloglik_KwL, y = y,
        method = "L-BFGS-B",
        lower = c(1e-3, 1e-3, 1e-3), upper = c(50, 50, 50)
      ),
      silent = TRUE
    )
    if (!inherits(opt, "try-error") && opt$convergence == 0 && is.finite(opt$value)) {
      if (is.null(best) || opt$value < best$value) best <- opt
    }
  }
  if (is.null(best)) return(NULL)
  list(theta = best$par[1], a = best$par[2], b = best$par[3], value = best$value)
}

# ---- Algorithm 1 analogue: generate a bootstrap subgroup -------------------
# (rKwL already implements this via inversion)

# ---- Algorithm 2 analogue: bootstrap percentile control chart -------------
#
#' Bootstrap percentile control chart for the KwL distribution
#'
#' Implements the parametric-bootstrap percentile control-chart algorithm of
#' Al-Lami, Jabbari Khamnei & Heydari (2025, PLOS ONE, for the
#' Lindley-Geometric distribution), adapted to the three-parameter
#' Kumaraswamy-Lindley distribution: (1) draw a bootstrap subgroup of size
#' \code{n} from KwL(theta,a,b); (2) refit the KwL parameters by MLE; (3)
#' compute the bootstrap percentile estimate \eqn{Q^*(u)}; repeat \code{B}
#' times and take the \eqn{\alpha/2} and \eqn{1-\alpha/2} order statistics as
#' the control limits. The whole procedure is repeated \code{MC} times to
#' obtain the mean and Monte-Carlo standard deviation of the limits.
#'
#' @param u Percentile of interest, \eqn{0<u<1}.
#' @param alpha Nominal false-alarm probability (e.g. 0.0027).
#' @param theta,a,b KwL parameters used to generate the Phase-I / bootstrap
#'   subgroups (typically the Phase-I MLEs).
#' @param n Subgroup size.
#' @param B Number of bootstrap iterations per Monte-Carlo repetition.
#' @param MC Number of outer Monte-Carlo repetitions.
#' @param start Starting values passed to \code{\link{fit_KwL}}.
#' @param verbose If \code{TRUE}, prints progress for each Monte-Carlo
#'   repetition.
#' @return A list with \code{MLCL}, \code{MUCL} (mean control limits),
#'   \code{SDLCL}, \code{SDUCL} (their Monte-Carlo standard deviations), and
#'   \code{LCL_all}, \code{UCL_all} (the per-repetition limits).
#' @examples
#' \donttest{
#' set.seed(1)
#' KwLBootCC(u = 0.05, alpha = 0.0027, theta = 1, a = 1.5, b = 1.5,
#'           n = 5, B = 200, MC = 5)
#' }
#' @export
KwLBootCC <- function(u, alpha, theta, a, b, n = 5, B = 2000, MC = 50,
                       start = c(theta = 1, a = 1, b = 1), verbose = FALSE) {
  LCL <- numeric(MC); UCL <- numeric(MC)

  for (m in seq_len(MC)) {
    QS <- rep(NA_real_, B)
    for (i in seq_len(B)) {
      samp <- rKwL(n, theta, a, b)
      fit  <- fit_KwL(samp, start = start)
      if (!is.null(fit)) {
        QS[i] <- qKwL(u, fit$theta, fit$a, fit$b)
      }
    }
    QS <- sort(QS[!is.na(QS) & is.finite(QS)])
    nb <- length(QS)
    lo <- max(1, floor(nb * alpha / 2))
    hi <- min(nb, ceiling(nb * (1 - alpha / 2)))
    LCL[m] <- QS[lo]
    UCL[m] <- QS[hi]
    if (verbose) cat(sprintf("MC %d/%d: LCL=%.5f UCL=%.5f (n_valid=%d)\n",
                              m, MC, LCL[m], UCL[m], nb))
  }

  list(
    MLCL = mean(LCL), MUCL = mean(UCL),
    SDLCL = stats::sd(LCL),  SDUCL = stats::sd(UCL),
    LCL_all = LCL, UCL_all = UCL
  )
}

# ---- ARL simulation (in-control or shifted parameters) --------------------
#
#' Average run length simulation for the KwL bootstrap control chart
#'
#' Simulates Phase-II monitoring: repeatedly draws subgroups of size
#' \code{n} from KwL(theta1,a1,b1), re-fits the parameters, computes
#' \eqn{\hat Q(u)}, and counts the number of in-control subgroups until
#' \eqn{\hat Q(u)} falls outside \code{[LCL, UCL]} (the run length). Use
#' \code{theta1 = theta}, \code{a1 = a}, \code{b1 = b} (the Phase-I
#' parameters) for the in-control ARL, or different values for the
#' out-of-control ARL under a specified parameter shift.
#'
#' @param u Percentile of interest.
#' @param LCL,UCL Control limits (typically from \code{\link{KwLBootCC}}).
#' @param theta1,a1,b1 KwL parameters generating the Phase-II subgroups.
#' @param n Subgroup size.
#' @param n_runs Number of independent run lengths to simulate.
#' @param max_run Cap on the run length (to guarantee termination).
#' @param start Starting values passed to \code{\link{fit_KwL}}.
#' @return A list with \code{ARL} (mean run length), \code{SDRL} (its
#'   standard deviation), and \code{RL_all} (the individual run lengths).
#' @examples
#' \donttest{
#' set.seed(1)
#' lim <- KwLBootCC(0.05, 0.0027, 1, 1.5, 1.5, n = 5, B = 200, MC = 5)
#' KwL_ARL(0.05, lim$MLCL, lim$MUCL, 1, 1.5, 1.5, n_runs = 10, max_run = 300)
#' }
#' @export
KwL_ARL <- function(u, LCL, UCL, theta1, a1, b1, n = 5, n_runs = 200,
                     max_run = 2000, start = c(theta = 1, a = 1, b = 1)) {
  RL <- numeric(n_runs)
  for (r in seq_len(n_runs)) {
    count <- 0
    inctrl <- TRUE
    while (inctrl && count < max_run) {
      count <- count + 1
      samp <- rKwL(n, theta1, a1, b1)
      fit  <- fit_KwL(samp, start = start)
      if (is.null(fit)) next
      Qhat <- qKwL(u, fit$theta, fit$a, fit$b)
      if (!is.finite(Qhat) || Qhat < LCL || Qhat > UCL) inctrl <- FALSE
    }
    RL[r] <- count
  }
  list(ARL = mean(RL), SDRL = stats::sd(RL), RL_all = RL)
}
