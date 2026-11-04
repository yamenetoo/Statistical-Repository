# Validation of bootqcc against exact results for the exponential distribution
# --------------------------------------------------------------------------
# For X ~ Exp(rate = 1) the MLE of the rate from a sample of size n is 1/mean(x),
# so T = Q(u; rate_hat) = -log(1 - u) * mean(x), and mean(x) ~ Gamma(n, scale = 1/n).
# Everything below therefore has a closed form to compare against.
#
# Usage: Rscript inst/validation/exponential-validation.R   (about 15 minutes)
library(bootqcc)
cores <- max(1, parallel::detectCores() - 1)
u <- 0.5; n <- 5; cc <- -log(1 - u)
fit <- function(x) list(rate = 1 / mean(x))

## Part A: known parameters, ARL0 against 1 / alpha ------------------------
exact_limits <- function(alpha) cc * qgamma(c(alpha / 2, 1 - alpha / 2), n, scale = 1 / n)
partA <- do.call(rbind, lapply(c(0.1, 0.01, 0.0027), function(alpha) {
  do.call(rbind, lapply(c(100, 1000, 10000, 100000), function(B) {
    ch <- suppressWarnings(bootqcc("exp", list(rate = 1), u = u, n = n, alpha = alpha, B = B,
                                   M = 10, arl_rep = 2000, fit_fn = fit, seed = 11))
    data.frame(alpha = alpha, B = B, B_alpha_2 = B * alpha / 2,
               LCL = ch$limits[["LCL"]], UCL = ch$limits[["UCL"]],
               LCL_exact = exact_limits(alpha)[1], UCL_exact = exact_limits(alpha)[2],
               ARL0 = ch$arl$mean, ARL0_se = ch$arl$sd / sqrt(2000), ARL0_exact = 1 / alpha)
  }))
}))

## Part B: unknown parameters (Phase I data), conditional false-alarm rate --
alpha <- 0.0027; m <- 20; N <- m * n; prob <- 0.9; S <- 60
qL <- function(a) qgamma(a / 2, n, scale = 1 / n)
qU <- function(a) qgamma(1 - a / 2, n, scale = 1 / n)
# exact conditional FAP of limits built from a Phase I mean mu_hat (true mean 1)
fap_exact <- function(lim) pgamma(lim[1] / cc, n, scale = 1 / n) + 1 - pgamma(lim[2] / cc, n, scale = 1 / n)
fap_rho <- function(rho, a) pgamma(rho * qL(a), n, scale = 1 / n) + 1 - pgamma(rho * qU(a), n, scale = 1 / n)
mean_fap <- function(a) integrate(function(r) fap_rho(r, a) * dgamma(r, N, scale = 1 / N), 0, Inf,
                                  rel.tol = 1e-10)$value
exact_mean_nominal <- mean_fap(alpha)
exact_alpha_mean <- uniroot(function(a) mean_fap(a) - alpha, c(alpha / 100, alpha), tol = 1e-12)$root
# exact distribution of the conditional FAP and its 90% quantile at the nominal level
rho_q <- qgamma(1 - prob, N, scale = 1 / N)      # FAP decreases in rho, so the 90% FAP quantile is at rho_q
exact_q90_nominal <- fap_rho(rho_q, alpha)
exact_alpha_q <- uniroot(function(a) fap_rho(rho_q, a) - alpha, c(alpha / 100, alpha), tol = 1e-12)$root

one <- function(s) {
  set.seed(1000 + s)
  ph <- matrix(rexp(N), ncol = n)
  ch <- suppressWarnings(bootqcc_phase1(ph, "exp", u = u, alpha = alpha, fit_fn = fit,
                                        B = 4000, R = 100, N_ref = 20000, adjust = "mean",
                                        prob = prob))
  lim <- function(a) cc * 0 + stats::quantile(ch$boot_stats, c(a / 2, 1 - a / 2), names = FALSE)
  c(fap_none = fap_exact(lim(alpha)), fap_mean = fap_exact(lim(ch$alpha_mean)),
    fap_quant = fap_exact(lim(ch$alpha_quantile)),
    alpha_mean = ch$alpha_mean, alpha_quantile = ch$alpha_quantile)
}
res <- do.call(rbind, parallel::mclapply(seq_len(S), one, mc.cores = cores))
summ <- data.frame(
  method = c("none (nominal alpha)", "mean-adjusted", "quantile-adjusted (90%)"),
  mean_FAP = colMeans(res[, c("fap_none", "fap_mean", "fap_quant")]),
  prop_FAP_le_alpha = colMeans(res[, c("fap_none", "fap_mean", "fap_quant")] <= alpha))
partB <- list(S = S, m = m, n = n, alpha = alpha, summary = summ,
              exact_mean_FAP_nominal = exact_mean_nominal,
              exact_alpha_mean = exact_alpha_mean, exact_alpha_quantile = exact_alpha_q,
              exact_q90_FAP_nominal = exact_q90_nominal,
              est_alpha_mean = mean(res[, "alpha_mean"]), est_alpha_quantile = mean(res[, "alpha_quantile"]),
              raw = res)
saveRDS(list(partA = partA, partB = partB), "inst/validation/exponential-validation.rds", version = 2)
print(partA, digits = 4); print(summ, digits = 4)
str(partB[-c(5, 11)])
