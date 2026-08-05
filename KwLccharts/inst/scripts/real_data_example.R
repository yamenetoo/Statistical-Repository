library(KwLccharts)
set.seed(777)

servtime <- c(
1.326, 0.841, 0.282, 2.830, 0.121,
0.644, 0.197, 1.581, 2.178, 1.553,
1.447, 2.343, 0.863, 3.658, 0.132,
4.033, 3.978, 2.416, 0.534, 0.501,
0.458, 4.003, 0.260, 1.099, 0.696,
0.164, 1.271, 0.641, 3.743, 0.395,
0.203, 0.296, 0.529, 1.485, 2.825,
0.115, 1.589, 2.444, 1.219, 3.578,
0.540, 0.507, 0.466, 0.047, 0.334)

fit <- fit_KwL(servtime, start = c(theta = 1, a = 1, b = 1))
cat("MLE (theta, a, b) =", fit$theta, fit$a, fit$b, "\n")
cat("log-likelihood =", -fit$value, "\n")

k <- length(fit_KwL) # dummy
aic <- 2*3 - 2*(-fit$value)
n_obs <- length(servtime)
bic <- 3*log(n_obs) - 2*(-fit$value)
cat("AIC =", aic, " BIC =", bic, "\n")

# Kolmogorov-Smirnov goodness of fit
ks <- ks.test(servtime, function(x) pKwL(x, fit$theta, fit$a, fit$b))
print(ks)

# subgroups of size 5 (same grouping as reference paper Table 10)
subgroups <- matrix(servtime, nrow = 9, ncol = 5, byrow = TRUE)
u_target <- 0.05
Qhat_sub <- numeric(9)
for (i in 1:9) {
  f <- fit_KwL(subgroups[i, ], start = c(theta = fit$theta, a = fit$a, b = fit$b))
  Qhat_sub[i] <- if (!is.null(f)) qKwL(u_target, f$theta, f$a, f$b) else NA
}
print(data.frame(subgroup = 1:9, Qhat = Qhat_sub))
write.csv(data.frame(subgroup = 1:9, subgroups, Qhat = Qhat_sub),
          "real_data_subgroups.csv", row.names = FALSE)

# Phase-I control limits from the fitted KwL parameters, n = 5, alpha = 0.0027
res <- KwLBootCC(u = u_target, alpha = 0.0027, theta = fit$theta, a = fit$a, b = fit$b,
                  n = 5, B = 200, MC = 8, start = c(theta = fit$theta, a = fit$a, b = fit$b))
cat(sprintf("Phase I limits: MLCL=%.6f SDLCL=%.6f MUCL=%.6f SDUCL=%.6f\n",
            res$MLCL, res$SDLCL, res$MUCL, res$SDUCL))
saveRDS(list(fit = fit, ks = ks, aic = aic, bic = bic, Qhat_sub = Qhat_sub, limits = res),
        "real_data_results.rds")
write.csv(data.frame(MLCL = res$MLCL, SDLCL = res$SDLCL, MUCL = res$MUCL, SDUCL = res$SDUCL),
          "real_data_limits.csv", row.names = FALSE)
