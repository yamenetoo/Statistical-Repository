library(KwLccharts)
set.seed(20260805)

B_lim  <- 300   # bootstrap iterations for control-limit estimation
MC_lim <- 10    # outer Monte Carlo repetitions for control-limit estimation
alpha  <- 0.0027
us     <- c(0.05, 0.10)
ns     <- c(4, 5, 6)

param_main <- c(theta = 1.0, a = 1.5, b = 1.5)
param_alt  <- c(theta = 1.5, a = 2.0, b = 1.0)

cat("=== PART 1: Average control limits (Table analogue) ===\n")
limits_results <- list()
idx <- 1
for (n in ns) {
  for (u in us) {
    res <- KwLBootCC(u = u, alpha = alpha,
                      theta = param_main["theta"], a = param_main["a"], b = param_main["b"],
                      n = n, B = B_lim, MC = MC_lim)
    limits_results[[idx]] <- data.frame(
      n = n, theta = param_main["theta"], a = param_main["a"], b = param_main["b"], u = u,
      MLCL = res$MLCL, SDLCL = res$SDLCL, MUCL = res$MUCL, SDUCL = res$SDUCL
    )
    cat(sprintf("n=%d u=%.2f (theta,a,b)=(%.2f,%.2f,%.2f)  MLCL=%.5f SDLCL=%.5f MUCL=%.5f SDUCL=%.5f\n",
                n, u, param_main["theta"], param_main["a"], param_main["b"],
                res$MLCL, res$SDLCL, res$MUCL, res$SDUCL))
    idx <- idx + 1
  }
}

# Alternative parameter set, only n = 5, to illustrate effect of KwL parameters
for (u in us) {
  res <- KwLBootCC(u = u, alpha = alpha,
                    theta = param_alt["theta"], a = param_alt["a"], b = param_alt["b"],
                    n = 5, B = B_lim, MC = MC_lim)
  limits_results[[idx]] <- data.frame(
    n = 5, theta = param_alt["theta"], a = param_alt["a"], b = param_alt["b"], u = u,
    MLCL = res$MLCL, SDLCL = res$SDLCL, MUCL = res$MUCL, SDUCL = res$SDUCL
  )
  cat(sprintf("n=%d u=%.2f (theta,a,b)=(%.2f,%.2f,%.2f)  MLCL=%.5f SDLCL=%.5f MUCL=%.5f SDUCL=%.5f\n",
              5, u, param_alt["theta"], param_alt["a"], param_alt["b"],
              res$MLCL, res$SDLCL, res$MUCL, res$SDUCL))
  idx <- idx + 1
}

limits_df <- do.call(rbind, limits_results)
saveRDS(limits_df, "limits_df.rds")
write.csv(limits_df, "limits_table.csv", row.names = FALSE)
cat("\nSaved limits_table.csv\n")

cat("\n=== PART 2: ARL simulation (in-control & out-of-control) ===\n")
n_runs  <- 25
max_run <- 600
arl_results <- list()
idx <- 1

for (u in us) {
  # Phase-I limits from the main parameter set at n = 5
  lim_row <- limits_df[limits_df$n == 5 & limits_df$u == u &
                          abs(limits_df$theta - param_main["theta"]) < 1e-9, ]
  LCL <- lim_row$MLCL; UCL <- lim_row$MUCL

  # In-control ARL (Phase-II generated from the SAME parameters)
  res_ic <- KwL_ARL(u, LCL, UCL, param_main["theta"], param_main["a"], param_main["b"],
                     n = 5, n_runs = n_runs, max_run = max_run)
  arl_results[[idx]] <- data.frame(u = u, shift = "in-control",
                                    theta1 = param_main["theta"], a1 = param_main["a"], b1 = param_main["b"],
                                    ARL = res_ic$ARL, SDRL = res_ic$SDRL)
  cat(sprintf("u=%.2f IN-CONTROL   ARL=%.3f SDRL=%.3f\n", u, res_ic$ARL, res_ic$SDRL))
  idx <- idx + 1

  # Out-of-control ARL (Phase-II generated from the ALTERNATIVE parameters)
  res_oc <- KwL_ARL(u, LCL, UCL, param_alt["theta"], param_alt["a"], param_alt["b"],
                     n = 5, n_runs = n_runs, max_run = max_run)
  arl_results[[idx]] <- data.frame(u = u, shift = "shifted",
                                    theta1 = param_alt["theta"], a1 = param_alt["a"], b1 = param_alt["b"],
                                    ARL = res_oc$ARL, SDRL = res_oc$SDRL)
  cat(sprintf("u=%.2f SHIFTED     ARL=%.3f SDRL=%.3f\n", u, res_oc$ARL, res_oc$SDRL))
  idx <- idx + 1
}

arl_df <- do.call(rbind, arl_results)
saveRDS(arl_df, "arl_df.rds")
write.csv(arl_df, "arl_table.csv", row.names = FALSE)
cat("\nSaved arl_table.csv\n")

cat("\n=== DONE ===\n")
