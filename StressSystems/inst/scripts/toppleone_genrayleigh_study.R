# Monte Carlo study: Topp-Leone strength, generalized Rayleigh stress.
# Same design as the original TPL_GRL_SS_MODEL scripts.
library(StressSystems)

X <- ss_dist("toppleone")
Y <- ss_dist("genrayleigh")

grid_x <- data.frame(shape = c(0.3, 0.5, 0.7, 0.9))

# The original grid was (scale, shape) in VGAM's rgenray parametrisation.
# Here the generalized Rayleigh uses (shape, rate) with rate = 1 / scale.
orig <- matrix(c(0.5, 0.5,
                 1,   1,
                 0.5, 1,
                 2,   2,
                 3,   3.5,
                 4.5, 1), ncol = 2, byrow = TRUE)
grid_y <- data.frame(shape = orig[, 2], rate = 1 / orig[, 1])

res <- ss_simulate_grid(X, Y, grid_x, grid_y,
                        n = c(10, 25, 50), MC = 1000,
                        methods = c("mle", "jackknife", "bayes"),
                        seed = 2026, cores = 1)   # cores > 1 on Linux/macOS

write.csv(res$summary, "R_summary.csv", row.names = FALSE)
write.csv(res$parameters, "parameter_summary.csv", row.names = FALSE)

# Systems with the same families
ss_multistress(X, Y, par_x = c(shape = 0.7),
               par_y = data.frame(shape = c(1, 2), rate = c(2, 3)))
ss_R3(X, Y, X, par_z = c(shape = 0.3), par_y = c(shape = 1, rate = 2),
      par_x = c(shape = 0.9))
ss_sk(X, Y, s = 2, k = 3, par_x = c(shape = 0.7), par_y = c(shape = 1, rate = 2))
ss_cascade(X, Y, n = 3, c = 0.8, par_x = c(shape = 0.7),
           par_y = c(shape = 1, rate = 2))
