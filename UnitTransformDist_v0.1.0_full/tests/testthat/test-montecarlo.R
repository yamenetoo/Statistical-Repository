test_that("run_monte_carlo recovers approximately unbiased estimates for a known case", {
  set.seed(11)
  d0 <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 3))
  mc <- run_monte_carlo(d0, true_params = c(shape = 3), ddist = dgamma,
                         start = c(shape = 1), lower = c(shape = 1e-6),
                         n = 150, n_sim = 100, seed = 11)
  expect_s3_class(mc, "utd_mc_sim")
  expect_true(abs(mc$bias[["shape"]]) < 0.5)
  expect_equal(mc$n_sim, 100)
})

test_that("summary.utd_mc_sim returns a printable summary object without side-effect printing", {
  set.seed(5)
  d0 <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2))
  mc <- run_monte_carlo(d0, true_params = c(shape = 2), ddist = dgamma,
                         start = c(shape = 1), lower = c(shape = 1e-6),
                         n = 80, n_sim = 40, seed = 5)
  s <- summary(mc)
  expect_s3_class(s, "summary.utd_mc_sim")
  expect_true(is.data.frame(s$table))
})

test_that("compare_estimators_mc returns one row per estimator per parameter", {
  set.seed(3)
  d0 <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2))
  specs <- list(mle = list(ddist = dgamma, start = c(shape = 1), lower = c(shape = 1e-6)))
  res <- compare_estimators_mc(d0, c(shape = 2), specs, n = 60, n_sim = 30, seed = 3)
  expect_equal(nrow(res), 1)
  expect_true(all(c("estimator", "param", "bias", "rmse") %in% names(res)))
})
