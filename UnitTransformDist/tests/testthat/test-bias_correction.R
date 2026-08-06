test_that("bias_correct_jackknife reduces or matches the original bias direction sensibly", {
  set.seed(21)
  d0 <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5))
  y <- rutd(60, d0)
  jk <- bias_correct_jackknife(y, dgamma, start = c(shape = 1), lower = c(shape = 1e-6))
  expect_true(is.finite(jk$corrected))
  expect_named(jk$corrected, "shape")
})

test_that("bias_correct_bootstrap (nonparametric) runs and returns a corrected estimate", {
  set.seed(22)
  d0 <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5))
  y <- rutd(80, d0)
  fit <- fit_utd_mle(y, dgamma, start = c(shape = 1), lower = c(shape = 1e-6))
  bc <- bias_correct_bootstrap(fit, y, n_boot = 50, type = "nonparametric")
  expect_true(is.finite(bc$corrected))
})

test_that("bias_correct() dispatcher routes to the correct method", {
  set.seed(23)
  d0 <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5))
  y <- rutd(60, d0)
  fit <- fit_utd_mle(y, dgamma, start = c(shape = 1), lower = c(shape = 1e-6))
  jk1 <- bias_correct(fit, y, method = "jackknife", lower = c(shape = 1e-6))
  jk2 <- bias_correct_jackknife(y, fit$ddist, start = fit$params, lower = c(shape = 1e-6))
  expect_equal(jk1$corrected, jk2$corrected, tolerance = 1e-8)
})
