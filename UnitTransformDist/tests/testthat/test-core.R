test_that("dutd/putd/qutd are mutually consistent (analytic contract)", {
  d0 <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5))
  y <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  # F(Q(p)) == p for the CDF/quantile pair
  p <- c(0.1, 0.5, 0.9)
  expect_equal(putd(qutd(p, d0), d0), p, tolerance = 1e-6)
  # density integrates to 1 over (0,1)
  total <- stats::integrate(function(y) dutd(y, d0), 0, 1)$value
  expect_equal(total, 1, tolerance = 1e-4)
})

test_that("boundary behaviour: dutd is 0 (not NaN) outside (0,1)", {
  d0 <- utd(dgamma, pgamma, qgamma, params = list(shape = 2.5))
  expect_equal(dutd(c(-1, 0, 1, 2), d0), c(0, 0, 0, 0))
})

test_that("qutd works via numeric inversion when qdist is omitted", {
  d0 <- utd(dgamma, pgamma, params = list(shape = 2.5))  # no qdist
  q_analytic <- utd(dgamma, pgamma, qgamma, params = list(shape = 2.5))
  expect_equal(qutd(0.4, d0), qutd(0.4, q_analytic), tolerance = 1e-4)
})

test_that("utd() validates its inputs", {
  expect_error(utd("not a function", pgamma), "must be a function")
  expect_error(utd(dgamma, pgamma, params = "not a list"), "must be a named list")
})

test_that("fit_utd_mle rejects out-of-range data", {
  expect_error(fit_utd_mle(c(0.5, 1.2), dgamma, start = c(shape = 1)), "strictly inside")
  expect_error(fit_utd_mle(c(0.5, NA), dgamma, start = c(shape = 1)), "NA")
})

test_that("fit_utd_mle recovers a known parameter reasonably well", {
  set.seed(123)
  d0 <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 3))
  y <- rutd(500, d0)
  fit <- fit_utd_mle(y, dgamma, start = c(shape = 1), lower = c(shape = 1e-6))
  expect_equal(unname(fit$params$shape), 3, tolerance = 0.3)
})

test_that("profile_likelihood_ci does not collapse to the MLE (regression test for the historical bug)", {
  set.seed(7)
  d0 <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5))
  y <- rutd(300, d0)
  fit <- fit_utd_mle(y, dgamma, start = c(shape = 1), lower = c(shape = 1e-6))
  ci <- profile_likelihood_ci(fit, "shape")
  expect_true(is.finite(ci$lower) && is.finite(ci$upper))
  expect_true(ci$upper - ci$lower > 1e-6)   # not a single-point collapse
  expect_true(ci$lower < ci$mle && ci$mle < ci$upper)
})
