test_that("safe_profile_likelihood returns a finite, non-degenerate interval", {
  set.seed(2)
  d0 <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5, rate = 1))
  y <- rutd(80, d0)
  fit <- fit_utd_mle(y, dgamma, start = c(shape = 1, rate = 1),
                      lower = c(shape = 1e-6, rate = 1e-6))
  mle <- unlist(fit$params)
  nll <- function(theta) fit$negloglik(theta)
  ci <- safe_profile_likelihood(nll, mle, "shape", param_lower = 0)
  expect_true(is.finite(ci$lower) && is.finite(ci$upper))
  expect_true(ci$lower < ci$mle && ci$mle < ci$upper)
})

test_that("the corrected profile likelihood achieves closer-to-nominal coverage than the naive fixed-nuisance version", {
  skip_on_cran()  # this is a Monte Carlo demonstration, not a fast unit check
  set.seed(1)
  true_dist <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5, rate = 1))
  res <- compare_profile_likelihood_coverage(
    true_dist, c(shape = 2.5, rate = 1), dgamma,
    start = c(shape = 1, rate = 1), param_name = "shape",
    n = 60, n_sim = 150, lower = c(shape = 1e-6, rate = 1e-6), seed = 1)
  naive_cov <- res$coverage_95[res$method == "naive_fixed_nuisance"]
  corrected_cov <- res$coverage_95[res$method == "corrected_reoptimized"]
  expect_true(corrected_cov > naive_cov)
  expect_true(corrected_cov > 0.85)   # close to nominal 0.95
})
