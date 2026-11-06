test_that("moments agree with direct integration", {
  m1 <- mbkwg(1, 2, 1, 1.5, 2, par = c(2, 1))
  d1 <- integrate(function(x) x * dbkwg(x, 2, 1, 1.5, 2, par = c(2, 1)),
                  0, Inf)$value
  expect_equal(m1, d1, tolerance = 1e-6)
  ## baseline case: Weibull moments
  expect_equal(mbkwg(1:2, par = c(2, 3)),
               3^(1:2) * gamma(1 + (1:2) / 2), tolerance = 1e-6)
  s <- bkwg_stats(par = c(1, 2))           # exponential with mean 2
  expect_equal(unname(s[c("mean", "variance", "skewness", "kurtosis")]),
               c(2, 4, 2, 9), tolerance = 1e-5)
})

test_that("mgf, entropy and mean residual life of the exponential case", {
  expect_equal(mgf_bkwg(0.25, baseline = "exp", par = c(rate = 1)),
               1 / (1 - 0.25), tolerance = 1e-6)
  expect_equal(entropy_bkwg(baseline = "exp", par = c(rate = 2)),
               1 - log(2), tolerance = 1e-6)
  expect_equal(entropy_bkwg(baseline = "norm", par = c(0, 1), type = "renyi",
                            delta = 2), log(2 * sqrt(pi)), tolerance = 1e-6)
  expect_equal(mrl_bkwg(c(0.5, 3), baseline = "exp", par = c(rate = 0.5)),
               c(2, 2), tolerance = 1e-5)
})

test_that("order statistics", {
  val <- integrate(dos_bkwg, 0, Inf, i = 3, N = 7, a = 2, b = 1, alpha = 1.5,
                   beta = 2, par = c(2, 1))$value
  expect_equal(val, 1, tolerance = 1e-6)
  expect_equal(pos_bkwg(1, 7, 7, 2, 1, 1.5, 2, par = c(2, 1)),
               pbkwg(1, 2, 1, 1.5, 2, par = c(2, 1))^7)
})
