test_that("maximum likelihood recovers the baseline model", {
  set.seed(11)
  x <- rweibull(400, shape = 2, scale = 3)
  f <- fit_bkwg(x, baseline = "weibull", submodel = "G")
  ref <- optim(c(1, 1), function(p) -sum(dweibull(x, exp(p[1]), exp(p[2]),
                                                  log = TRUE)))
  expect_equal(f$loglik, -ref$value, tolerance = 1e-4)
  expect_equal(unname(coef(f)), c(2, 3), tolerance = 0.15)
  expect_true(all(is.finite(f$se)))
})

test_that("sub-models are nested in terms of log-likelihood", {
  set.seed(12)
  x <- rbkwg(300, a = 3, b = 0.5, alpha = 2, beta = 0.6, par = c(1, 1))
  fG <- fit_bkwg(x, submodel = "G", hessian = FALSE)
  fK <- fit_bkwg(x, submodel = "KwG", hessian = FALSE)
  fB <- fit_bkwg(x, submodel = "BKwG", hessian = FALSE, nstart = 5,
                 start = coef(fK))
  expect_gte(fK$loglik, fG$loglik - 1e-6)
  expect_gte(fB$loglik, fK$loglik - 1e-3)
  expect_equal(fB$k, 6)
  expect_equal(fK$k, 4)
  expect_equal(unname(fK$estimate[c("alpha", "beta")]), c(1, 1))
})

test_that("all estimation methods run and give sensible fits", {
  set.seed(13)
  x <- rbkwg(200, a = 2, b = 1.5, alpha = 1, beta = 1, baseline = "exp",
             par = c(rate = 1))
  for (m in c("mle", "mps", "ad", "cvm", "lse", "wlse")) {
    f <- fit_bkwg(x, baseline = "exp", submodel = "KwG", method = m,
                  fixed = c(rate = 1))
    expect_s3_class(f, "bkwgfit")
    expect_equal(unname(coef(f)["a"]), 2, tolerance = 0.5, info = m)
  }
})

test_that("methods and inference functions work", {
  set.seed(14)
  x <- rbkwg(150, a = 2, b = 1, alpha = 3, beta = 0.8, par = c(1.2, 1))
  f <- fit_bkwg(x, submodel = "BKwG", fixed = c(scale = 1))
  expect_s3_class(summary(f), "summary.bkwgfit")
  expect_equal(nrow(confint(f)), f$k)
  expect_equal(attr(logLik(f), "df"), f$k)
  expect_equal(AIC(f), f$AIC)
  expect_length(residuals(f), length(x))
  expect_length(quantile(f), 3)
  g <- gof_bkwg(f)
  expect_s3_class(g, "bkwggof")
  expect_equal(nrow(g$table), 3)
  l <- lrt_bkwg(f, "KwG")
  expect_s3_class(l, "htest")
  expect_equal(unname(l$parameter), 2)
  w <- wald_bkwg(f, c(alpha = 1, beta = 1))
  expect_s3_class(w, "htest")
  tab <- compare_bkwg(x, submodels = c("KwG", "G"), fixed = c(scale = 1))
  expect_equal(nrow(tab), 2)
  pdf(NULL)
  expect_silent(plot(f))
  dev.off()
})

test_that("uniform baseline (no baseline parameters) can be fitted", {
  set.seed(15)
  x <- rbkwg(200, a = 2, b = 3, alpha = 1.5, beta = 1, baseline = "unif")
  f <- fit_bkwg(x, baseline = "unif", submodel = "EKwG")
  expect_equal(f$k, 3)
  expect_true(is.finite(f$loglik))
})
