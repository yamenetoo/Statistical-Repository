kw <- list(
  name = "kumaraswamy",
  d = function(x, a, b) a * b * x^(a - 1) * (1 - x^a)^(b - 1),
  q = function(p, a, b) (1 - (1 - p)^(1 / b))^(1 / a),
  r = function(n, a, b) (1 - (1 - runif(n))^(1 / b))^(1 / a)
)

test_that("works with a built-in distribution", {
  ch <- bootqcc("weibull", list(shape = 2, scale = 1), u = 0.5, n = 8,
                alpha = 0.05, B = 100, M = 3, arl_rep = 20, seed = 1)
  expect_s3_class(ch, "bootqcc")
  expect_lt(ch$limits[["LCL"]], ch$limits[["UCL"]])
  # The true median of Weibull(2, 1) lies inside the limits
  med <- qweibull(0.5, 2, 1)
  expect_true(med > ch$limits[["LCL"]] && med < ch$limits[["UCL"]])
  expect_length(ch$arl$runs, 20)
  expect_output(print(ch), "Parametric bootstrap control chart")
})

test_that("works with a user-defined distribution", {
  ch <- bootqcc(kw, list(a = 2, b = 3), u = 0.5, n = 15,
                alpha = 0.05, B = 100, M = 3, arl_rep = 20, seed = 2)
  expect_equal(ch$distr_name, "kumaraswamy")
  expect_lt(ch$limits[["LCL"]], ch$limits[["UCL"]])
  expect_true(all(ch$limits > 0 & ch$limits < 1))
})

test_that("in-control ARL is close to 1/alpha for a moderate alpha", {
  skip_on_cran()
  ch <- bootqcc("exp", list(rate = 1), u = 0.5, n = 10, alpha = 0.1,
                B = 300, M = 5, arl_rep = 300, seed = 3)
  expect_gt(ch$arl$mean, 5)
  expect_lt(ch$arl$mean, 20)
})

test_that("a scale shift shortens the ARL", {
  ch <- bootqcc("weibull", list(shape = 2, scale = 1), u = 0.5, n = 8,
                alpha = 0.05, B = 100, M = 3, arl_rep = 50, seed = 4)
  shifted <- arl(ch, params = list(shape = 2, scale = 2), rep = 50)
  expect_lt(shifted$mean, ch$arl$mean)
})

test_that("custom statistic and custom fitter are supported", {
  mean_stat <- function(est) est$scale * gamma(1 + 1 / est$shape)
  ch <- bootqcc("weibull", list(shape = 2, scale = 1), statistic = mean_stat,
                n = 8, alpha = 0.05, B = 100, M = 2, arl_rep = 10, seed = 5)
  expect_true(is.na(ch$u))
  expect_lt(ch$limits[["LCL"]], ch$limits[["UCL"]])

  exp_fit <- function(x) list(rate = 1 / mean(x))
  ch2 <- bootqcc("exp", list(rate = 1), u = 0.5, n = 10, fit_fn = exp_fit,
                 alpha = 0.05, B = 100, M = 2, arl_rep = 10, seed = 6)
  expect_lt(ch2$limits[["LCL"]], ch2$limits[["UCL"]])
})

test_that("monitoring flags clearly shifted samples", {
  ch <- bootqcc("weibull", list(shape = 2, scale = 1), u = 0.5, n = 8,
                alpha = 0.05, B = 100, M = 3, arl_rep = 10, seed = 7)
  set.seed(8)
  bad <- matrix(rweibull(8 * 5, 2, 5), ncol = 8)
  res <- bootqcc_monitor(ch, bad)
  expect_s3_class(res, "bootqcc_monitor")
  expect_true(all(res$status == "out of control"))
})

test_that("grid returns one row per setting", {
  g <- expand.grid(shape = c(1.5, 2), scale = 1)
  res <- bootqcc_grid("weibull", g, u_values = c(0.1, 0.9), n_values = 6,
                      alpha = 0.1, B = 50, M = 2, arl_rep = 5, seed = 9)
  expect_equal(nrow(res), 4)
  expect_true(all(c("MLCL", "MUCL", "MARL", "SDARL") %in% names(res)))
})

test_that("bad input gives clear errors", {
  expect_error(bootqcc("nosuchdist", list(a = 1)), "Could not find")
  expect_error(bootqcc("weibull", c(2, 1)), "named list")
  expect_error(bootqcc("weibull", list(shape = 2, scale = 1), u = 2), "'u'")
  expect_error(bootqcc(42, list(a = 1)), "distribution name")
})
