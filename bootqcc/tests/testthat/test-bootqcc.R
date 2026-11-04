kw <- list(
  name = "kumaraswamy",
  d = function(x, a, b) a * b * x^(a - 1) * (1 - x^a)^(b - 1),
  q = function(p, a, b) (1 - (1 - p)^(1 / b))^(1 / a),
  r = function(n, a, b) (1 - (1 - runif(n))^(1 / b))^(1 / a)
)

test_that("works with a built-in distribution", {
  ch <- bootqcc("weibull", list(shape = 2, scale = 1), u = 0.5, n = 8,
                alpha = 0.05, B = 400, M = 3, arl_rep = 20, seed = 1)
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
                alpha = 0.05, B = 400, M = 3, arl_rep = 20, seed = 2)
  expect_equal(ch$distr_name, "kumaraswamy")
  expect_lt(ch$limits[["LCL"]], ch$limits[["UCL"]])
  expect_true(all(ch$limits > 0 & ch$limits < 1))
})

test_that("in-control ARL is close to 1/alpha for a moderate alpha", {
  skip_on_cran()
  ch <- bootqcc("exp", list(rate = 1), u = 0.5, n = 10, alpha = 0.1,
                B = 600, M = 5, arl_rep = 300, seed = 3)
  expect_gt(ch$arl$mean, 5)
  expect_lt(ch$arl$mean, 20)
})

test_that("a scale shift shortens the ARL", {
  ch <- bootqcc("weibull", list(shape = 2, scale = 1), u = 0.5, n = 8,
                alpha = 0.05, B = 400, M = 3, arl_rep = 50, seed = 4)
  shifted <- arl(ch, params = list(shape = 2, scale = 2), rep = 50)
  expect_lt(shifted$mean, ch$arl$mean)
})

test_that("custom statistic and custom fitter are supported", {
  mean_stat <- function(est) est$scale * gamma(1 + 1 / est$shape)
  ch <- bootqcc("weibull", list(shape = 2, scale = 1), statistic = mean_stat,
                n = 8, alpha = 0.05, B = 400, M = 2, arl_rep = 10, seed = 5)
  expect_true(is.na(ch$u))
  expect_lt(ch$limits[["LCL"]], ch$limits[["UCL"]])

  exp_fit <- function(x) list(rate = 1 / mean(x))
  ch2 <- bootqcc("exp", list(rate = 1), u = 0.5, n = 10, fit_fn = exp_fit,
                 alpha = 0.05, B = 400, M = 2, arl_rep = 10, seed = 6)
  expect_lt(ch2$limits[["LCL"]], ch2$limits[["UCL"]])
})

test_that("monitoring flags clearly shifted samples", {
  ch <- bootqcc("weibull", list(shape = 2, scale = 1), u = 0.5, n = 8,
                alpha = 0.05, B = 400, M = 3, arl_rep = 10, seed = 7)
  set.seed(8)
  bad <- matrix(rweibull(8 * 5, 2, 5), ncol = 8)
  res <- bootqcc_monitor(ch, bad)
  expect_s3_class(res, "bootqcc_monitor")
  expect_true(all(res$status == "out of control"))
})

test_that("grid returns one row per setting", {
  g <- expand.grid(shape = c(1.5, 2), scale = 1)
  res <- bootqcc_grid("weibull", g, u_values = c(0.1, 0.9), n_values = 6,
                      alpha = 0.1, B = 200, M = 2, arl_rep = 5, seed = 9)
  expect_equal(nrow(res), 4)
  expect_true(all(c("MLCL", "MUCL", "MARL", "SDARL") %in% names(res)))
})

test_that("bad input gives clear errors", {
  expect_error(bootqcc("nosuchdist", list(a = 1)), "Could not find")
  expect_error(bootqcc("weibull", c(2, 1)), "named list")
  expect_error(bootqcc("weibull", list(shape = 2, scale = 1), u = 2), "'u'")
  expect_error(bootqcc(42, list(a = 1)), "distribution name")
})

test_that("small B gives a warning about tail resolution", {
  expect_warning(bootqcc("exp", list(rate = 1), n = 5, alpha = 0.0027, B = 100, M = 1,
                         arl_rep = 2, fit_fn = function(x) list(rate = 1 / mean(x))),
                 "too few bootstrap replicates")
})

test_that("limits approach the exact exponential limits", {
  fit <- function(x) list(rate = 1 / mean(x))
  ch <- bootqcc("exp", list(rate = 1), u = 0.5, n = 5, alpha = 0.1, B = 4000, M = 5,
                arl_rep = 2, fit_fn = fit, seed = 21)
  ex <- -log(0.5) * qgamma(c(0.05, 0.95), 5, scale = 1 / 5)
  expect_equal(unname(ch$limits), ex, tolerance = 0.05)
})

test_that("bootqcc_phase1 estimates parameters and adjusts the level", {
  set.seed(31)
  fit <- function(x) list(rate = 1 / mean(x))
  ph <- matrix(rexp(10 * 5), ncol = 5)
  ch <- bootqcc_phase1(ph, "exp", u = 0.5, alpha = 0.1, fit_fn = fit, B = 400, R = 30,
                       N_ref = 2000, adjust = "mean", seed = 32)
  expect_s3_class(ch, c("bootqcc_phase1", "bootqcc"))
  expect_equal(ch$params_hat$rate, 1 / mean(ph), tolerance = 1e-8)
  expect_lte(ch$alpha_adj, 0.1)
  expect_lte(ch$limits[["LCL"]], ch$limits_unadjusted[["LCL"]])
  expect_gte(ch$limits[["UCL"]], ch$limits_unadjusted[["UCL"]])
  expect_output(print(ch), "Phase I data")
  # works with the Phase II tools
  res <- bootqcc_monitor(ch, matrix(rexp(5 * 3, rate = 0.05), ncol = 5))
  expect_s3_class(res, "bootqcc_monitor")
  expect_error(bootqcc_phase1(rexp(20), "exp", fit_fn = fit), "'n'")
})
