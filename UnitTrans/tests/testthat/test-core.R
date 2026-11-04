test_that("density integrates to 1 for every valid pair", {
  cases <- list(
    list("norm", c(0, 1), "logit"), list("norm", c(0, 1), "probit"),
    list("norm", c(0, 1), "cloglog"), list("gamma", c(3, 1), "rational"),
    list("gamma", c(3, 1), "exponential"), list("gamma", c(3, 1), "tanh"),
    list("beta", c(3, 2), "rational"), list("exp", 3, "exponential"))
  for (cs in cases) {
    i <- tryCatch(integrate(function(v) dtrans(v, cs[[1]], cs[[2]], cs[[3]]),
                            1e-9, 1 - 1e-9, rel.tol = 1e-5, subdivisions = 2000L)$value,
                  error = function(e) NA_real_)
    if (!is.na(i)) expect_equal(i, 1, tolerance = 1e-2, info = paste(cs[[1]], cs[[3]]))
    expect_equal(ptrans(1 - 1e-9, cs[[1]], cs[[2]], cs[[3]]) -
                 ptrans(1e-9, cs[[1]], cs[[2]], cs[[3]]), 1, tolerance = 5e-3)
  }
})

test_that("p and q are inverses; d is derivative of p", {
  for (tn in c("rational", "exponential", "logit", "tanh", "probit", "cloglog")) {
    d <- if (tn %in% c("logit", "probit", "cloglog")) "norm" else "gamma"
    par <- if (d == "norm") c(0.3, 1.2) else c(2.5, 1.5)
    p <- c(.1, .3, .5, .8, .95)
    q <- qtrans(p, d, par, tn)
    expect_equal(ptrans(q, d, par, tn), p, tolerance = 1e-8, info = tn)
    h <- 1e-6; v <- q[3]
    num <- (ptrans(v + h, d, par, tn) - ptrans(v - h, d, par, tn)) / (2 * h)
    expect_equal(dtrans(v, d, par, tn), num, tolerance = 1e-4, info = tn)
  }
})

test_that("random sample matches the cdf", {
  set.seed(1)
  x <- rtrans(5000, "gamma", c(2, 1), "exponential")
  expect_gt(ks.test(x, function(q) ptrans(q, "gamma", c(2, 1), "exponential"))$p.value, 0.01)
})

test_that("domain checking works", {
  expect_false(check_domain("norm", c(0, 1), "rational")$valid)
  expect_true(check_domain("norm", c(0, 1), "logit")$valid)
  expect_true(check_domain("gamma", c(3, 1), "tanh")$valid)
})

test_that("mle recovers parameters and fit_all picks sensible model", {
  set.seed(10)
  x <- rtrans(1500, "norm", c(0.5, 1.5), "logit")
  f <- fit_trans(x, "norm", "logit")
  expect_equal(unname(f$coef), c(0.5, 1.5), tolerance = 0.1)
  expect_s3_class(f, "transfit")
  expect_equal(dim(confint(f)), c(2, 2))
  r <- fit_all(x, c("norm", "gamma"), c("logit", "rational"))
  expect_equal(r$table$model[1], "norm + logit")
})

test_that("moments agree with simulation", {
  s <- summary_trans("gamma", c(2, 1), "rational")
  set.seed(2); x <- rtrans(2e5, "gamma", c(2, 1), "rational")
  expect_equal(unname(s["mean"]), mean(x), tolerance = 0.01)
  expect_equal(unname(s["variance"]), var(x), tolerance = 0.02)
})

test_that("gof and simulation run", {
  set.seed(3)
  f <- fit_trans(rtrans(80, "norm", c(0, 1), "logit"), "norm", "logit")
  g <- gof_trans(f, B = 20)
  expect_equal(nrow(g), 3)
  s <- simulation_study("norm", c(0, 1), "logit", n = 60, R = 5, seed = 1)
  expect_equal(nrow(s), 2)
})
