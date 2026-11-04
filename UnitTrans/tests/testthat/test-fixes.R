# ---- transformations --------------------------------------------------------
test_that("registry: tanh replaces arcsech, alias still works", {
  tl <- trans_list()
  expect_equal(nrow(tl), 6)
  expect_false("arcsech" %in% tl$name)
  expect_true("tanh" %in% tl$name)
  y <- c(0.1, 1, 3)
  expect_equal(forward_transformation(y, "tanh"), 1 - tanh(y / 2))
  expect_equal(inverse_transformation(c(.2, .7), "tanh"), 2 * atanh(1 - c(.2, .7)))
  expect_warning(a <- dtrans(.4, "gamma", c(2, 1), "arcsech"), "renamed")
  expect_equal(a, dtrans(.4, "gamma", c(2, 1), "tanh"))
})

test_that("log-Jacobians equal the numerical derivative of the inverse map", {
  v <- c(.05, .2, .5, .8, .95); h <- 1e-6
  for (tn in trans_list()$name) {
    tr <- get_transformation(tn)
    num <- abs((tr$ginv(v + h) - tr$ginv(v - h)) / (2 * h))
    expect_equal(exp(tr$logjac(v)), num, tolerance = 1e-5, info = tn)
  }
  expect_error(get_transformation("nope"), "Unknown transformation")
})

# ---- closed-form special cases ----------------------------------------------
test_that("known unit distributions are reproduced", {
  v <- c(.1, .3, .6, .9)
  # Exp(lambda) + exponential  ->  Beta(lambda, 1)
  expect_equal(dtrans(v, "exp", 3, "exponential"), dbeta(v, 3, 1))
  expect_equal(ptrans(v, "exp", 3, "exponential"), pbeta(v, 3, 1))
  expect_equal(qtrans(v, "exp", 3, "exponential"), qbeta(v, 3, 1), tolerance = 1e-8)
  # Gamma(a, b) + exponential -> unit-gamma
  a <- 2.5; b <- 1.7
  expect_equal(dtrans(v, "gamma", c(a, b), "exponential"),
               b^a / gamma(a) * (-log(v))^(a - 1) * v^(b - 1))
  # Normal + logit -> logit-normal
  mu <- .4; s <- 1.3
  expect_equal(dtrans(v, "norm", c(mu, s), "logit"),
               dnorm(qlogis(v), mu, s) / (v * (1 - v)))
  # probit with N(0,1), logit with logistic, exponential with Exp(1) -> uniform
  expect_equal(dtrans(v, "norm", c(0, 1), "probit"), rep(1, 4))
  expect_equal(dtrans(v, "logis", c(0, 1), "logit"), rep(1, 4))
  expect_equal(dtrans(v, "exp", 1, "exponential"), rep(1, 4))
})

# ---- tails and log.p --------------------------------------------------------
test_that("p/q are consistent for lower.tail and log.p, for all transformations", {
  for (tn in trans_list()$name) {
    d <- if (tn %in% c("logit", "probit", "cloglog")) "norm" else "gamma"
    par <- if (d == "norm") c(0.3, 1.2) else c(2.5, 1.5)
    p <- c(.01, .2, .5, .9, .99)
    q <- qtrans(p, d, par, tn)
    expect_equal(ptrans(q, d, par, tn, lower.tail = FALSE), 1 - p, tolerance = 1e-8, info = tn)
    expect_equal(ptrans(q, d, par, tn, log.p = TRUE), log(p), tolerance = 1e-8, info = tn)
    expect_equal(ptrans(q, d, par, tn, lower.tail = FALSE, log.p = TRUE), log1p(-p),
                 tolerance = 1e-8, info = tn)
    expect_equal(qtrans(log(p), d, par, tn, log.p = TRUE), q, tolerance = 1e-8, info = tn)
    expect_equal(qtrans(1 - p, d, par, tn, lower.tail = FALSE), q, tolerance = 1e-8, info = tn)
    expect_equal(ptrans(c(-1, 0, 1, 2), d, par, tn), c(0, 0, 1, 1), info = tn)
    expect_equal(ptrans(c(-1, 2), d, par, tn, lower.tail = FALSE), c(1, 0), info = tn)
  }
})

test_that("upper tail keeps precision where 1 - p is exactly 0", {
  v <- 1 - 1e-12
  y <- v / (1 - v)
  s <- ptrans(v, "lnorm", c(0, 1), "rational", lower.tail = FALSE)
  expect_gt(s, 0)
  expect_equal(s, plnorm(y, 0, 1, lower.tail = FALSE))
  expect_true(is.finite(hazard_trans(v, "lnorm", c(0, 1), "rational")))
})

test_that("hazard agrees with the closed form of Beta(lambda, 1)", {
  v <- c(.2, .5, .9)
  expect_equal(hazard_trans(v, "exp", 3, "exponential"), 3 * v^2 / (1 - v^3))
  expect_equal(survival_trans(v, "exp", 3, "exponential"), 1 - v^3)
})

# ---- moments ----------------------------------------------------------------
test_that("moments match Beta(3, 1)", {
  s <- summary_trans("exp", 3, "exponential")
  expect_equal(unname(s["mean"]), 3 / 4, tolerance = 1e-6)
  expect_equal(unname(s["variance"]), 3 / 80, tolerance = 1e-6)
  expect_equal(unname(moments_trans("exp", 3, "exponential", 2)), 3 / 5, tolerance = 1e-6)
})

# ---- baseline domain --------------------------------------------------------
test_that("check_domain measures the mass outside the support", {
  d <- check_domain("norm", c(0, 1), "rational")
  expect_false(d$valid); expect_equal(d$mass_outside, 0.5)
  expect_true(check_domain("norm", c(10, 1), "rational")$valid)   # mass ~ 1e-23
  expect_true(check_domain("gamma", c(0.2, 1), "exponential")$valid)  # tiny shape
  expect_true(check_domain("norm", c(0, 1), "logit")$valid)
  expect_error(fit_trans(runif(20, .1, .9), "norm", "rational"), "Incompatible")
  det <- detect_transformations("norm", c(0, 1))
  expect_equal(det$valid, c(FALSE, FALSE, TRUE, FALSE, TRUE, TRUE))
})

# ---- custom baselines without log / lower.tail / log.p ----------------------
test_that("a minimal user-defined baseline works", {
  assign("dminexp", function(x, lambda) lambda * exp(-lambda * x), envir = globalenv())
  assign("pminexp", function(q, lambda) 1 - exp(-lambda * q), envir = globalenv())
  assign("qminexp", function(p, lambda) -log(1 - p) / lambda, envir = globalenv())
  on.exit(rm(list = c("dminexp", "pminexp", "qminexp"), envir = globalenv()), add = TRUE)
  v <- c(.2, .5, .8)
  expect_equal(baseline_params("minexp"), "lambda")
  expect_equal(dtrans(v, "minexp", c(lambda = 2), "rational"),
               dtrans(v, "exp", 2, "rational"))
  expect_equal(ptrans(v, "minexp", 2, "exponential", lower.tail = FALSE, log.p = TRUE),
               ptrans(v, "exp", 2, "exponential", lower.tail = FALSE, log.p = TRUE))
  expect_equal(qtrans(c(.1, .9), "minexp", 2, "exponential", lower.tail = FALSE),
               qtrans(c(.1, .9), "exp", 2, "exponential", lower.tail = FALSE), tolerance = 1e-8)
  set.seed(1); x <- rtrans(300, "minexp", 2, "rational")
  f <- fit_trans(x, "minexp", "rational")
  expect_s3_class(f, "transfit")
})

# ---- fitting ----------------------------------------------------------------
test_that("fit_trans validates data and supports adjust", {
  x <- c(0, .2, .4, .6, 1)
  expect_error(fit_trans(x, "norm", "logit"), "strictly inside")
  set.seed(1); y <- c(0, 1, rtrans(100, "norm", c(0, 1), "logit"))
  f <- fit_trans(y, "norm", "logit", adjust = TRUE)
  expect_true(f$adjusted)
  expect_true(all(f$data > 0 & f$data < 1))
  expect_error(fit_trans(c(.1, .2), "norm", "logit"), "at least 3")
})

test_that("standard errors, methods and quantile residuals", {
  set.seed(11)
  x <- rtrans(2000, "norm", c(0.5, 1.5), "logit")
  f <- fit_trans(x, "norm", "logit")
  expect_true(f$hessian_ok)
  expect_true(all(is.finite(f$se)))
  expect_equal(nobs(f), 2000L)
  expect_equal(AIC(f), f$AIC)
  expect_equal(BIC(f), f$BIC, tolerance = 1e-8)
  r <- residuals(f)
  expect_lt(abs(mean(r)), 0.1); expect_equal(sd(r), 1, tolerance = 0.1)
  expect_equal(predict(f, c(.2, .5), "cdf"), ptrans(c(.2, .5), "norm", coef(f), "logit"))
  expect_equal(dim(confint(f, parm = "mean")), c(1, 2))
  expect_equal(rownames(confint(f, parm = "sd")), "sd")
})

test_that("positive-bounded parameters recover on the log scale and k = 1 works", {
  set.seed(12)
  x <- rtrans(1500, "exp", 2.5, "exponential")
  f <- fit_trans(x, "exp", "exponential")
  expect_equal(unname(f$coef), 2.5, tolerance = 0.1)
  x2 <- rtrans(1500, "gamma", c(0.6, 2), "exponential")
  f2 <- fit_trans(x2, "gamma", "exponential")
  expect_equal(unname(f2$coef), c(0.6, 2), tolerance = 0.2)
  expect_true(all(f2$coef > 0))
})

test_that("profile and Wald intervals agree in large samples and contain the truth", {
  set.seed(13)
  x <- rtrans(1500, "weibull", c(2, 1.5), "exponential")
  f <- fit_trans(x, "weibull", "exponential")
  w <- confint(f); p <- confint(f, method = "profile")
  expect_equal(unname(p), unname(w), tolerance = 0.05)
  expect_true(all(p[, 1] < f$coef & f$coef < p[, 2]))
  expect_true(2 > p[1, 1] && 2 < p[1, 2])
  expect_true(1.5 > p[2, 1] && 1.5 < p[2, 2])
})

test_that("seeds are reproducible and the caller's RNG stream is restored", {
  set.seed(99); a <- runif(1)
  set.seed(99)
  x <- rtrans(120, "norm", c(0, 1), "logit")
  set.seed(99); invisible(runif(0))
  f1 <- fit_trans(x, "norm", "logit", seed = 5)
  g1 <- gof_trans(f1, B = 10, seed = 3)
  s1 <- simulation_study("norm", c(0, 1), "logit", n = 60, R = 3, seed = 4)
  b1 <- confint(f1, method = "bootstrap", B = 10, seed = 2)
  set.seed(99); b <- runif(1)
  expect_equal(a, b)                      # nothing above consumed the stream
  expect_equal(gof_trans(f1, B = 10, seed = 3), g1)
  expect_equal(confint(f1, method = "bootstrap", B = 10, seed = 2), b1)
  expect_equal(simulation_study("norm", c(0, 1), "logit", n = 60, R = 3, seed = 4), s1)
})

test_that("bootstrap interval reports failures and has the right shape", {
  set.seed(14)
  f <- fit_trans(rtrans(150, "norm", c(0, 1), "logit"), "norm", "logit")
  ci <- confint(f, method = "bootstrap", B = 20, seed = 1)
  expect_equal(dim(ci), c(2, 2))
  expect_true(!is.null(attr(ci, "failed")))
  expect_true(all(ci[, 1] < ci[, 2]))
})

# ---- model search -----------------------------------------------------------
test_that("fit_all ranks with weights, is reproducible and parallel-safe", {
  set.seed(15)
  x <- rtrans(400, "gamma", c(3, 2), "rational")
  r1 <- fit_all(x, c("gamma", "lnorm"), c("rational", "exponential"), seed = 1)
  expect_equal(sum(r1$table$weight), 1)
  expect_true(all(r1$table$converged))
  expect_equal(r1$table$delta[1], 0)
  expect_equal(r1$table$model[1], "gamma + rational")
  expect_equal(best_model(r1)$dist, "gamma")
  r2 <- fit_all(x, c("gamma", "lnorm"), c("rational", "exponential"), seed = 1)
  expect_equal(r1$table, r2$table)
  skip_on_os("windows")
  r3 <- fit_all(x, c("gamma", "lnorm"), c("rational", "exponential"), seed = 1, cores = 2)
  expect_equal(r1$table, r3$table, tolerance = 1e-8)
  expect_error(fit_all(x, "norm", "rational"), "No valid combination")
})

# ---- goodness of fit --------------------------------------------------------
test_that("gof p-values are never zero and reject a wrong model", {
  set.seed(16)
  x <- rtrans(400, "norm", c(0.5, 1.5), "logit")
  wrong <- fit_trans(x, "cauchy", "probit")
  g <- gof_trans(wrong, B = 29, seed = 1)
  expect_true(all(g$p.value >= 1 / 30))
  expect_true(all(g$p.value <= 1))
  expect_lt(min(g$p.value), 0.1)
  good <- fit_trans(x, "norm", "logit")
  g2 <- gof_trans(good, B = 29, seed = 1)
  expect_gt(g2$p.value[g2$test == "CvM"], 0.05)
})

# ---- simulation study -------------------------------------------------------
test_that("simulation study shows small bias and sensible coverage", {
  s <- simulation_study("norm", c(mean = 0.3, sd = 1.2), "logit", n = 300, R = 60, seed = 1)
  expect_equal(nrow(s), 2)
  expect_true(all(abs(s$bias) < 0.1))
  expect_true(all(s$coverage > 0.85))
  expect_true(all(s$failed == 0))
})

# ---- plot -------------------------------------------------------------------
test_that("plot works for all panels and subsets", {
  set.seed(17)
  f <- fit_trans(rtrans(100, "norm", c(0, 1), "logit"), "norm", "logit")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_silent(plot(f)); expect_silent(plot(f, which = 3))
})
