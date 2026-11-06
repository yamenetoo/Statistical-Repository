th <- list(a = 2, b = 0.7, alpha = 1.6, beta = 2.3)

test_that("density integrates to one for every built-in baseline", {
  for (bl in bkwg_baselines()$baseline) {
    B <- BKwG:::.get_baseline(bl)
    lo <- B$support[1]; up <- B$support[2]
    val <- integrate(function(x) dbkwg(x, th$a, th$b, th$alpha, th$beta,
                                       baseline = bl),
                     lo, up, rel.tol = 1e-8, subdivisions = 2000L)$value
    expect_equal(val, 1, tolerance = 1e-5, info = bl)
  }
})

test_that("cdf is the integral of the density", {
  x <- c(0.3, 1, 2.5)
  F <- pbkwg(x, th$a, th$b, th$alpha, th$beta, par = c(1.5, 2))
  I <- vapply(x, function(t) integrate(dbkwg, 0, t, a = th$a, b = th$b,
                                       alpha = th$alpha, beta = th$beta,
                                       par = c(1.5, 2))$value, 0)
  expect_equal(F, I, tolerance = 1e-7)
})

test_that("quantile inverts the cdf and tails are consistent", {
  u <- c(1e-6, 0.01, 0.3, 0.5, 0.9, 0.999)
  for (bl in c("weibull", "norm", "gamma", "lomax", "frechet", "gumbel")) {
    q <- qbkwg(u, th$a, th$b, th$alpha, th$beta, baseline = bl)
    expect_equal(pbkwg(q, th$a, th$b, th$alpha, th$beta, baseline = bl), u,
                 tolerance = 1e-6, info = bl)
  }
  x <- c(0.5, 2)
  expect_equal(pbkwg(x, 2, 1, 2, 3, lower.tail = FALSE),
               1 - pbkwg(x, 2, 1, 2, 3))
  expect_equal(pbkwg(x, 2, 1, 2, 3, log.p = TRUE), log(pbkwg(x, 2, 1, 2, 3)))
  expect_equal(qbkwg(log(0.3), 2, 1, 2, 3, log.p = TRUE), qbkwg(0.3, 2, 1, 2, 3))
  expect_equal(qbkwg(0.3, 2, 1, 2, 3, lower.tail = FALSE), qbkwg(0.7, 2, 1, 2, 3))
})

test_that("special cases reduce to known distributions", {
  x <- seq(0.1, 4, by = 0.3)
  ## G itself
  expect_equal(dbkwg(x, par = c(2, 3)), dweibull(x, 2, 3))
  expect_equal(pbkwg(x, baseline = "gamma", par = c(2, 1)), pgamma(x, 2, 1))
  ## Kumaraswamy-G
  G <- pweibull(x, 1.5, 1); g <- dweibull(x, 1.5, 1)
  expect_equal(dbkwg(x, 2, 3, 1, 1, par = c(1.5, 1)),
               2 * 3 * g * G^(2 - 1) * (1 - G^2)^(3 - 1))
  expect_equal(pbkwg(x, 2, 3, 1, 1, par = c(1.5, 1)), 1 - (1 - G^2)^3)
  ## beta-G
  expect_equal(pbkwg(x, 1, 1, 2.5, 0.7, par = c(1.5, 1)), pbeta(G, 2.5, 0.7))
  ## exponentiated Kumaraswamy-G
  expect_equal(pbkwg(x, 2, 3, 1.7, 1, par = c(1.5, 1)), (1 - (1 - G^2)^3)^1.7)
  ## beta-Kumaraswamy on (0,1) with uniform baseline
  u <- c(0.2, 0.5, 0.8)
  expect_equal(pbkwg(u, 2, 3, 1.5, 2, baseline = "unif"),
               pbeta(1 - (1 - u^2)^3, 1.5, 2))
})

test_that("closed-form density formula is reproduced", {
  x <- c(0.4, 1.1, 2.2)
  G <- plnorm(x, 0.2, 0.8); g <- dlnorm(x, 0.2, 0.8)
  K <- 1 - (1 - G^1.3)^0.6
  f <- 1.3 * 0.6 * g * G^0.3 * (1 - G^1.3)^(0.6 * 2.1 - 1) * K^(0.8 - 1) /
    beta(0.8, 2.1)
  expect_equal(dbkwg(x, 1.3, 0.6, 0.8, 2.1, baseline = "lnorm",
                     par = c(0.2, 0.8)), f)
})

test_that("survival, hazard and cumulative hazard agree", {
  x <- c(0.2, 1, 3)
  S <- sbkwg(x, 2, 1, 2, 3)
  expect_equal(hbkwg(x, 2, 1, 2, 3), dbkwg(x, 2, 1, 2, 3) / S)
  expect_equal(chbkwg(x, 2, 1, 2, 3), -log(S))
})

test_that("random generation follows the distribution", {
  set.seed(42)
  y <- rbkwg(4000, th$a, th$b, th$alpha, th$beta, par = c(1.5, 2))
  p <- ks.test(y, function(q) pbkwg(q, th$a, th$b, th$alpha, th$beta,
                                    par = c(1.5, 2)))$p.value
  expect_gt(p, 0.001)
  expect_length(rbkwg(1:7), 7)
})

test_that("analytical score matches a numerical gradient", {
  set.seed(5)
  x <- rbkwg(50, 1.7, 0.9, 1.4, 2.2, par = c(1.3, 1))
  ll <- function(v) sum(dbkwg(x, v[1], v[2], v[3], v[4], par = c(1.3, 1),
                              log = TRUE))
  v <- c(1.7, 0.9, 1.4, 2.2)
  h <- 1e-6
  num <- vapply(1:4, function(j) {
    e <- replace(numeric(4), j, h)
    (ll(v + e) - ll(v - e)) / (2 * h)
  }, 0)
  expect_equal(unname(score_bkwg(x, 1.7, 0.9, 1.4, 2.2, par = c(1.3, 1))),
               num, tolerance = 1e-5)
})

test_that("mixture representation reproduces the density", {
  w <- bkwg_mixture(2.5, 1.5, J = 400)
  expect_equal(sum(w$weight), 1, tolerance = 1e-3)
  x <- c(0.5, 1, 2)
  expect_equal(dbkwg_series(x, 2, 1, 3, 1.5, par = c(2, 1), J = 10),
               dbkwg(x, 2, 1, 3, 1.5, par = c(2, 1)))
})

test_that("input checks work", {
  expect_error(dbkwg(1, a = -1))
  expect_error(dbkwg(1, a = c(1, 2)))
  expect_error(dbkwg(1, baseline = "nope"))
  expect_error(dbkwg(1, par = c(foo = 1)))
  expect_equal(dbkwg(-1), 0)
})

test_that("user-defined baselines work, with and without a quantile function", {
  ray <- bkwg_baseline("rayleigh",
    pG = function(q, par) ifelse(q > 0, 1 - exp(-q^2 / (2 * par[["s"]]^2)), 0),
    dG = function(x, par) ifelse(x > 0, x / par[["s"]]^2 *
                                   exp(-x^2 / (2 * par[["s"]]^2)), 0),
    par = c(s = 1), lower = 0, upper = Inf, support = c(0, Inf))
  val <- integrate(dbkwg, 0, Inf, a = 2, b = 1, alpha = 1.5, beta = 0.8,
                   baseline = ray, par = c(s = 2))$value
  expect_equal(val, 1, tolerance = 1e-6)
  q <- qbkwg(0.4, 2, 1, 1.5, 0.8, baseline = ray, par = c(s = 2))
  expect_equal(pbkwg(q, 2, 1, 1.5, 0.8, baseline = ray, par = c(s = 2)), 0.4,
               tolerance = 1e-8)
})
