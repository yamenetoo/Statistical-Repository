test_that("built-in families are proper distributions", {
  for (v in c(0.3, 0.9, 2)) {
    expect_equal(integrate(dtoppleone, 0, 1, shape = v)$value, 1, tolerance = 1e-6)
    u <- c(0.1, 0.5, 0.9)
    expect_equal(ptoppleone(qtoppleone(u, v), v), u, tolerance = 1e-10)
  }
  expect_equal(integrate(dgenrayleigh, 0, Inf, shape = 1.5, rate = 2)$value, 1,
               tolerance = 1e-6)
  u <- c(0.1, 0.5, 0.9)
  expect_equal(pgenrayleigh(qgenrayleigh(u, 0.7, 3), 0.7, 3), u, tolerance = 1e-10)
})

test_that("parameters are discovered from the density", {
  W <- suppressMessages(ss_dist("weibull"))
  expect_equal(W$par_names, c("shape", "scale"))
  G <- suppressMessages(ss_dist("gamma"))
  expect_equal(G$par_names, c("shape", "rate"))
  N <- suppressMessages(ss_dist("norm"))
  expect_equal(N$constraints$mean, c(-Inf, Inf))
  expect_equal(N$lower, -Inf)
})

test_that("R matches closed forms", {
  X <- suppressMessages(ss_dist("exp", par = c(rate = 2)))
  Y <- suppressMessages(ss_dist("exp", par = c(rate = 3)))
  expect_equal(ss_R(X, Y), 3 / 5, tolerance = 1e-7)
  expect_equal(ss_R(X, Y, method = "density"), 3 / 5, tolerance = 1e-7)
  # P(cY < X): cY ~ exp(3 / c)
  expect_equal(ss_R(X, Y, c = 2), 1.5 / (2 + 1.5), tolerance = 1e-7)

  Xn <- suppressMessages(ss_dist("norm", par = c(mean = 1, sd = 1)))
  Yn <- suppressMessages(ss_dist("norm", par = c(mean = 0, sd = 2)))
  expect_equal(ss_R(Xn, Yn), pnorm(1 / sqrt(5)), tolerance = 1e-7)
})

test_that("systems match closed forms for exponentials", {
  lx <- 1; a <- 2; b <- 3
  X <- suppressMessages(ss_dist("exp", par = c(rate = lx)))
  Ya <- suppressMessages(ss_dist("exp", par = c(rate = a)))
  Yb <- suppressMessages(ss_dist("exp", par = c(rate = b)))
  exact <- 1 - lx / (lx + a) - lx / (lx + b) + lx / (lx + a + b)
  expect_equal(ss_multistress(X, list(Ya, Yb)), exact, tolerance = 1e-7)
  expect_equal(ss_multistress(X, list(Ya, Yb), components = "independent"),
               ss_R(X, Ya) * ss_R(X, Yb), tolerance = 1e-7)

  # P(Z < Y < X) with Z ~ exp(l), Y ~ exp(mu), X ~ exp(u)
  l <- 1; mu <- 2; u <- 0.5
  Z <- suppressMessages(ss_dist("exp", par = c(rate = l)))
  Y <- suppressMessages(ss_dist("exp", par = c(rate = mu)))
  U <- suppressMessages(ss_dist("exp", par = c(rate = u)))
  expect_equal(ss_R3(Z, Y, U), mu / (mu + u) - mu / (mu + l + u), tolerance = 1e-7)

  expect_equal(ss_sk(X, Ya, 1, 1), ss_R(X, Ya), tolerance = 1e-7)
  expect_equal(ss_cascade(X, Ya, n = 1)$R, ss_R(X, Ya), tolerance = 1e-7)
  P <- c(ss_R(X, Ya), ss_R(X, Ya, c = 0.5))
  expect_equal(ss_cascade(X, Ya, n = 2, c = 0.5)$R, P[1] + (1 - P[1]) * P[2],
               tolerance = 1e-7)
})

test_that("estimators recover parameters", {
  set.seed(1)
  E <- suppressMessages(ss_dist("exp"))
  x <- rexp(200, 2)
  f <- ss_fit(x, E, "mle")
  expect_equal(unname(f$estimate), 1 / mean(x), tolerance = 1e-4)
  expect_equal(unname(sqrt(f$vcov[1, 1])), (1 / mean(x)) / sqrt(200), tolerance = 1e-2)
  for (m in c("mps", "lse", "wlse", "cvm", "ad"))
    expect_lt(abs(ss_fit(x, E, m)$estimate - 2), 0.5)

  TL <- ss_dist("toppleone")
  xt <- rtoppleone(300, 0.6)
  expect_lt(abs(ss_fit(xt, TL, "mle")$estimate - 0.6), 0.1)
  GR <- ss_dist("genrayleigh")
  yg <- rgenrayleigh(300, 1.5, 2)
  expect_lt(max(abs(ss_fit(yg, GR, "mle")$estimate - c(1.5, 2))), 0.5)
})

test_that("jackknife, bootstrap and Bayes run and give sensible R", {
  skip_on_cran()
  set.seed(2)
  TL <- ss_dist("toppleone")
  GR <- ss_dist("genrayleigh")
  x <- rtoppleone(40, 0.7)
  y <- rgenrayleigh(40, 1, 2)
  true_R <- ss_R(TL, GR, c(shape = 0.7), c(shape = 1, rate = 2))
  r1 <- ss_estimate_R(x, y, TL, GR, "mle")
  expect_true(r1$lower < r1$estimate && r1$estimate < r1$upper)
  r2 <- ss_estimate_R(x, y, TL, GR, "jackknife")
  r3 <- ss_estimate_R(x, y, TL, GR, "bayes", iter = 2000, burnin = 500, ndraw = 200)
  for (r in list(r1, r2, r3)) expect_lt(abs(r$estimate - true_R), 0.2)
})
