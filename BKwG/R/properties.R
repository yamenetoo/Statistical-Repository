## Structural properties of the BKw-G family -----------------------------------

.Qfun <- function(s) {
  function(u) .qbkwg_raw(u, s$a, s$b, s$alpha, s$beta, s$B, s$par)
}

## r-th raw moment  E[X^r] = int_0^1 Q(u)^r du
mbkwg <- function(r = 1, a = 1, b = 1, alpha = 1, beta = 1,
                  baseline = "weibull", par = NULL, ...) {
  s <- .prep(a, b, alpha, beta, baseline, par)
  Q <- .Qfun(s)
  vapply(r, function(rr) .integrate(function(u) Q(u)^rr, 0, 1, ...),
         numeric(1))
}

## moment generating function  M(t) = int_0^1 exp(t Q(u)) du
mgf_bkwg <- function(t, a = 1, b = 1, alpha = 1, beta = 1,
                     baseline = "weibull", par = NULL, ...) {
  s <- .prep(a, b, alpha, beta, baseline, par)
  Q <- .Qfun(s)
  vapply(t, function(tt) .integrate(function(u) exp(tt * Q(u)), 0, 1, ...),
         numeric(1))
}

## mean, variance, moment and quantile based skewness / kurtosis
bkwg_stats <- function(a = 1, b = 1, alpha = 1, beta = 1,
                       baseline = "weibull", par = NULL, ...) {
  s <- .prep(a, b, alpha, beta, baseline, par)
  Q <- .Qfun(s)
  m <- vapply(1:4, function(rr) .integrate(function(u) Q(u)^rr, 0, 1, ...),
              numeric(1))
  v <- m[2] - m[1]^2
  mu3 <- m[3] - 3 * m[1] * m[2] + 2 * m[1]^3
  mu4 <- m[4] - 4 * m[1] * m[3] + 6 * m[1]^2 * m[2] - 3 * m[1]^4
  q <- Q(c(1, 2, 3, 4, 5, 6, 7) / 8)
  c(mean = m[1], variance = v, sd = sqrt(max(v, 0)),
    skewness = mu3 / v^1.5, kurtosis = mu4 / v^2,
    median = q[4],
    bowley = (q[6] + q[2] - 2 * q[4]) / (q[6] - q[2]),
    moors = (q[7] - q[5] + q[3] - q[1]) / (q[6] - q[2]))
}

## Shannon and Renyi entropies
entropy_bkwg <- function(a = 1, b = 1, alpha = 1, beta = 1,
                         baseline = "weibull", par = NULL,
                         type = c("shannon", "renyi"), delta = 2, ...) {
  type <- match.arg(type)
  s <- .prep(a, b, alpha, beta, baseline, par)
  Q <- .Qfun(s)
  lf <- function(u) .ldens(Q(u), s$a, s$b, s$alpha, s$beta, s$B, s$par)
  if (type == "shannon") {
    return(-.integrate(lf, 0, 1, ...))
  }
  if (length(delta) != 1L || delta <= 0 || delta == 1)
    stop("'delta' must be a positive number different from 1")
  I <- .integrate(function(u) exp((delta - 1) * lf(u)), 0, 1, ...)
  log(I) / (1 - delta)
}

## mean residual life  m(t) = E[X - t | X > t]
mrl_bkwg <- function(t, a = 1, b = 1, alpha = 1, beta = 1,
                     baseline = "weibull", par = NULL, ...) {
  s <- .prep(a, b, alpha, beta, baseline, par)
  Q <- .Qfun(s)
  vapply(t, function(tt) {
    Ft <- .pbkwg_raw(tt, s$a, s$b, s$alpha, s$beta, s$B, s$par)
    if (!is.finite(Ft) || Ft >= 1) return(NA_real_)
    .integrate(Q, Ft, 1, ...) / (1 - Ft) - tt
  }, numeric(1))
}

## density and cdf of the i-th order statistic in a sample of size N
.check_os <- function(i, N) {
  if (length(i) != 1L || length(N) != 1L || i < 1 || N < i ||
      i != round(i) || N != round(N))
    stop("'i' and 'N' must be integers with 1 <= i <= N")
}

dos_bkwg <- function(x, i, N, a = 1, b = 1, alpha = 1, beta = 1,
                     baseline = "weibull", par = NULL, log = FALSE) {
  .check_os(i, N)
  s <- .prep(a, b, alpha, beta, baseline, par)
  ld <- .ldens(x, s$a, s$b, s$alpha, s$beta, s$B, s$par)
  lF <- .pbkwg_raw(x, s$a, s$b, s$alpha, s$beta, s$B, s$par, log.p = TRUE)
  lS <- .pbkwg_raw(x, s$a, s$b, s$alpha, s$beta, s$B, s$par,
                   lower.tail = FALSE, log.p = TRUE)
  out <- lgamma(N + 1) - lgamma(i) - lgamma(N - i + 1) + ld +
    .mul(i - 1, lF) + .mul(N - i, lS)
  out[ld == -Inf] <- -Inf
  if (log) out else exp(out)
}

pos_bkwg <- function(q, i, N, a = 1, b = 1, alpha = 1, beta = 1,
                     baseline = "weibull", par = NULL) {
  .check_os(i, N)
  s <- .prep(a, b, alpha, beta, baseline, par)
  Fq <- .pbkwg_raw(q, s$a, s$b, s$alpha, s$beta, s$B, s$par)
  stats::pbeta(Fq, i, N - i + 1)
}

## Linear representation: f(x) = sum_j w_j * KwG(x; a, b (beta + j))
bkwg_mixture <- function(alpha, beta, J = 100) {
  if (!.shape_ok(1, 1, alpha, beta))
    stop("'alpha' and 'beta' must be positive")
  j <- 0:J
  w <- (-1)^j * choose(alpha - 1, j) * exp(-lbeta(alpha, beta)) / (beta + j)
  data.frame(j = j, weight = w, b_multiplier = beta + j)
}

dbkwg_series <- function(x, a = 1, b = 1, alpha = 1, beta = 1,
                         baseline = "weibull", par = NULL, J = 100) {
  s <- .prep(a, b, alpha, beta, baseline, par)
  mix <- bkwg_mixture(s$alpha, s$beta, J)
  out <- numeric(length(x))
  for (k in seq_len(nrow(mix))) {
    if (mix$weight[k] == 0) next
    out <- out + mix$weight[k] *
      exp(.ldens(x, s$a, s$b * mix$b_multiplier[k], 1, 1, s$B, s$par))
  }
  out
}

## Quick plots of pdf, cdf, survival or hazard
plot_bkwg <- function(a = 1, b = 1, alpha = 1, beta = 1,
                      baseline = "weibull", par = NULL,
                      what = c("pdf", "cdf", "survival", "hazard"),
                      from = NULL, to = NULL, add = FALSE, n = 401, ...) {
  what <- match.arg(what)
  s <- .prep(a, b, alpha, beta, baseline, par)
  if (is.null(from))
    from <- .qbkwg_raw(0.001, s$a, s$b, s$alpha, s$beta, s$B, s$par)
  if (is.null(to))
    to <- .qbkwg_raw(0.999, s$a, s$b, s$alpha, s$beta, s$B, s$par)
  xx <- seq(from, to, length.out = n)
  yy <- switch(what,
    pdf = exp(.ldens(xx, s$a, s$b, s$alpha, s$beta, s$B, s$par)),
    cdf = .pbkwg_raw(xx, s$a, s$b, s$alpha, s$beta, s$B, s$par),
    survival = .pbkwg_raw(xx, s$a, s$b, s$alpha, s$beta, s$B, s$par,
                          lower.tail = FALSE),
    hazard = exp(.ldens(xx, s$a, s$b, s$alpha, s$beta, s$B, s$par) -
                   .pbkwg_raw(xx, s$a, s$b, s$alpha, s$beta, s$B, s$par,
                              lower.tail = FALSE, log.p = TRUE)))
  if (add) {
    graphics::lines(xx, yy, ...)
  } else {
    args <- list(...)
    if (is.null(args$xlab)) args$xlab <- "x"
    if (is.null(args$ylab))
      args$ylab <- switch(what, pdf = "f(x)", cdf = "F(x)",
                          survival = "S(x)", hazard = "h(x)")
    if (is.null(args$type)) args$type <- "l"
    do.call(graphics::plot, c(list(x = xx, y = yy), args))
  }
  invisible(data.frame(x = xx, y = yy))
}
