## The beta Kumaraswamy-G (BKw-G) family ----------------------------------------
##
## Kw-G cdf:   K(x) = 1 - [1 - G(x)^a]^b
## BKw-G cdf:  F(x) = I_{K(x)}(alpha, beta)
## BKw-G pdf:  f(x) = a b g(x) G^{a-1} (1 - G^a)^{b beta - 1} K^{alpha - 1} / B(alpha, beta)
##
## (Handique, Chakraborty and Ali, 2017; alpha = m and beta = n in their notation)

## log G, log(1 - G^a), log(1 - K) and log K, all on the log scale
.parts <- function(x, a, b, B, par) {
  logG <- B$pG(x, par, log.p = TRUE)
  log1mGa <- .log1mexp(-a * logG)
  logSK <- b * log1mGa
  logK <- .log1mexp(-logSK)
  list(logG = logG, log1mGa = log1mGa, logSK = logSK, logK = logK)
}

.ldens <- function(x, a, b, alpha, beta, B, par, check = FALSE) {
  P <- .parts(x, a, b, B, par)
  lg <- B$dG(x, par, log = TRUE)
  t2 <- .mul(a - 1, P$logG)
  t3 <- .mul(b * beta - 1, P$log1mGa)
  t4 <- .mul(alpha - 1, P$logK)
  ld <- log(a) + log(b) + lg + t2 + t3 + t4 - lbeta(alpha, beta)
  if (check) {
    ## for extreme parameter values the terms can be huge (say 1e26) and
    ## cancel, leaving a log-density that is pure rounding error; the
    ## optimiser must not be rewarded for such points
    err <- 8 * .Machine$double.eps *
      (abs(lg) + abs(t2) + abs(t3) + abs(t4) + abs(log(a)) + abs(log(b)) +
         abs(lbeta(alpha, beta)))
    ld[!is.na(err) & err > 1e-3] <- -Inf
  }
  ld[!is.na(lg) & lg == -Inf] <- -Inf
  ## G(x) rounded to 1 inside the support (e.g. a user pG without
  ## lower.tail/log.p): (1 - G^a)^(b*beta - 1) is then Inf * g(x), but the
  ## true density there is negligible, so return 0 instead of Inf
  ld[!is.na(ld) & ld == Inf & is.finite(lg) & P$log1mGa == -Inf &
       x < B$support[2L]] <- -Inf
  ld[is.nan(ld)] <- -Inf
  ld
}

.pbkwg_raw <- function(q, a, b, alpha, beta, B, par, lower.tail = TRUE,
                       log.p = FALSE) {
  P <- .parts(q, a, b, B, par)
  K <- exp(P$logK)
  SK <- exp(P$logSK)
  out <- rep(NA_real_, length(q))
  i <- !is.na(K) & K <= 0.5
  j <- !is.na(K) & !i
  out[i] <- stats::pbeta(K[i], alpha, beta, lower.tail = lower.tail,
                         log.p = log.p)
  out[j] <- stats::pbeta(SK[j], beta, alpha, lower.tail = !lower.tail,
                         log.p = log.p)
  out
}

.qbkwg_raw <- function(p, a, b, alpha, beta, B, par, lower.tail = TRUE,
                       log.p = FALSE) {
  ## 1 - Z where Z is the p-quantile of Beta(alpha, beta)
  SK <- stats::qbeta(p, beta, alpha, lower.tail = !lower.tail, log.p = log.p)
  logGa <- .log1mexp(-log(SK) / b)      # log(1 - SK^(1/b)) = a log G
  B$qG(exp(logGa / a), par)
}

dbkwg <- function(x, a = 1, b = 1, alpha = 1, beta = 1,
                  baseline = "weibull", par = NULL, log = FALSE) {
  s <- .prep(a, b, alpha, beta, baseline, par)
  ld <- .ldens(x, s$a, s$b, s$alpha, s$beta, s$B, s$par)
  if (log) ld else exp(ld)
}

pbkwg <- function(q, a = 1, b = 1, alpha = 1, beta = 1,
                  baseline = "weibull", par = NULL,
                  lower.tail = TRUE, log.p = FALSE) {
  s <- .prep(a, b, alpha, beta, baseline, par)
  .pbkwg_raw(q, s$a, s$b, s$alpha, s$beta, s$B, s$par,
             lower.tail = lower.tail, log.p = log.p)
}

qbkwg <- function(p, a = 1, b = 1, alpha = 1, beta = 1,
                  baseline = "weibull", par = NULL,
                  lower.tail = TRUE, log.p = FALSE) {
  s <- .prep(a, b, alpha, beta, baseline, par)
  .qbkwg_raw(p, s$a, s$b, s$alpha, s$beta, s$B, s$par,
             lower.tail = lower.tail, log.p = log.p)
}

rbkwg <- function(n, a = 1, b = 1, alpha = 1, beta = 1,
                  baseline = "weibull", par = NULL) {
  if (length(n) > 1L) n <- length(n)
  s <- .prep(a, b, alpha, beta, baseline, par)
  SK <- stats::rbeta(n, s$beta, s$alpha)     # 1 - Z,  Z ~ Beta(alpha, beta)
  G <- exp(.log1mexp(-log(SK) / s$b) / s$a)
  s$B$qG(G, s$par)
}

sbkwg <- function(x, a = 1, b = 1, alpha = 1, beta = 1,
                  baseline = "weibull", par = NULL, log = FALSE) {
  s <- .prep(a, b, alpha, beta, baseline, par)
  .pbkwg_raw(x, s$a, s$b, s$alpha, s$beta, s$B, s$par,
             lower.tail = FALSE, log.p = log)
}

hbkwg <- function(x, a = 1, b = 1, alpha = 1, beta = 1,
                  baseline = "weibull", par = NULL, log = FALSE) {
  s <- .prep(a, b, alpha, beta, baseline, par)
  ld <- .ldens(x, s$a, s$b, s$alpha, s$beta, s$B, s$par)
  lS <- .pbkwg_raw(x, s$a, s$b, s$alpha, s$beta, s$B, s$par,
                   lower.tail = FALSE, log.p = TRUE)
  lh <- ld - lS
  if (log) lh else exp(lh)
}

chbkwg <- function(x, a = 1, b = 1, alpha = 1, beta = 1,
                   baseline = "weibull", par = NULL) {
  s <- .prep(a, b, alpha, beta, baseline, par)
  -.pbkwg_raw(x, s$a, s$b, s$alpha, s$beta, s$B, s$par,
              lower.tail = FALSE, log.p = TRUE)
}

## Analytical score of (a, b, alpha, beta) for fixed baseline parameters
score_bkwg <- function(x, a = 1, b = 1, alpha = 1, beta = 1,
                       baseline = "weibull", par = NULL, sum = TRUE) {
  s <- .prep(a, b, alpha, beta, baseline, par)
  P <- .parts(x, s$a, s$b, s$B, s$par)
  Ta <- exp(s$a * P$logG)                 # G^a
  oneMT <- exp(P$log1mGa)                 # 1 - G^a
  U <- exp(P$logSK)                       # 1 - K
  K <- exp(P$logK)
  dKda <- s$b * exp((s$b - 1) * P$log1mGa) * Ta * P$logG
  s_a <- 1 / s$a + P$logG - (s$b * s$beta - 1) * Ta * P$logG / oneMT +
    (s$alpha - 1) * dKda / K
  s_b <- 1 / s$b + s$beta * P$log1mGa - (s$alpha - 1) * U * P$log1mGa / K
  s_alpha <- P$logK - digamma(s$alpha) + digamma(s$alpha + s$beta)
  s_beta <- s$b * P$log1mGa - digamma(s$beta) + digamma(s$alpha + s$beta)
  S <- cbind(a = s_a, b = s_b, alpha = s_alpha, beta = s_beta)
  if (sum) colSums(S) else S
}
