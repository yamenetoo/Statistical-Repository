## Baseline (parent) distributions G for the BKw-G family ---------------------

## make a cdf function accept lower.tail / log.p
.wrap_p <- function(f) {
  if (all(c("lower.tail", "log.p") %in% names(formals(f)))) return(f)
  force(f)
  function(q, par, lower.tail = TRUE, log.p = FALSE) {
    p <- f(q, par)
    p <- pmin(pmax(p, 0), 1)
    if (!lower.tail) p <- 1 - p
    if (log.p) log(p) else p
  }
}

## make a density function accept log
.wrap_d <- function(f) {
  if ("log" %in% names(formals(f))) return(f)
  force(f)
  function(x, par, log = FALSE) {
    d <- f(x, par)
    if (log) log(d) else d
  }
}

## numerical inverse of a cdf when no quantile function is supplied
.numeric_q <- function(pG, support) {
  force(pG)
  force(support)
  function(p, par) {
    vapply(p, function(pp) {
      if (is.na(pp)) return(NA_real_)
      if (pp <= 0) return(support[1L])
      if (pp >= 1) return(support[2L])
      f <- function(x) pG(x, par) - pp
      lo <- if (is.finite(support[1L])) support[1L] else -1
      hi <- if (is.finite(support[2L])) support[2L] else 1
      it <- 0L
      while (f(lo) > 0 && it < 2000L) {
        lo <- if (lo < 0) 2 * lo else lo - 1
        it <- it + 1L
      }
      while (f(hi) < 0 && it < 4000L) {
        hi <- if (hi > 0) 2 * hi else hi + 1
        it <- it + 1L
      }
      stats::uniroot(f, c(lo, hi), tol = 1e-12)$root
    }, numeric(1))
  }
}

bkwg_baseline <- function(name, pG, dG, qG = NULL, par = numeric(0),
                          lower = NULL, upper = NULL,
                          support = c(-Inf, Inf), start = NULL) {
  if (!is.character(name) || length(name) != 1L)
    stop("'name' must be a single character string")
  if (!is.function(pG) || !is.function(dG))
    stop("'pG' and 'dG' must be functions of (x, par)")
  par <- unlist(par)
  if (length(par)) {
    if (!is.numeric(par) || is.null(names(par)) || any(names(par) == ""))
      stop("'par' must be a named numeric vector of default parameter values")
  } else {
    par <- numeric(0)
  }
  np <- length(par)
  if (is.null(lower)) lower <- rep(-Inf, np)
  if (is.null(upper)) upper <- rep(Inf, np)
  lower <- rep_len(as.numeric(lower), np)
  upper <- rep_len(as.numeric(upper), np)
  if (np) names(lower) <- names(upper) <- names(par)
  if (length(support) != 2L || support[1L] >= support[2L])
    stop("'support' must be c(lower, upper) with lower < upper")
  pG <- .wrap_p(pG)
  dG <- .wrap_d(dG)
  if (is.null(qG)) {
    qG <- .numeric_q(pG, support)
  } else if (!is.function(qG)) {
    stop("'qG' must be NULL or a function of (p, par)")
  }
  if (is.null(start)) {
    start <- function(x) par
  } else if (!is.function(start)) {
    sv <- unlist(start)
    start <- function(x) sv
  }
  structure(list(name = name, pG = pG, dG = dG, qG = qG, par = par,
                 lower = lower, upper = upper, support = support,
                 start = start),
            class = "bkwg_baseline")
}

print.bkwg_baseline <- function(x, ...) {
  cat("Baseline distribution G: ", x$name, "\n", sep = "")
  if (length(x$par)) {
    cat("Parameters (default values):\n")
    print(x$par)
  } else {
    cat("No parameters\n")
  }
  cat("Support: (", x$support[1L], ", ", x$support[2L], ")\n", sep = "")
  invisible(x)
}

## Built-in baselines ---------------------------------------------------------

.bl_exp <- function() bkwg_baseline(
  "exp",
  pG = function(q, par, lower.tail = TRUE, log.p = FALSE)
    stats::pexp(q, rate = par[["rate"]], lower.tail = lower.tail,
                log.p = log.p),
  dG = function(x, par, log = FALSE)
    stats::dexp(x, rate = par[["rate"]], log = log),
  qG = function(p, par) stats::qexp(p, rate = par[["rate"]]),
  par = c(rate = 1), lower = 0, upper = Inf, support = c(0, Inf),
  start = function(x) c(rate = 1 / mean(x)))

.bl_weibull <- function() bkwg_baseline(
  "weibull",
  pG = function(q, par, lower.tail = TRUE, log.p = FALSE)
    stats::pweibull(q, shape = par[["shape"]], scale = par[["scale"]],
                    lower.tail = lower.tail, log.p = log.p),
  dG = function(x, par, log = FALSE)
    stats::dweibull(x, shape = par[["shape"]], scale = par[["scale"]],
                    log = log),
  qG = function(p, par)
    stats::qweibull(p, shape = par[["shape"]], scale = par[["scale"]]),
  par = c(shape = 1, scale = 1), lower = c(0, 0), upper = c(Inf, Inf),
  support = c(0, Inf),
  start = function(x) {
    lx <- log(x[x > 0])
    k <- 1.2825 / stats::sd(lx)
    if (!is.finite(k) || k <= 0) k <- 1
    c(shape = k, scale = exp(mean(lx) + 0.5772 / k))
  })

.bl_norm <- function() bkwg_baseline(
  "norm",
  pG = function(q, par, lower.tail = TRUE, log.p = FALSE)
    stats::pnorm(q, mean = par[["mean"]], sd = par[["sd"]],
                 lower.tail = lower.tail, log.p = log.p),
  dG = function(x, par, log = FALSE)
    stats::dnorm(x, mean = par[["mean"]], sd = par[["sd"]], log = log),
  qG = function(p, par) stats::qnorm(p, mean = par[["mean"]], sd = par[["sd"]]),
  par = c(mean = 0, sd = 1), lower = c(-Inf, 0), upper = c(Inf, Inf),
  support = c(-Inf, Inf),
  start = function(x) c(mean = mean(x), sd = stats::sd(x)))

.bl_lnorm <- function() bkwg_baseline(
  "lnorm",
  pG = function(q, par, lower.tail = TRUE, log.p = FALSE)
    stats::plnorm(q, meanlog = par[["meanlog"]], sdlog = par[["sdlog"]],
                  lower.tail = lower.tail, log.p = log.p),
  dG = function(x, par, log = FALSE)
    stats::dlnorm(x, meanlog = par[["meanlog"]], sdlog = par[["sdlog"]],
                  log = log),
  qG = function(p, par)
    stats::qlnorm(p, meanlog = par[["meanlog"]], sdlog = par[["sdlog"]]),
  par = c(meanlog = 0, sdlog = 1), lower = c(-Inf, 0), upper = c(Inf, Inf),
  support = c(0, Inf),
  start = function(x) {
    lx <- log(x[x > 0])
    c(meanlog = mean(lx), sdlog = stats::sd(lx))
  })

.bl_gamma <- function() bkwg_baseline(
  "gamma",
  pG = function(q, par, lower.tail = TRUE, log.p = FALSE)
    stats::pgamma(q, shape = par[["shape"]], rate = par[["rate"]],
                  lower.tail = lower.tail, log.p = log.p),
  dG = function(x, par, log = FALSE)
    stats::dgamma(x, shape = par[["shape"]], rate = par[["rate"]], log = log),
  qG = function(p, par)
    stats::qgamma(p, shape = par[["shape"]], rate = par[["rate"]]),
  par = c(shape = 1, rate = 1), lower = c(0, 0), upper = c(Inf, Inf),
  support = c(0, Inf),
  start = function(x) {
    m <- mean(x)
    v <- stats::var(x)
    c(shape = m^2 / v, rate = m / v)
  })

.bl_logis <- function() bkwg_baseline(
  "logis",
  pG = function(q, par, lower.tail = TRUE, log.p = FALSE)
    stats::plogis(q, location = par[["location"]], scale = par[["scale"]],
                  lower.tail = lower.tail, log.p = log.p),
  dG = function(x, par, log = FALSE)
    stats::dlogis(x, location = par[["location"]], scale = par[["scale"]],
                  log = log),
  qG = function(p, par)
    stats::qlogis(p, location = par[["location"]], scale = par[["scale"]]),
  par = c(location = 0, scale = 1), lower = c(-Inf, 0), upper = c(Inf, Inf),
  support = c(-Inf, Inf),
  start = function(x)
    c(location = stats::median(x), scale = stats::sd(x) * sqrt(3) / pi))

.bl_gumbel <- function() bkwg_baseline(
  "gumbel",
  pG = function(q, par, lower.tail = TRUE, log.p = FALSE) {
    z <- (q - par[["location"]]) / par[["scale"]]
    .from_logF(-exp(-z), lower.tail, log.p)
  },
  dG = function(x, par, log = FALSE) {
    z <- (x - par[["location"]]) / par[["scale"]]
    ld <- -log(par[["scale"]]) - z - exp(-z)
    if (log) ld else exp(ld)
  },
  qG = function(p, par) par[["location"]] - par[["scale"]] * log(-log(p)),
  par = c(location = 0, scale = 1), lower = c(-Inf, 0), upper = c(Inf, Inf),
  support = c(-Inf, Inf),
  start = function(x) {
    s <- stats::sd(x) * sqrt(6) / pi
    c(location = mean(x) - 0.5772157 * s, scale = s)
  })

.bl_lomax <- function() bkwg_baseline(
  "lomax",
  pG = function(q, par, lower.tail = TRUE, log.p = FALSE) {
    q <- pmax(q, 0)
    .from_logS(-par[["shape"]] * log1p(q / par[["scale"]]), lower.tail, log.p)
  },
  dG = function(x, par, log = FALSE) {
    al <- par[["shape"]]
    la <- par[["scale"]]
    ld <- log(al) - log(la) - (al + 1) * log1p(pmax(x, 0) / la)
    ld[!is.na(x) & x < 0] <- -Inf
    if (log) ld else exp(ld)
  },
  qG = function(p, par)
    par[["scale"]] * ((1 - p)^(-1 / par[["shape"]]) - 1),
  par = c(shape = 2, scale = 1), lower = c(0, 0), upper = c(Inf, Inf),
  support = c(0, Inf),
  start = function(x) c(shape = 3, scale = 2 * mean(x)))

.bl_llogis <- function() bkwg_baseline(
  "llogis",
  pG = function(q, par, lower.tail = TRUE, log.p = FALSE) {
    lz <- par[["shape"]] * (log(pmax(q, 0)) - log(par[["scale"]]))
    .from_logS(-.log1pexp(lz), lower.tail, log.p)
  },
  dG = function(x, par, log = FALSE) {
    k <- par[["shape"]]
    s <- par[["scale"]]
    lr <- log(pmax(x, 0)) - log(s)
    ld <- log(k) - log(s) + .mul(k - 1, lr) - 2 * .log1pexp(k * lr)
    ld[!is.na(x) & x <= 0] <- -Inf
    if (log) ld else exp(ld)
  },
  qG = function(p, par) par[["scale"]] * (p / (1 - p))^(1 / par[["shape"]]),
  par = c(shape = 1, scale = 1), lower = c(0, 0), upper = c(Inf, Inf),
  support = c(0, Inf),
  start = function(x) {
    lx <- log(x[x > 0])
    k <- pi / (sqrt(3) * stats::sd(lx))
    if (!is.finite(k) || k <= 0) k <- 1
    c(shape = k, scale = exp(stats::median(lx)))
  })

.bl_frechet <- function() bkwg_baseline(
  "frechet",
  pG = function(q, par, lower.tail = TRUE, log.p = FALSE) {
    z <- exp(-par[["shape"]] * (log(pmax(q, 0)) - log(par[["scale"]])))
    .from_logF(-z, lower.tail, log.p)
  },
  dG = function(x, par, log = FALSE) {
    k <- par[["shape"]]
    s <- par[["scale"]]
    lr <- log(pmax(x, 0)) - log(s)
    ld <- log(k) - log(s) - (k + 1) * lr - exp(-k * lr)
    ld[!is.na(x) & x <= 0] <- -Inf
    if (log) ld else exp(ld)
  },
  qG = function(p, par) par[["scale"]] * (-log(p))^(-1 / par[["shape"]]),
  par = c(shape = 1, scale = 1), lower = c(0, 0), upper = c(Inf, Inf),
  support = c(0, Inf),
  start = function(x) {
    lx <- log(x[x > 0])
    k <- pi / (sqrt(6) * stats::sd(lx))
    if (!is.finite(k) || k <= 0) k <- 1
    c(shape = k, scale = exp(mean(lx) - 0.5772157 / k))
  })

.bl_burrxii <- function() bkwg_baseline(
  "burrxii",
  pG = function(q, par, lower.tail = TRUE, log.p = FALSE) {
    lr <- log(pmax(q, 0)) - log(par[["scale"]])
    .from_logS(-par[["k"]] * .log1pexp(par[["c"]] * lr), lower.tail, log.p)
  },
  dG = function(x, par, log = FALSE) {
    cc <- par[["c"]]
    k <- par[["k"]]
    s <- par[["scale"]]
    lr <- log(pmax(x, 0)) - log(s)
    ld <- log(cc) + log(k) - log(s) + .mul(cc - 1, lr) -
      (k + 1) * .log1pexp(cc * lr)
    ld[!is.na(x) & x <= 0] <- -Inf
    if (log) ld else exp(ld)
  },
  qG = function(p, par)
    par[["scale"]] * ((1 - p)^(-1 / par[["k"]]) - 1)^(1 / par[["c"]]),
  par = c(c = 1, k = 1, scale = 1), lower = c(0, 0, 0),
  upper = c(Inf, Inf, Inf), support = c(0, Inf),
  start = function(x) {
    lx <- log(x[x > 0])
    cc <- pi / (sqrt(3) * stats::sd(lx))
    if (!is.finite(cc) || cc <= 0) cc <- 1
    c(c = cc, k = 1, scale = exp(stats::median(lx)))
  })

.bl_unif <- function() bkwg_baseline(
  "unif",
  pG = function(q, par, lower.tail = TRUE, log.p = FALSE) {
    p <- pmin(pmax(q, 0), 1)
    if (!lower.tail) p <- 1 - p
    if (log.p) log(p) else p
  },
  dG = function(x, par, log = FALSE) {
    d <- ifelse(x >= 0 & x <= 1, 1, 0)
    if (log) log(d) else d
  },
  qG = function(p, par) p,
  par = numeric(0), support = c(0, 1))

.baseline_names <- c("exp", "weibull", "norm", "lnorm", "gamma", "logis",
                     "gumbel", "lomax", "llogis", "frechet", "burrxii",
                     "unif")

.builtin_baseline <- function(name) {
  switch(name,
         exp = .bl_exp(), weibull = .bl_weibull(), norm = .bl_norm(),
         lnorm = .bl_lnorm(), gamma = .bl_gamma(), logis = .bl_logis(),
         gumbel = .bl_gumbel(), lomax = .bl_lomax(), llogis = .bl_llogis(),
         frechet = .bl_frechet(), burrxii = .bl_burrxii(), unif = .bl_unif())
}

.get_baseline <- function(baseline) {
  if (inherits(baseline, "bkwg_baseline")) return(baseline)
  if (is.character(baseline) && length(baseline) == 1L) {
    nm <- tolower(baseline)
    if (!nm %in% .baseline_names)
      stop("unknown baseline '", baseline, "'; available: ",
           paste(.baseline_names, collapse = ", "),
           ", or build one with bkwg_baseline()", call. = FALSE)
    return(.builtin_baseline(nm))
  }
  stop("'baseline' must be a name of a built-in baseline or an object ",
       "created by bkwg_baseline()", call. = FALSE)
}

bkwg_baselines <- function() {
  out <- lapply(.baseline_names, function(nm) {
    B <- .builtin_baseline(nm)
    data.frame(baseline = nm,
               parameters = if (length(B$par)) paste(names(B$par), collapse = ", ") else "(none)",
               support = sprintf("(%s, %s)", B$support[1L], B$support[2L]),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}
