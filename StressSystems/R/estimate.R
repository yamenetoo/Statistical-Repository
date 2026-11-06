.point_methods <- c("mle", "mps", "lse", "wlse", "cvm", "ad")

#' Estimate the parameters of a distribution
#'
#' Works for any \code{\link{ss_dist}}. Optimisation is done on an
#' unconstrained scale (log for positive parameters, logit for bounded ones),
#' starting from a screened set of random starting points, so no starting
#' values are needed.
#'
#' Point methods:
#' \describe{
#'   \item{\code{"mle"}}{maximum likelihood (with observed-information covariance).}
#'   \item{\code{"mps"}}{maximum product of spacings.}
#'   \item{\code{"lse"}, \code{"wlse"}}{(weighted) least squares on the CDF.}
#'   \item{\code{"cvm"}}{Cramer-von Mises minimum distance.}
#'   \item{\code{"ad"}}{Anderson-Darling minimum distance.}
#' }
#' Resampling and Bayesian methods: \code{"jackknife"}, \code{"bootstrap"}
#' (bias-corrected, built on \code{base}) and \code{"bayes"}; see
#' \code{\link{ss_jackknife}}, \code{\link{ss_bootstrap}}, \code{\link{ss_bayes}}.
#'
#' @param x Data.
#' @param dist An \code{ss_dist}.
#' @param method Estimation method.
#' @param base Point method used inside the jackknife and bootstrap.
#' @param start Optional starting values (named or in parameter order).
#' @param ... Passed to \code{ss_bootstrap} or \code{ss_bayes}.
#' @return An object of class \code{"ss_fit"}.
#' @export
ss_fit <- function(x, dist, method = c(.point_methods, "jackknife", "bootstrap", "bayes"),
                   base = "mle", start = NULL, ...) {
  method <- match.arg(method)
  x <- .check_data(x, dist)
  switch(method,
         jackknife = ss_jackknife(x, dist, base = base, start = start),
         bootstrap = ss_bootstrap(x, dist, base = base, start = start, ...),
         bayes = ss_bayes(x, dist, start = start, ...),
         .fit_point(x, dist, method, start))
}

#' Jackknife bias-corrected estimator
#'
#' Computes \eqn{\tilde\theta = n\hat\theta - (n-1)\bar\theta_{(\cdot)}}, where
#' \eqn{\bar\theta_{(\cdot)}} is the mean of the leave-one-out estimates, and
#' the jackknife standard error. If the corrected value leaves the parameter
#' space, the mean of the leave-one-out estimates is returned with a warning.
#'
#' @inheritParams ss_fit
#' @export
ss_jackknife <- function(x, dist, base = "mle", start = NULL) {
  base <- match.arg(base, .point_methods)
  x <- .check_data(x, dist)
  full <- .fit_point(x, dist, base, start)
  n <- length(x)
  loo <- do.call(rbind, lapply(seq_len(n), function(i)
    tryCatch(.fit_point(x[-i], dist, base, start = full$estimate, screen = FALSE)$estimate,
             error = function(e) setNames(rep(NA_real_, length(dist$par_names)),
                                          dist$par_names))))
  m <- colMeans(loo, na.rm = TRUE)
  est <- n * full$estimate - (n - 1) * m
  if (!.inside(est, dist)) {
    warning("jackknife correction left the parameter space; using the leave-one-out mean")
    est <- m
  }
  se <- sqrt((n - 1) / n * colSums(sweep(loo, 2, m)^2, na.rm = TRUE))
  structure(list(estimate = est, method = "jackknife", base = base, se = se,
                 full = full, replicates = loo, dist = dist, n = n),
            class = "ss_fit")
}

#' Bootstrap bias-corrected estimator
#'
#' Computes \eqn{2\hat\theta - \bar\theta^{*}} from \code{B} bootstrap
#' re-estimates, together with bootstrap standard errors.
#'
#' @inheritParams ss_fit
#' @param B Number of bootstrap samples.
#' @param type \code{"parametric"} or \code{"nonparametric"}.
#' @param ... Ignored.
#' @export
ss_bootstrap <- function(x, dist, base = "mle", start = NULL, B = 200,
                         type = c("parametric", "nonparametric"), ...) {
  base <- match.arg(base, .point_methods)
  type <- match.arg(type)
  x <- .check_data(x, dist)
  full <- .fit_point(x, dist, base, start)
  n <- length(x)
  reps <- do.call(rbind, lapply(seq_len(B), function(b) tryCatch({
    xb <- if (type == "parametric") .rand(dist, n, full$estimate) else sample(x, n, TRUE)
    .fit_point(xb, dist, base, start = full$estimate, screen = FALSE)$estimate
  }, error = function(e) setNames(rep(NA_real_, length(dist$par_names)), dist$par_names))))
  m <- colMeans(reps, na.rm = TRUE)
  est <- 2 * full$estimate - m
  if (!.inside(est, dist)) {
    warning("bootstrap correction left the parameter space; using the uncorrected estimate")
    est <- full$estimate
  }
  se <- apply(reps, 2, sd, na.rm = TRUE)
  structure(list(estimate = est, method = "bootstrap", base = base, type = type,
                 se = se, full = full, replicates = reps, dist = dist, n = n),
            class = "ss_fit")
}

#' Bayesian estimation by Metropolis-Hastings
#'
#' Random-walk Metropolis-Hastings on the unconstrained scale with a
#' multivariate normal proposal whose covariance starts from the inverse
#' observed information of the MLE and whose scale is tuned during burn-in.
#' The Jacobian of the transformation is included, so the prior is on the
#' natural parameter scale.
#'
#' @inheritParams ss_fit
#' @param prior Function of the named parameter vector returning the log prior
#'   density. \code{NULL} gives a flat prior on the natural scale. For example,
#'   independent Gamma(1, 1) priors:
#'   \code{function(th) sum(dgamma(th, 1, 1, log = TRUE))}.
#' @param iter Number of draws kept after burn-in (before thinning).
#' @param burnin Number of burn-in iterations.
#' @param thin Thinning interval.
#' @param proposal Optional proposal covariance on the unconstrained scale.
#' @param ... Ignored.
#' @return An \code{ss_fit} with the posterior mean as \code{estimate}, the
#'   posterior median, the draws and the acceptance rate.
#' @export
ss_bayes <- function(x, dist, prior = NULL, iter = 10000, burnin = 2000, thin = 1,
                     start = NULL, proposal = NULL, ...) {
  x <- .check_data(x, dist)
  init <- tryCatch(.fit_point(x, dist, "mle", start), error = function(e) NULL)
  eta <- if (!is.null(init)) init$eta else
    .to_free(if (is.null(start)) dist$default else start, dist)
  k <- length(eta)
  S <- proposal
  if (is.null(S)) {
    S <- init$vcov_free
    if (is.null(S) || any(!is.finite(S))) S <- diag(0.1, k)
  }
  S <- (2.38^2 / k) * as.matrix(S)
  L <- tryCatch(t(chol(S)), error = function(e) diag(sqrt(pmax(diag(S), 1e-6)), k))

  lp <- function(e) {
    th <- .from_free(e, dist)
    ll <- tryCatch(suppressWarnings(sum(.dens(dist, x, th, log = TRUE))),
                   error = function(err) -Inf)
    pr <- if (is.null(prior)) 0 else prior(th)
    v <- ll + pr + .logjac(e, dist)
    if (is.finite(v)) v else -Inf
  }

  target <- if (k == 1) 0.44 else 0.234
  cur <- lp(eta)
  scl <- 1
  acc <- 0
  win <- 0
  out <- matrix(NA_real_, floor(iter / thin), k,
                dimnames = list(NULL, dist$par_names))
  j <- 0
  for (t in seq_len(burnin + iter)) {
    prop <- eta + scl * drop(L %*% rnorm(k))
    lpp <- lp(prop)
    if (is.finite(lpp) && (!is.finite(cur) || log(runif(1)) < lpp - cur)) {
      eta <- prop
      cur <- lpp
      win <- win + 1
      if (t > burnin) acc <- acc + 1
    }
    if (t <= burnin && t %% 100 == 0) {
      scl <- scl * exp(win / 100 - target)
      win <- 0
    }
    if (t > burnin && (t - burnin) %% thin == 0 && j < nrow(out)) {
      j <- j + 1
      out[j, ] <- .from_free(eta, dist)
    }
  }
  structure(list(estimate = colMeans(out), median = apply(out, 2, stats::median),
                 se = apply(out, 2, sd), method = "bayes", draws = out,
                 acceptance = acc / iter, mle = init, dist = dist, n = length(x)),
            class = "ss_fit")
}

#' @export
print.ss_fit <- function(x, ...) {
  cat("Method:", x$method, if (!is.null(x$base)) paste0("(base: ", x$base, ")"),
      " n =", x$n, "\n")
  tab <- data.frame(estimate = x$estimate)
  if (!is.null(x$se)) tab$se <- x$se[names(x$estimate)]
  else if (!is.null(x$vcov)) tab$se <- sqrt(pmax(diag(x$vcov), 0))
  print(tab)
  if (!is.null(x$acceptance)) cat("Acceptance rate:", round(x$acceptance, 3), "\n")
  invisible(x)
}

# ---------------------------------------------------------------------------

.check_data <- function(x, dist) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (any(x < dist$lower | x > dist$upper))
    stop("data outside the support (", dist$lower, ", ", dist$upper, ")")
  if (length(x) <= length(dist$par_names)) stop("too few observations")
  x
}

.objective <- function(method, x, dist) {
  xs <- sort(x)
  n <- length(xs)
  i <- seq_len(n)
  switch(method,
    mle = function(th) -sum(.dens(dist, xs, th, log = TRUE)),
    mps = function(th) {
      D <- diff(c(0, .cdf(dist, xs, th), 1))
      -sum(log(pmax(D, 1e-300)))
    },
    lse = function(th) sum((.cdf(dist, xs, th) - i / (n + 1))^2),
    wlse = function(th) {
      w <- (n + 1)^2 * (n + 2) / (i * (n - i + 1))
      sum(w * (.cdf(dist, xs, th) - i / (n + 1))^2)
    },
    cvm = function(th) 1 / (12 * n) + sum((.cdf(dist, xs, th) - (2 * i - 1) / (2 * n))^2),
    ad = function(th) {
      Fx <- pmin(pmax(.cdf(dist, xs, th), 1e-300), 1 - 1e-15)
      -n - mean((2 * i - 1) * (log(Fx) + log(1 - rev(Fx))))
    },
    stop("unknown method: ", method))
}

.free_obj <- function(obj, dist) {
  function(eta) {
    th <- .from_free(eta, dist)
    v <- tryCatch(suppressWarnings(obj(th)), error = function(e) Inf)
    if (length(v) != 1 || !is.finite(v)) 1e100 else v
  }
}

.optimize <- function(fn, eta0) {
  if (length(eta0) == 1) {
    return(optim(eta0, fn, method = "Brent", lower = eta0 - 15, upper = eta0 + 15,
                 control = list(reltol = 1e-12)))
  }
  o <- optim(eta0, fn, method = "Nelder-Mead", control = list(maxit = 3000, reltol = 1e-11))
  o2 <- tryCatch(optim(o$par, fn, method = "BFGS", control = list(maxit = 500, reltol = 1e-12)),
                 error = function(e) NULL)
  if (!is.null(o2) && o2$value <= o$value) o <- o2
  o
}

.fit_point <- function(x, dist, method, start = NULL, screen = TRUE, n_random = 30) {
  method <- match.arg(method, .point_methods)
  obj <- .free_obj(.objective(method, x, dist), dist)
  if (!is.null(start)) start <- .name_par(start, dist)
  if (screen || is.null(start)) {
    base <- if (!is.null(start)) start else dist$default
    e0 <- .to_free(base, dist)
    k <- length(e0)
    cand <- rbind(e0, matrix(rep(e0, each = n_random) + rnorm(n_random * k, sd = 2),
                             ncol = k))
    scr <- .free_obj(.objective("lse", x, dist), dist)
    sv <- apply(cand, 1, scr)
    cand <- cand[order(sv)[seq_len(min(3, nrow(cand)))], , drop = FALSE]
    if (!is.null(start)) cand <- rbind(.to_free(start, dist), cand)
  } else {
    cand <- matrix(.to_free(start, dist), nrow = 1)
  }
  fits <- lapply(seq_len(nrow(cand)), function(j)
    tryCatch(.optimize(obj, cand[j, ]), error = function(e) NULL))
  fits <- Filter(Negate(is.null), fits)
  if (!length(fits)) stop("optimisation failed")
  o <- fits[[which.min(vapply(fits, function(f) f$value, numeric(1)))]]
  if (o$value >= 1e100) stop("no valid parameter value found")
  theta <- .from_free(o$par, dist)
  V <- Vf <- NULL
  if (method == "mle") {
    H <- tryCatch(optimHess(o$par, obj), error = function(e) NULL)
    Vf <- if (!is.null(H)) tryCatch(solve(H), error = function(e) NULL) else NULL
    if (!is.null(Vf)) {
      J <- .jac(o$par, dist)
      V <- Vf * outer(J, J)
      dimnames(V) <- dimnames(Vf) <- list(dist$par_names, dist$par_names)
    }
  }
  structure(list(estimate = theta, method = method, eta = o$par, value = o$value,
                 convergence = o$convergence, vcov = V, vcov_free = Vf,
                 dist = dist, n = length(x)),
            class = "ss_fit")
}
