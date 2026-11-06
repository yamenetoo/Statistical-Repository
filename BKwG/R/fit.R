## Estimation of BKw-G models ---------------------------------------------------

bkwg_submodels <- function() {
  data.frame(
    submodel = c("BKwG", "KwG", "BetaG", "EKwG", "BEG", "EG", "G"),
    name = c("Beta Kumaraswamy-G", "Kumaraswamy-G", "Beta-G",
             "Exponentiated Kumaraswamy-G", "Beta exponentiated-G",
             "Exponentiated-G", "Baseline G"),
    restriction = c("none", "alpha = beta = 1", "a = b = 1", "beta = 1",
                    "b = 1", "b = alpha = beta = 1",
                    "a = b = alpha = beta = 1"),
    free_shapes = c(4L, 2L, 2L, 3L, 3L, 1L, 0L),
    stringsAsFactors = FALSE)
}

.submodel_fixed <- function(submodel) {
  switch(submodel,
         BKwG = numeric(0),
         KwG = c(alpha = 1, beta = 1),
         BetaG = c(a = 1, b = 1),
         EKwG = c(beta = 1),
         BEG = c(b = 1),
         EG = c(b = 1, alpha = 1, beta = 1),
         G = c(a = 1, b = 1, alpha = 1, beta = 1),
         stop("unknown submodel '", submodel, "'", call. = FALSE))
}

.method_label <- function(method) {
  switch(method,
         mle = "maximum likelihood",
         mps = "maximum product of spacings",
         ad = "minimum Anderson-Darling distance",
         cvm = "minimum Cramer-von Mises distance",
         lse = "ordinary least squares",
         wlse = "weighted least squares")
}

.assemble <- function(th, fx, allnames, free) {
  v <- rep(NA_real_, length(allnames))
  names(v) <- allnames
  if (length(fx)) v[names(fx)] <- fx
  if (length(free)) v[free] <- th
  v
}

.make_objective <- function(x, B, fx, free, allnames, method) {
  xs <- sort(x)
  n <- length(x)
  i <- seq_len(n)
  bnames <- names(B$par)
  function(th) {
    v <- .assemble(th, fx, allnames, free)
    a <- v[["a"]]; b <- v[["b"]]; al <- v[["alpha"]]; be <- v[["beta"]]
    if (!.shape_ok(a, b, al, be)) return(Inf)
    bp <- v[bnames]
    if (length(bp) &&
        (any(!is.finite(bp)) || any(bp <= B$lower) || any(bp >= B$upper)))
      return(Inf)
    Fx <- function(lower.tail = TRUE, log.p = FALSE)
      .pbkwg_raw(xs, a, b, al, be, B, bp, lower.tail = lower.tail,
                 log.p = log.p)
    val <- switch(method,
      mle = -sum(.ldens(x, a, b, al, be, B, bp)),
      mps = {
        d <- diff(c(0, Fx(), 1))
        -sum(log(pmax(d, 1e-300)))
      },
      cvm = sum((Fx() - (2 * i - 1) / (2 * n))^2) + 1 / (12 * n),
      ad = -n - mean((2 * i - 1) *
                       (Fx(log.p = TRUE) + rev(Fx(FALSE, log.p = TRUE)))),
      lse = sum((Fx() - i / (n + 1))^2),
      wlse = sum((n + 1)^2 * (n + 2) / (i * (n - i + 1)) *
                   (Fx() - i / (n + 1))^2))
    if (is.finite(val)) val else Inf
  }
}

.candidates <- function(x, B, fx, free, allnames, obj, start, nstart) {
  bpar <- B$par
  bs <- tryCatch(unlist(B$start(x)), error = function(e) NULL)
  if (length(bpar) && !is.null(bs)) {
    nm <- intersect(names(bs), names(bpar))
    nm <- nm[is.finite(bs[nm]) & bs[nm] > B$lower[nm] & bs[nm] < B$upper[nm]]
    bpar[nm] <- bs[nm]
  }
  base <- c(a = 1, b = 1, alpha = 1, beta = 1, bpar)
  cand <- list()
  if (!is.null(start)) {
    st <- unlist(start)
    if (is.null(names(st)))
      stop("'start' must be a named vector or list", call. = FALSE)
    bad <- setdiff(names(st), free)
    if (length(bad))
      warning("ignoring start values for non-free parameter(s): ",
              paste(bad, collapse = ", "), call. = FALSE)
    use <- intersect(names(st), free)
    s0 <- base
    s0[use] <- st[use]
    cand[[1L]] <- s0[free]
  }
  grid <- expand.grid(a = c(0.5, 1, 2, 4), b = c(0.5, 1, 2, 4),
                      alpha = c(0.5, 1, 2), beta = c(0.5, 1, 2))
  pool <- lapply(seq_len(nrow(grid)), function(k) {
    s <- base
    s[.shape_names] <- unlist(grid[k, .shape_names], use.names = FALSE)
    s[free]
  })
  pool <- unique(pool)
  v <- vapply(pool, obj, numeric(1))
  ord <- order(v)
  ord <- ord[is.finite(v[ord])]
  cand <- c(cand, pool[ord[seq_len(min(nstart, length(ord)))]])
  if (!length(cand))
    stop("could not find finite starting values; try supplying 'start'",
         call. = FALSE)
  cand
}

.optimize <- function(z0, fz, control) {
  ctrl <- list(maxit = 5000, reltol = 1e-10)
  ctrl[names(control)] <- control
  r1 <- NULL
  if (length(z0) > 1L) {
    r1 <- tryCatch(stats::optim(z0, fz, method = "Nelder-Mead",
                                control = ctrl),
                   error = function(e) NULL)
    if (!is.null(r1) && is.finite(r1$value)) z0 <- r1$par
  }
  r2 <- tryCatch(stats::optim(z0, fz, method = "BFGS", control = ctrl),
                 error = function(e) NULL)
  best <- r2
  if (is.null(best) || (!is.null(r1) && r1$value < best$value)) best <- r1
  if (is.null(best) || !is.finite(best$value)) return(NULL)
  best
}

fit_bkwg <- function(x, baseline = "weibull", submodel = "BKwG",
                     fixed = NULL, start = NULL,
                     method = c("mle", "mps", "ad", "cvm", "lse", "wlse"),
                     nstart = 3, hessian = TRUE, control = list()) {
  cl <- match.call()
  dname <- paste(deparse(substitute(x)), collapse = " ")
  method <- match.arg(method)
  submodel <- match.arg(submodel, bkwg_submodels()$submodel)
  if (!is.numeric(x)) stop("'x' must be numeric")
  if (anyNA(x)) {
    warning("missing values removed")
    x <- x[!is.na(x)]
  }
  if (any(!is.finite(x))) stop("'x' must contain finite values only")
  if (length(x) < 5L) stop("at least 5 observations are needed")
  B <- .get_baseline(baseline)
  if (any(x < B$support[1L] | x > B$support[2L]))
    stop("some observations lie outside the support of the '", B$name,
         "' baseline", call. = FALSE)

  allnames <- c(.shape_names, names(B$par))
  fx <- .submodel_fixed(submodel)
  if (!is.null(fixed)) {
    fixed <- unlist(fixed)
    if (!is.numeric(fixed) || is.null(names(fixed)) ||
        any(!names(fixed) %in% allnames))
      stop("'fixed' must be a named vector/list with names among: ",
           paste(allnames, collapse = ", "), call. = FALSE)
    fx[names(fixed)] <- fixed
  }
  free <- setdiff(allnames, names(fx))
  lo_all <- c(a = 0, b = 0, alpha = 0, beta = 0, B$lower)
  up_all <- c(a = Inf, b = Inf, alpha = Inf, beta = Inf, B$upper)
  lo <- lo_all[free]
  up <- up_all[free]

  obj <- .make_objective(x, B, fx, free, allnames, method)

  if (!length(free)) {
    est <- numeric(0)
    val <- obj(est)
    conv <- 0L
    if (!is.finite(val)) stop("the fully fixed model has non-finite objective")
  } else {
    cands <- .candidates(x, B, fx, free, allnames, obj, start, nstart)
    best <- NULL
    for (s0 in cands) {
      if (any(s0 <= lo | s0 >= up)) next
      fz <- function(z) obj(.from_z(z, lo, up))
      r <- .optimize(.to_z(s0, lo, up), fz, control)
      if (!is.null(r) && (is.null(best) || r$value < best$value)) best <- r
    }
    if (is.null(best)) stop("optimisation failed for all starting values")
    est <- .from_z(best$par, lo, up)
    names(est) <- free
    val <- best$value
    conv <- best$convergence
    if (conv != 0L) warning("optimiser did not report convergence (code ",
                            conv, ")", call. = FALSE)
  }

  full <- .assemble(est, fx, allnames, free)
  bp <- full[names(B$par)]
  ll <- sum(.ldens(x, full[["a"]], full[["b"]], full[["alpha"]],
                   full[["beta"]], B, bp))
  k <- length(free)
  n <- length(x)

  V <- matrix(NA_real_, k, k, dimnames = list(free, free))
  if (hessian && k > 0L && method %in% c("mle", "mps")) {
    nll <- function(th) {
      v <- .assemble(th, fx, allnames, free)
      b0 <- v[names(B$par)]
      if (!.shape_ok(v[["a"]], v[["b"]], v[["alpha"]], v[["beta"]]) ||
          (length(b0) && (any(b0 <= B$lower) || any(b0 >= B$upper))))
        return(NA_real_)
      -sum(.ldens(x, v[["a"]], v[["b"]], v[["alpha"]], v[["beta"]], B, b0))
    }
    H <- tryCatch({
      if (requireNamespace("numDeriv", quietly = TRUE)) {
        numDeriv::hessian(nll, est)
      } else {
        stats::optimHess(est, nll)
      }
    }, error = function(e) NULL)
    if (!is.null(H) && all(is.finite(H))) {
      Vi <- tryCatch(solve(H), error = function(e) NULL)
      if (!is.null(Vi)) V[] <- Vi
    }
    if (all(is.na(V))) {
      warning("the observed information matrix could not be inverted; ",
              "standard errors are not available", call. = FALSE)
    } else if (any(diag(V) <= 0)) {
      warning("the observed information matrix is not positive definite; ",
              "standard errors are unreliable", call. = FALSE)
    }
  }
  se <- suppressWarnings(sqrt(diag(V)))
  names(se) <- free

  structure(list(
    estimate = full, coefficients = est, free = free, fixed = fx,
    fixed_user = fixed, vcov = V, se = se, loglik = ll, objective = val,
    n = n, k = k,
    AIC = -2 * ll + 2 * k,
    CAIC = -2 * ll + 2 * k * n / (n - k - 1),
    BIC = -2 * ll + k * log(n),
    HQIC = -2 * ll + 2 * k * log(log(n)),
    method = method, baseline = B, submodel = submodel,
    convergence = conv, x = x, data.name = dname, call = cl),
    class = "bkwgfit")
}

## arguments for the d/p/q functions at the fitted values
.fit_args <- function(fit) {
  e <- fit$estimate
  list(a = e[["a"]], b = e[["b"]], alpha = e[["alpha"]], beta = e[["beta"]],
       baseline = fit$baseline, par = e[names(fit$baseline$par)])
}
