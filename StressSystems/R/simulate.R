#' Monte Carlo study of estimators of R
#'
#' Repeatedly draws strength and stress samples, estimates \eqn{R} with every
#' method in \code{methods} and reports bias, mean squared error, interval
#' coverage and average interval length for \eqn{R}, plus bias and MSE of each
#' parameter. Failed replications are counted, not silently set to zero.
#'
#' @param strength,stress Objects from \code{\link{ss_dist}}.
#' @param par_x,par_y True parameter vectors.
#' @param n Sample sizes: a vector (equal sizes for X and Y) or a two-column
#'   matrix of \code{(n_x, n_y)} pairs.
#' @param MC Number of replications.
#' @param methods Estimation methods (see \code{\link{ss_fit}}).
#' @param c Stress scaling factor.
#' @param ci Interval type passed to \code{\link{ss_estimate_R}}.
#' @param level Confidence level.
#' @param B Bootstrap size (bootstrap estimator and \code{ci = "boot"}).
#' @param iter,burnin MCMC settings for \code{"bayes"}.
#' @param prior Log prior for \code{"bayes"}.
#' @param seed Optional seed.
#' @param cores Number of cores (forked processes; Linux and macOS only).
#' @param keep_raw Keep every replication in the result.
#' @return An object of class \code{"ss_sim"} with \code{summary} (for R),
#'   \code{parameters} (per-parameter bias and MSE), \code{true_R} and
#'   optionally \code{raw}.
#' @export
ss_simulate <- function(strength, stress, par_x, par_y, n, MC = 1000,
                        methods = c("mle", "jackknife", "bayes"), c = 1,
                        ci = "auto", level = 0.95, B = 200, iter = 5000,
                        burnin = 1000, prior = NULL, seed = NULL, cores = 1,
                        keep_raw = TRUE) {
  par_x <- .name_par(par_x, strength)
  par_y <- .name_par(par_y, stress)
  true_R <- ss_R(strength, stress, par_x, par_y, c = c)
  nn <- if (is.matrix(n)) n else cbind(n, n)
  if (!is.null(seed)) set.seed(seed)
  run <- if (cores > 1 && .Platform$OS.type == "unix") {
    function(f) mclapply(seq_len(MC), f, mc.cores = cores)
  } else {
    function(f) lapply(seq_len(MC), f)
  }

  raw <- do.call(rbind, lapply(seq_len(nrow(nn)), function(s) {
    nx <- nn[s, 1]
    ny <- nn[s, 2]
    one_rep <- function(r) {
      x <- .rand(strength, nx, par_x)
      y <- .rand(stress, ny, par_y)
      do.call(rbind, lapply(methods, function(m) {
        res <- tryCatch(suppressWarnings(
          ss_estimate_R(x, y, strength, stress, method = m, ci = ci, level = level,
                        c = c, B = B, iter = iter, burnin = burnin, prior = prior)),
          error = function(e) NULL)
        .sim_row(r, nx, ny, m, res, strength, stress)
      }))
    }
    reps <- run(one_rep)
    reps <- reps[vapply(reps, is.data.frame, logical(1))]
    do.call(rbind, reps)
  }))

  truth <- c(setNames(par_x, paste0("x_", names(par_x))),
             setNames(par_y, paste0("y_", names(par_y))))
  sm <- .summarise_sim(raw, true_R, truth)
  structure(list(summary = sm$R, parameters = sm$par, true_R = true_R,
                 par_x = par_x, par_y = par_y, c = c,
                 raw = if (keep_raw) raw else NULL),
            class = "ss_sim")
}

#' Monte Carlo study over a grid of parameter values
#'
#' Runs \code{\link{ss_simulate}} for every combination of a row of
#' \code{grid_x} with a row of \code{grid_y} and stacks the summaries, with the
#' true parameters added as columns.
#'
#' @inheritParams ss_simulate
#' @param grid_x,grid_y Data frames (or matrices with column names) of
#'   strength and stress parameter values, one scenario per row.
#' @param ... Passed to \code{\link{ss_simulate}}.
#' @export
ss_simulate_grid <- function(strength, stress, grid_x, grid_y, n, ...,
                             keep_raw = FALSE) {
  grid_x <- as.data.frame(grid_x)
  grid_y <- as.data.frame(grid_y)
  combos <- expand.grid(ix = seq_len(nrow(grid_x)), iy = seq_len(nrow(grid_y)))
  res <- lapply(seq_len(nrow(combos)), function(j) {
    px <- unlist(grid_x[combos$ix[j], , drop = FALSE])
    py <- unlist(grid_y[combos$iy[j], , drop = FALSE])
    s <- ss_simulate(strength, stress, px, py, n, ..., keep_raw = keep_raw)
    lab <- as.data.frame(as.list(c(setNames(px, paste0("x_", names(px))),
                                   setNames(py, paste0("y_", names(py))))))
    list(R = cbind(lab[rep(1, nrow(s$summary)), , drop = FALSE], s$summary),
         par = cbind(lab[rep(1, nrow(s$parameters)), , drop = FALSE], s$parameters),
         raw = if (keep_raw) cbind(lab[rep(1, nrow(s$raw)), , drop = FALSE], s$raw))
  })
  out <- list(summary = do.call(rbind, lapply(res, `[[`, "R")),
              parameters = do.call(rbind, lapply(res, `[[`, "par")))
  rownames(out$summary) <- rownames(out$parameters) <- NULL
  if (keep_raw) out$raw <- do.call(rbind, lapply(res, `[[`, "raw"))
  out
}

#' @export
print.ss_sim <- function(x, ...) {
  cat("Monte Carlo study of R = P(cY < X), c =", x$c, ", true R =",
      format(x$true_R, digits = 6), "\n")
  print(x$summary, row.names = FALSE, digits = 4)
  invisible(x)
}

.sim_row <- function(r, nx, ny, m, res, strength, stress) {
  px <- setNames(rep(NA_real_, length(strength$par_names)),
                 paste0("x_", strength$par_names))
  py <- setNames(rep(NA_real_, length(stress$par_names)),
                 paste0("y_", stress$par_names))
  Rh <- lo <- up <- NA_real_
  if (!is.null(res)) {
    Rh <- res$estimate
    lo <- res$lower
    up <- res$upper
    px[] <- res$fit_x$estimate[strength$par_names]
    py[] <- res$fit_y$estimate[stress$par_names]
  }
  data.frame(rep = r, n_x = nx, n_y = ny, method = m, R_hat = Rh, lower = lo,
             upper = up, as.list(px), as.list(py), stringsAsFactors = FALSE)
}

.summarise_sim <- function(raw, true_R, truth) {
  groups <- split(raw, list(raw$n_x, raw$n_y, raw$method), drop = TRUE)
  R <- do.call(rbind, lapply(groups, function(d) {
    ok <- is.finite(d$R_hat)
    e <- d$R_hat[ok]
    ci_ok <- is.finite(d$lower) & is.finite(d$upper)
    data.frame(n_x = d$n_x[1], n_y = d$n_y[1], method = d$method[1], R = true_R,
               mean = mean(e), bias = mean(e) - true_R, mse = mean((e - true_R)^2),
               coverage = if (any(ci_ok)) mean(d$lower[ci_ok] <= true_R &
                                                 true_R <= d$upper[ci_ok]) else NA_real_,
               length = if (any(ci_ok)) mean(d$upper[ci_ok] - d$lower[ci_ok]) else NA_real_,
               failed = sum(!ok), stringsAsFactors = FALSE)
  }))
  par <- do.call(rbind, lapply(groups, function(d) {
    do.call(rbind, lapply(names(truth), function(p) {
      e <- d[[p]][is.finite(d[[p]])]
      data.frame(n_x = d$n_x[1], n_y = d$n_y[1], method = d$method[1],
                 parameter = p, true = truth[[p]], mean = mean(e),
                 bias = mean(e) - truth[[p]], mse = mean((e - truth[[p]])^2),
                 stringsAsFactors = FALSE)
    }))
  }))
  R <- R[order(R$n_x, R$n_y, R$method), ]
  par <- par[order(par$n_x, par$n_y, par$method, par$parameter), ]
  rownames(R) <- rownames(par) <- NULL
  list(R = R, par = par)
}
