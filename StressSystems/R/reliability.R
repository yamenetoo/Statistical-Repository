#' Stress-strength reliability R = P(cY < X)
#'
#' Computes \eqn{R = P(cY < X) = E_X[F_Y(X / c)]} for independent strength
#' \eqn{X} and stress \eqn{Y}. With \code{c = 1} this is the usual
#' \eqn{P(Y < X)}; \code{c} scales the stress (amplification \code{c > 1},
#' attenuation \code{c < 1}).
#'
#' By default the expectation is computed on the probability scale,
#' \eqn{R = \int_0^1 F_Y(Q_X(u) / c)\,du}. The integrand is bounded in [0, 1],
#' so densities that are unbounded at the edge of the support and heavy tails
#' do not cause trouble. Without a quantile function the density is
#' integrated piecewise between quantiles of \eqn{X}. If integration fails,
#' a Monte Carlo estimate is used and a warning is given.
#'
#' @param strength,stress Objects from \code{\link{ss_dist}}.
#' @param par_x,par_y Parameter vectors; default to those stored in the objects.
#' @param c Positive stress scaling factor.
#' @param method \code{"auto"}, \code{"quantile"}, \code{"density"} or \code{"mc"}.
#' @param rel.tol Relative tolerance for \code{integrate}.
#' @return A number in [0, 1].
#' @examples
#' X <- ss_dist("exp", par = c(rate = 2))
#' Y <- ss_dist("exp", par = c(rate = 3))
#' ss_R(X, Y)        # exact value 3 / 5
#' @export
ss_R <- function(strength, stress, par_x = NULL, par_y = NULL, c = 1,
                 method = "auto", rel.tol = 1e-8) {
  par_x <- .need_par(strength, par_x)
  par_y <- .need_par(stress, par_y)
  .check_c(c)
  g <- function(x) .cdf(stress, x / c, par_y)
  .expect(strength, par_x, g, method, rel.tol)
}

.check_c <- function(c) {
  if (length(c) != 1 || !is.finite(c) || c <= 0) stop("'c' must be a single positive number")
}

# E[g(X)] for X ~ dist(par)
.expect <- function(dist, par, g, method = "auto", rel.tol = 1e-8, mc_n = 1e5) {
  method <- match.arg(method, c("auto", "quantile", "density", "mc"))
  if (method == "auto") method <- if (!is.null(dist$q)) "quantile" else "density"
  mc <- function() mean(g(.rand(dist, mc_n, par)))
  if (method == "mc") return(mc())
  run <- function(tol) {
    if (method == "quantile") {
      integrate(function(u) g(.quant(dist, u, par)), 0, 1,
                rel.tol = tol, subdivisions = 2000L)$value
    } else {
      .integrate_density(dist, par, g, tol)
    }
  }
  val <- tryCatch(run(rel.tol),
                  error = function(e) tryCatch(run(1e-5), error = function(e2) NA_real_))
  if (!is.finite(val)) {
    warning("numerical integration failed; using a Monte Carlo estimate (",
            mc_n, " draws)")
    val <- mc()
  }
  min(max(val, 0), 1)
}

.integrate_density <- function(dist, par, g, rel.tol) {
  br <- .quant(dist, c(1e-6, 0.01, 0.25, 0.5, 0.75, 0.99, 1 - 1e-6), par)
  br <- sort(unique(c(dist$lower, br[is.finite(br)], dist$upper)))
  f <- function(x) g(x) * .dens(dist, x, par)
  sum(vapply(seq_len(length(br) - 1), function(i)
    integrate(f, br[i], br[i + 1], rel.tol = rel.tol, subdivisions = 2000L)$value,
    numeric(1)))
}
