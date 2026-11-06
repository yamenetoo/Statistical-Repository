#' Reliability in a multi-stress environment
#'
#' A component of strength \eqn{X} is exposed to \eqn{k} independent stresses
#' \eqn{Y_1, \dots, Y_k}.
#'
#' With \code{components = "shared"} (default) one component must withstand
#' all stresses: \eqn{R = P(cY_j < X \ \forall j) = E_X[\prod_j F_{Y_j}(X/c)]}.
#' With \code{components = "independent"} each stress acts on its own
#' independent copy of the component and \eqn{R = \prod_j P(cY_j < X_j)}.
#' The two differ because in the shared case the events share the same
#' \eqn{X}.
#'
#' @param strength An \code{ss_dist}.
#' @param stresses Either one \code{ss_dist} (same family for every stress,
#'   with one row of \code{par_y} per stress) or a list of \code{ss_dist}.
#' @param par_x Strength parameters.
#' @param par_y Matrix or data frame with one row per stress, or a list of
#'   parameter vectors. Defaults to the parameters stored in \code{stresses}.
#' @param c Stress scaling factor.
#' @param components \code{"shared"} or \code{"independent"}.
#' @param ... Passed to the integration (\code{method}, \code{rel.tol}).
#' @return A number in [0, 1].
#' @examples
#' X <- ss_dist("toppleone", par = c(shape = 0.7))
#' Y <- ss_dist("genrayleigh")
#' st <- data.frame(shape = c(1, 2), rate = c(2, 3))
#' ss_multistress(X, Y, par_y = st)
#' ss_multistress(X, Y, par_y = st, components = "independent")
#' @export
ss_multistress <- function(strength, stresses, par_x = NULL, par_y = NULL, c = 1,
                           components = c("shared", "independent"), ...) {
  components <- match.arg(components)
  par_x <- .need_par(strength, par_x)
  .check_c(c)
  st <- .stage_list(stresses, par_y)
  if (components == "independent")
    return(prod(vapply(st, function(s) ss_R(strength, s$dist, par_x, s$par, c = c, ...),
                       numeric(1))))
  g <- function(x) {
    out <- rep(1, length(x))
    for (s in st) out <- out * .cdf(s$dist, x / c, s$par)
    out
  }
  .expect(strength, par_x, g, ...)
}

#' Three-variable reliability P(Z < Y < X)
#'
#' Probability that the stress \eqn{Y} lies between a lower strength \eqn{Z}
#' and an upper strength \eqn{X}:
#' \eqn{R_3 = E_Y[F_Z(Y)\,(1 - F_X(Y))]}. Note that this is not the product
#' \eqn{P(Y < X)\,P(Z < Y)}, since both events involve the same \eqn{Y}.
#'
#' @param lower Distribution of the lower bound \eqn{Z}.
#' @param stress Distribution of \eqn{Y}.
#' @param upper Distribution of the upper bound \eqn{X}.
#' @param par_z,par_y,par_x Parameter vectors.
#' @param ... Passed to the integration.
#' @return A number in [0, 1].
#' @examples
#' Z <- ss_dist("exp", par = c(rate = 1))
#' Y <- ss_dist("exp", par = c(rate = 2))
#' X <- ss_dist("exp", par = c(rate = 0.5))
#' ss_R3(Z, Y, X)          # exact value 2/2.5 - 2/3.5
#' @export
ss_R3 <- function(lower, stress, upper, par_z = NULL, par_y = NULL, par_x = NULL, ...) {
  par_z <- .need_par(lower, par_z)
  par_y <- .need_par(stress, par_y)
  par_x <- .need_par(upper, par_x)
  g <- function(y) .cdf(lower, y, par_z) * (1 - .cdf(upper, y, par_x))
  .expect(stress, par_y, g, ...)
}

#' Reliability of an s-out-of-k system
#'
#' The system has \eqn{k} independent components with strengths
#' \eqn{X_1, \dots, X_k} from the same distribution and works if at least
#' \eqn{s} of them exceed the common stress \eqn{cY}:
#' \deqn{R_{s,k} = \sum_{i=s}^{k} {k \choose i}
#'   E_Y\big[(1 - F_X(cY))^{i} F_X(cY)^{k-i}\big].}
#'
#' @param strength,stress Objects from \code{\link{ss_dist}}.
#' @param s,k Integers with \eqn{1 \le s \le k}.
#' @param par_x,par_y Parameter vectors.
#' @param c Stress scaling factor.
#' @param ... Passed to the integration.
#' @return A number in [0, 1].
#' @references Bhattacharyya, G. K. and Johnson, R. A. (1974). Estimation of
#'   reliability in a multicomponent stress-strength model. \emph{Journal of
#'   the American Statistical Association}.
#'   \doi{10.2307/2286173}
#' @examples
#' X <- ss_dist("toppleone", par = c(shape = 0.7))
#' Y <- ss_dist("genrayleigh", par = c(shape = 1, rate = 2))
#' sapply(1:3, function(s) ss_sk(X, Y, s = s, k = 3))
#' @export
ss_sk <- function(strength, stress, s, k, par_x = NULL, par_y = NULL, c = 1, ...) {
  if (!(s >= 1 && k >= s)) stop("need 1 <= s <= k")
  par_x <- .need_par(strength, par_x)
  par_y <- .need_par(stress, par_y)
  .check_c(c)
  g <- function(y) pbinom(s - 1, k, 1 - .cdf(strength, c * y, par_x), lower.tail = FALSE)
  .expect(stress, par_y, g, ...)
}

#' Reliability of a cascade system
#'
#' An \eqn{n}-stage cascade (standby) system: when stage \eqn{i - 1} fails,
#' stage \eqn{i} takes over and its stress is scaled by \eqn{c^{i-1}}. With
#' \eqn{P_i = P(c^{i-1} Y_i < X_i)} and independent stages,
#' \deqn{R = \sum_{i=1}^{n} P_i \prod_{j < i} (1 - P_j).}
#'
#' @param strength,stress One \code{ss_dist} each (same family at every stage)
#'   or lists of \code{ss_dist} with one element per stage.
#' @param n Number of stages.
#' @param c Attenuation (\code{c < 1}) or amplification (\code{c > 1}) factor.
#' @param par_x,par_y One parameter vector (same at every stage), a matrix or
#'   data frame with one row per stage, or a list of vectors.
#' @param ... Passed to \code{\link{ss_R}}.
#' @return An object of class \code{"ss_cascade"} with the system reliability
#'   \code{R}, the stage probabilities \code{P} and each stage's contribution.
#' @examples
#' X <- ss_dist("toppleone", par = c(shape = 0.7))
#' Y <- ss_dist("genrayleigh", par = c(shape = 1, rate = 2))
#' ss_cascade(X, Y, n = 3, c = 0.8)
#' @export
ss_cascade <- function(strength, stress, n, c = 1, par_x = NULL, par_y = NULL, ...) {
  .check_c(c)
  sx <- .stage_list(strength, par_x, n)
  sy <- .stage_list(stress, par_y, n)
  P <- vapply(seq_len(n), function(i)
    ss_R(sx[[i]]$dist, sy[[i]]$dist, sx[[i]]$par, sy[[i]]$par, c = c^(i - 1), ...),
    numeric(1))
  survive_before <- cumprod(c(1, 1 - P))[seq_len(n)]
  contrib <- P * survive_before
  structure(list(R = sum(contrib), P = P, contribution = contrib, c = c),
            class = "ss_cascade")
}

#' @export
print.ss_cascade <- function(x, ...) {
  cat("Cascade system, c =", x$c, "\n")
  print(data.frame(stage = seq_along(x$P), P = x$P, contribution = x$contribution),
        row.names = FALSE)
  cat("System reliability R =", format(x$R, digits = 6), "\n")
  invisible(x)
}

# Turn (dist or list of dists) + (vector / matrix / list of pars) into a list
# of list(dist, par). If n is given, a single dist and single par are recycled.
.stage_list <- function(dists, pars, n = NULL) {
  single <- inherits(dists, "ss_dist")
  if (is.matrix(pars) || is.data.frame(pars)) {
    pars <- lapply(seq_len(nrow(pars)), function(i) {
      v <- as.numeric(unlist(pars[i, ]))
      names(v) <- colnames(pars)
      v
    })
  } else if (!is.null(pars) && !is.list(pars)) {
    pars <- list(pars)
  }
  if (single) {
    m <- if (!is.null(pars)) length(pars) else 1
    if (!is.null(n) && m == 1) m <- n
    if (is.null(pars)) pars <- list(NULL)
    dists <- rep(list(dists), m)
    if (length(pars) == 1) pars <- rep(pars, m)
  } else {
    if (is.null(pars)) pars <- rep(list(NULL), length(dists))
    if (length(pars) == 1 && length(dists) > 1) pars <- rep(pars, length(dists))
  }
  if (length(pars) != length(dists)) stop("number of parameter sets does not match")
  if (!is.null(n) && length(dists) != n) stop("expected ", n, " stages")
  lapply(seq_along(dists), function(i)
    list(dist = dists[[i]], par = .need_par(dists[[i]], pars[[i]])))
}
