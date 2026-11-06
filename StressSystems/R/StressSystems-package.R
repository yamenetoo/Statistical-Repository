#' StressSystems: Stress-Strength Reliability for Arbitrary Distributions
#'
#' Computes and estimates stress-strength reliability \eqn{R = P(Y < X)} for
#' any pair of continuous distributions, for single components and for
#' multi-stress, three-variable, s-out-of-k and cascade systems.
#'
#' Main functions:
#' \itemize{
#'   \item \code{\link{ss_dist}}: define a distribution (built-in, by name or
#'     user-defined).
#'   \item \code{\link{ss_R}}, \code{\link{ss_multistress}},
#'     \code{\link{ss_R3}}, \code{\link{ss_sk}}, \code{\link{ss_cascade}}:
#'     reliability measures.
#'   \item \code{\link{ss_fit}}, \code{\link{ss_jackknife}},
#'     \code{\link{ss_bootstrap}}, \code{\link{ss_bayes}}: estimation.
#'   \item \code{\link{ss_estimate_R}}: estimate of \eqn{R} with a
#'     confidence interval.
#'   \item \code{\link{ss_simulate}}, \code{\link{ss_simulate_grid}}: Monte
#'     Carlo studies.
#' }
#'
#' @references Kotz, S., Lumelskii, Y. and Pensky, M. (2003). \emph{The
#'   Stress-Strength Model and Its Generalizations: Theory and Applications}.
#'   World Scientific, Singapore.
#' @keywords internal
#' @importFrom stats integrate optim optimHess pbinom plogis qlogis qnorm
#'   quantile rnorm runif sd setNames uniroot
#' @importFrom utils modifyList
#' @importFrom parallel mclapply
"_PACKAGE"
