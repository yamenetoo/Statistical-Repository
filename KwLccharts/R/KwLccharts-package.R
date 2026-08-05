#' KwLccharts: Bootstrap Percentile Control Charts for the Kumaraswamy-Lindley Distribution
#'
#' Provides the density, cumulative distribution, closed-form quantile
#' (via the Lambert W function) and random-generation functions for the
#' Kumaraswamy-Lindley (KwL) distribution, together with maximum-likelihood
#' fitting (\code{\link{fit_KwL}}) and a parametric bootstrap percentile
#' control chart (\code{\link{KwLBootCC}}, \code{\link{KwL_ARL}}), following
#' the methodology of Al-Lami, Jabbari Khamnei and Heydari (2025, PLOS ONE)
#' for the Lindley-Geometric distribution, adapted here to the
#' three-parameter Kumaraswamy-Lindley family of Cordeiro and de Castro
#' (2011).
#'
#' @section Main functions:
#' \itemize{
#'   \item Distribution: \code{\link{dKwL}}, \code{\link{pKwL}},
#'     \code{\link{qKwL}}, \code{\link{rKwL}}
#'   \item Baseline Lindley: \code{\link{gLindley}}, \code{\link{GLindley}},
#'     \code{\link{QLindley}}
#'   \item Estimation: \code{\link{fit_KwL}}
#'   \item Control chart: \code{\link{KwLBootCC}}, \code{\link{KwL_ARL}}
#' }
#'
#' @keywords internal
"_PACKAGE"
