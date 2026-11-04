#' bootqcc: parametric bootstrap control charts for any distribution
#'
#' Phase I limits for a quantile (or any statistic) of an arbitrary parametric
#' distribution via parametric bootstrap and maximum likelihood, with Monte
#' Carlo ARL evaluation and Phase II monitoring.
#'
#' @section Main functions:
#' [bootqcc()] builds a chart, [arl()] re-evaluates its run length (for example
#' under a parameter shift), [bootqcc_monitor()] applies it to new samples and
#' [bootqcc_grid()] runs a simulation study over many settings.
#'
#' @references
#' Efron, B. and Tibshirani, R. J. (1993). \emph{An Introduction to the Bootstrap}.
#' Chapman and Hall, New York.
#'
#' Montgomery, D. C. (2019). \emph{Introduction to Statistical Quality Control},
#' 8th edition. Wiley, Hoboken, NJ.
#'
#' Al-Lami, M. A. H., Jabbari Khamnei, H. and Heydari, A. A. (2025). Bootstrap
#' control chart for the Lindley-Geometric distribution. \emph{PLOS ONE}.
#' \doi{10.1371/journal.pone.0316449}
#'
#' @importFrom stats optim quantile sd setNames runif
#' @importFrom graphics hist abline plot
#' @keywords internal
"_PACKAGE"
