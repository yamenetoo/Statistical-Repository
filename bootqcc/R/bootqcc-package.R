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
#' @importFrom stats optim quantile sd setNames runif
#' @importFrom graphics hist abline plot
#' @keywords internal
"_PACKAGE"
