#' UnitTrans: distributions on (0,1) via transformations of baseline laws
#'
#' Given a baseline variable Y with any `d/p/q/r` functions and a monotone
#' transformation V = g(Y) into (0,1), the package derives the density,
#' distribution and quantile functions of V, simulates from it, fits it by
#' maximum likelihood, ranks many candidate models, tests goodness of fit,
#' and runs simulation studies.
#' @importFrom stats coef vcov logLik confint
#' @keywords internal
"_PACKAGE"
