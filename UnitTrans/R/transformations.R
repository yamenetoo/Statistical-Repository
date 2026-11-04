# Registry of transformations V = g(Y) mapping a baseline variable Y into (0,1)

.make_trans <- function(name, formula, inverse, domain, decreasing,
                        g, ginv, logjac) {
  structure(list(name = name, formula = formula, inverse = inverse,
                 domain = domain, decreasing = decreasing,
                 g = g, ginv = ginv, logjac = logjac),
            class = "unitrans")
}

.registry <- function() {
  list(
    rational = .make_trans(
      "rational", "V = Y / (1 + Y)", "Y = V / (1 - V)", "nonneg", FALSE,
      g = function(y) y / (1 + y),
      ginv = function(v) v / (1 - v),
      logjac = function(v) -2 * log1p(-v)),
    exponential = .make_trans(
      "exponential", "V = exp(-Y)", "Y = -log(V)", "positive", TRUE,
      g = function(y) exp(-y),
      ginv = function(v) -log(v),
      logjac = function(v) -log(v)),
    logit = .make_trans(
      "logit", "V = 1 / (1 + exp(-Y))", "Y = log(V / (1 - V))", "real", FALSE,
      g = function(y) stats::plogis(y),
      ginv = function(v) stats::qlogis(v),
      logjac = function(v) -log(v) - log1p(-v)),
    # V = 2 / (1 + exp(Y)) = 1 - tanh(Y / 2); previously (wrongly) named "arcsech".
    tanh = .make_trans(
      "tanh", "V = 2 / (1 + exp(Y)) = 1 - tanh(Y / 2)",
      "Y = log(2 / V - 1) = 2 atanh(1 - V)", "positive", TRUE,
      g = function(y) 2 * stats::plogis(-y),
      ginv = function(v) log(2 - v) - log(v),
      logjac = function(v) log(2) - log(v) - log(2 - v)),
    probit = .make_trans(
      "probit", "V = pnorm(Y)", "Y = qnorm(V)", "real", FALSE,
      g = function(y) stats::pnorm(y),
      ginv = function(v) stats::qnorm(v),
      logjac = function(v) -stats::dnorm(stats::qnorm(v), log = TRUE)),
    cloglog = .make_trans(
      "cloglog", "V = 1 - exp(-exp(Y))", "Y = log(-log(1 - V))", "real", FALSE,
      g = function(y) -expm1(-exp(y)),
      ginv = function(v) log(-log1p(-v)),
      logjac = function(v) -log1p(-v) - log(-log1p(-v)))
  )
}

#' List the available transformations
#'
#' @return A data frame with the name, forward formula, inverse formula,
#'   admissible baseline domain and monotonicity of each transformation.
#' @examples
#' trans_list()
#' @export
trans_list <- function() {
  r <- .registry()
  data.frame(
    name = names(r),
    forward = vapply(r, `[[`, "", "formula"),
    inverse = vapply(r, `[[`, "", "inverse"),
    domain = vapply(r, `[[`, "", "domain"),
    decreasing = vapply(r, `[[`, TRUE, "decreasing"),
    row.names = NULL, stringsAsFactors = FALSE)
}

#' Get a transformation object by name
#'
#' @param name One of the names returned by [trans_list()]. The old name
#'   `"arcsech"` is still accepted (with a warning) and refers to `"tanh"`:
#'   the map `V = 2/(1 + exp(Y))` equals `1 - tanh(Y/2)` and is not the
#'   inverse hyperbolic secant.
#' @return An object of class `unitrans` with elements `g` (forward map),
#'   `ginv` (inverse map), `logjac` (log |dY/dV|), `domain` and `decreasing`.
#' @examples
#' tr <- get_transformation("logit")
#' tr$g(0)
#' @export
get_transformation <- function(name) {
  if (inherits(name, "unitrans")) return(name)
  r <- .registry()
  if (is.character(name) && length(name) == 1L && identical(name, "arcsech")) {
    warning("Transformation 'arcsech' was renamed 'tanh' (V = 1 - tanh(Y/2)); ",
            "using 'tanh'.", call. = FALSE)
    name <- "tanh"
  }
  if (!is.character(name) || length(name) != 1L || !name %in% names(r))
    stop("Unknown transformation '", format(name)[1], "'. Available: ",
         paste(names(r), collapse = ", "), call. = FALSE)
  r[[name]]
}

#' @export
print.unitrans <- function(x, ...) {
  cat("Transformation:", x$name, "\n  forward:", x$formula,
      "\n  inverse:", x$inverse, "\n  baseline domain:", x$domain,
      "\n  monotone:", if (x$decreasing) "decreasing" else "increasing", "\n")
  invisible(x)
}

#' Inverse transformation
#'
#' Maps values of V in (0,1) back to the baseline scale Y.
#' @param v Numeric vector in (0,1).
#' @param trans Transformation name or object.
#' @return Numeric vector Y = g^{-1}(V).
#' @examples
#' inverse_transformation(0.3, "logit")
#' @export
inverse_transformation <- function(v, trans) get_transformation(trans)$ginv(v)

#' Forward transformation
#'
#' @param y Numeric vector on the baseline scale.
#' @param trans Transformation name or object.
#' @return Numeric vector V = g(Y).
#' @examples
#' forward_transformation(1.2, "rational")
#' @export
forward_transformation <- function(y, trans) get_transformation(trans)$g(y)
