# ==== Design principles: see R/core.R ====
# registry.R is a pure convenience layer. No other file in this package
# calls into it -- every estimation/diagnostic function accepts ddist/pdist
# directly, so the package never restricts users to a fixed distribution
# list. This file only saves typing for a handful of common presets.

.utd_registry_env <- new.env(parent = emptyenv())

#' @title Register a base distribution under a short name
#' @param name character key
#' @param ddist,pdist,qdist,rdist distribution functions (qdist/rdist optional)
#' @export
#' @examples
#' utd_register("gamma", dgamma, pgamma, qgamma, rgamma)
#' utd_list_registered()
utd_register <- function(name, ddist, pdist, qdist = NULL, rdist = NULL) {
  if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
    stop("`name` must be a single non-empty character string.", call. = FALSE)
  }
  .utd_registry_env[[name]] <- list(ddist = ddist, pdist = pdist, qdist = qdist, rdist = rdist)
  invisible(name)
}

#' @title Retrieve a registered base distribution as a template
#' @param name character key previously registered
#' @return a list with elements ddist, pdist, qdist, rdist
#' @export
#' @examples
#' utd_register("gamma", dgamma, pgamma, qgamma, rgamma)
#' tpl <- utd_get("gamma")
#' d <- utd(tpl$ddist, tpl$pdist, tpl$qdist, tpl$rdist, params = list(shape = 2))
utd_get <- function(name) {
  if (!exists(name, envir = .utd_registry_env, inherits = FALSE)) {
    stop(sprintf("'%s' is not registered. See utd_list_registered().", name), call. = FALSE)
  }
  .utd_registry_env[[name]]
}

#' @title List currently registered distribution names
#' @export
#' @examples
#' utd_list_registered()
utd_list_registered <- function() {
  sort(ls(envir = .utd_registry_env))
}

#' @title Register a handful of common base distributions for convenience
#' @details Purely optional. Registers gamma, weibull, lognormal, and
#'   exponential using base R's stats distributions.
#' @export
#' @examples
#' utd_preset_common()
#' utd_list_registered()
utd_preset_common <- function() {
  utd_register("gamma", dgamma, pgamma, qgamma, rgamma)
  utd_register("weibull", dweibull, pweibull, qweibull, rweibull)
  utd_register("lnorm", dlnorm, plnorm, qlnorm, rlnorm)
  utd_register("exp", dexp, pexp, qexp, rexp)
  invisible(utd_list_registered())
}
