## Fuzzification (Ross 2000, Section 3) -------------------------------------

fuzzify <- function(value, spread, relative = FALSE,
                    shape = c("triangular", "trapezoidal", "gaussian"),
                    core = 0, name = NULL) {
  ## Turn crisp readings into fuzzy sets that express their imprecision,
  ## e.g. a voltmeter reading +/- 1% (Figure 1.6).
  shape <- match.arg(shape)
  value <- as.numeric(value)
  if (!length(value) || anyNA(value)) stop("'value' must be numeric without NA")
  spread <- rep_len(as.numeric(spread), length(value))
  core <- rep_len(as.numeric(core), length(value))
  if (relative) {
    spread <- spread * abs(value)
    core <- core * abs(value)
  }
  if (any(spread <= 0)) stop("'spread' must be positive")
  if (any(core < 0) || any(core >= spread)) stop("need 0 <= core < spread")
  nm <- if (is.null(name)) format(value) else rep_len(as.character(name), length(value))
  sets <- lapply(seq_along(value), function(i) {
    v <- value[i]
    s <- spread[i]
    switch(shape,
           triangular = fs_triangular(v - s, v, v + s, name = nm[i]),
           trapezoidal = fs_trapezoidal(v - s, v - core[i], v + core[i], v + s,
                                        name = nm[i]),
           gaussian = fs_gaussian(v, sd = s / 3, lower = v - s, upper = v + s,
                                  name = nm[i]))
  })
  if (length(sets) == 1L) sets[[1L]] else sets
}
