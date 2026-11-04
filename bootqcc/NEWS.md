# bootqcc 0.1.0

* New `bootqcc_phase1()` builds a chart from raw Phase I data: it estimates the
  parameters, builds the bootstrap limits, measures the conditional
  false-alarm probability caused by the estimation with an outer bootstrap and
  can widen the limits (`adjust = "mean"` or `"quantile"`).
* `bootqcc()` now warns when `B * alpha / 2 < 5`, because the tail quantiles
  are then poorly resolved and the in-control ARL is too small.
* New vignette with the validation against exact exponential results and a
  case study; validation script in `inst/validation`.
* Documentation: methodology figure, mathematical formulation and references.

* First release.
* `bootqcc()` builds Phase I control limits for a quantile (or any
  user-defined statistic) of an arbitrary parametric distribution by
  parametric bootstrap with maximum likelihood estimation, and estimates the
  in-control average run length (ARL) by Monte Carlo simulation.
* `arl()` re-evaluates the run length of a chart, for example under a
  parameter shift.
* `bootqcc_monitor()` applies a chart to new samples (Phase II).
* `bootqcc_grid()` runs a simulation study over a grid of settings.
