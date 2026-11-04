# bootqcc 0.1.0

* First release.
* `bootqcc()` builds Phase I control limits for a quantile (or any
  user-defined statistic) of an arbitrary parametric distribution by
  parametric bootstrap with maximum likelihood estimation, and estimates the
  in-control average run length (ARL) by Monte Carlo simulation.
* `arl()` re-evaluates the run length of a chart, for example under a
  parameter shift.
* `bootqcc_monitor()` applies a chart to new samples (Phase II).
* `bootqcc_grid()` runs a simulation study over a grid of settings.
