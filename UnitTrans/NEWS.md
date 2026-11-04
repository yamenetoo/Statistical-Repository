# UnitTrans 0.2.0

## Bug fixes
* The transformation `V = 2 / (1 + exp(Y))` was called `"arcsech"` but equals
  `1 - tanh(Y / 2)`. It is now `"tanh"`; `"arcsech"` still works as a
  deprecated alias and gives a warning.
* `ptrans(lower.tail = FALSE)` and `survival_trans()` no longer compute `1 - p`,
  which lost all precision in the upper tail; the tail is taken directly from
  the baseline. `hazard_trans()` is evaluated on the log scale. `ptrans()` and
  `qtrans()` gain `log.p`.
* `check_domain()` now uses the baseline cdf to measure the probability mass
  outside the admissible domain (e.g. `P(Y < 0)` for `"rational"`) instead of
  inspecting a quantile grid, and `fit_trans()` refuses incompatible
  baseline/transformation pairs (`check_support = FALSE` to override).
* `gof_trans()` p-values use `(1 + #{T* >= T}) / (B + 1)` and are never 0;
  failed bootstrap replications are reported.
* The bootstrap in `confint()` clips simulated data to (0,1) like
  `gof_trans()` and `simulation_study()`, and reports failed replications.
* Functions that take a `seed` no longer change the caller's random number
  stream.

## Improvements
* Optimisation runs on an unconstrained (log) scale with BFGS, a
  Nelder-Mead fallback and a quantile-matching starting point; the Hessian is
  checked for positive definiteness (`hessian_ok`).
* Baseline functions are resolved once per call instead of at every evaluation, and user
  baselines without `log`, `lower.tail` or `log.p` arguments are supported.
* `fit_all()` reports Akaike weights and a `converged` column, drops
  non-converged fits, accepts `seed` and `cores`.
* New: `confint(method = "profile")`, `residuals()`, `predict()`, `nobs()`,
  `fit_trans(adjust = TRUE)` for data containing 0 or 1, `plot(which = )`.
* `summary_trans()` integrates central moments directly (no cancellation).
