# BKwG 0.1.0

* First version: d/p/q/r, survival, hazard and cumulative hazard functions of
  the beta Kumaraswamy-G family for 12 built-in baselines and user-defined
  baselines; analytical score; moments, generating function, entropies, mean
  residual life, order statistics and mixture representation; estimation by
  MLE, MPS, AD, CvM, LSE and WLSE for the full model and six nested
  sub-models; likelihood ratio and Wald tests; goodness-of-fit statistics with
  parametric bootstrap p-values; model comparison tables.
* A vignette, `vignette("BKwG")`, introduces the package with worked
  examples.
* Fixes found before release: `dbkwg()` no longer returns `Inf` far in the
  right tail for user-defined baselines whose `pG` has no `lower.tail`
  argument; maximum likelihood no longer accepts parameter values where the
  log-density is pure rounding error (which could give a spurious
  log-likelihood of 0 in `compare_bkwg()`); harmless `NaNs produced` warnings
  from the optimiser are no longer shown.
