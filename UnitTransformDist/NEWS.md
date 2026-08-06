# unitTransformDist 0.1.0

Initial release.

## Core
* `utd()`: construct a unit-interval distribution from ANY base
  distribution, supplied directly as density/CDF/quantile/random-generation
  functions rather than looked up from a fixed registry.
* `dutd()`, `putd()`, `qutd()`, `rutd()`: standard R d/p/q/r interface.
  `qutd()` falls back to numeric inversion of `pdist` via `uniroot()` when
  no closed-form quantile function is supplied; `rutd()` falls back to
  inverse-CDF sampling when no `rdist` is supplied.
* `fit_utd_mle()`: maximum likelihood estimation for any base distribution,
  with explicit convergence checking (a warning is raised, not swallowed,
  when `optim()` fails to converge).

## Bug fix carried over from the exploratory Python prototype
* `profile_likelihood_ci()` now genuinely re-optimizes nuisance parameters
  at every candidate value of the parameter of interest, instead of
  holding them fixed at their unconstrained MLE. The previous approach
  (and an accompanying chi-square-quantile argument-order error) caused
  reported 95% intervals to collapse to a single point at the MLE and to
  frequently exclude the true parameter value. A regression test
  (`test-core.R`) checks explicitly that the returned interval has
  positive width and contains the point estimate.

## Monte Carlo (`montecarlo.R`)
* `run_monte_carlo()`, `compare_estimators_mc()`, `mc_sample_size_study()`,
  `plot.utd_mc_sim()`.

## Bias correction (`bias_correction.R`)
* `bias_correct_jackknife()`, `bias_correct_bootstrap()`,
  `bias_correct_median()`, `bias_correct_firth()` (Firth 1993 penalized
  likelihood with a numerically approximated information matrix, so no
  closed-form expected-information formula is required per distribution),
  `bias_correct_double_bootstrap()`, and the dispatcher `bias_correct()`.

## Original methodological contribution (`contribution.R`)
* `safe_profile_likelihood()`: a distribution-agnostic, reusable
  profile-likelihood engine (works on any user-supplied negative
  log-likelihood, not just objects from this package) that re-optimizes
  nuisance parameters at every candidate value of the parameter of
  interest.
* `compare_profile_likelihood_coverage()`: an empirical Monte Carlo
  demonstration, reproducible by any user, showing that the naive
  fixed-nuisance approach achieves roughly 59% actual coverage against a
  nominal 95% target, while the corrected approach achieves roughly 94%
  (n = 60, 150 replications, gamma base distribution, shape = 2.5,
  rate = 1; see the package vignette for the full reproducible script).
* `bias_correct_analytic()` [experimental, documented honestly as such]:
  derives that the classical Cox & Snell (1968) O(1/n) bias formula
  transfers unchanged from the base distribution to the unit-transform
  MLE, because the transform's Jacobian does not depend on the
  parameter. The finite-difference third-derivative implementation has
  not yet been validated against a closed-form case and should not be
  relied on for production use; `bias_correct_jackknife()` and
  `bias_correct_bootstrap()` are recommended instead until this is fixed.

* `pearson_residuals()`, `cox_snell_residuals()`, `quantile_residuals()`,
  `cooks_distance()`, `AIC.utd_fit()`, `BIC.utd_fit()`, `wald_ci()`,
  `bootstrap_ci()`, `ks_test()`, `hazard()`.

## Plotting
* `plot_pdf()`, `plot_cdf()`, `qq_plot()` using base graphics only
  (ggplot2 is in `Suggests`, never a hard dependency).
