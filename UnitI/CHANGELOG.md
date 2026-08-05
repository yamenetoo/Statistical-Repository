# Changelog

All notable changes to this project are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [0.3.0] - 2026-08-04

### Added
- `lmoments.py`: sample L-moments, L-moment ratios, and L-moment matching
  estimation (`fit_lmoments`) for arbitrary base distributions.
- `regression.py`: `UnitTransformRegression`, a beta-regression-style model
  linking the scale parameter to a linear predictor via a log or identity link.
- `mixture.py`: `UnitTransformMixture`, finite mixtures fit via EM, with
  `pdf`, `cdf`, `rvs`, `aic`, `bic`.
- `bayes.py`: `map_estimate`, `metropolis_hastings`, `posterior_summary`,
  and `PriorSpec` for dependency-free Bayesian estimation.
- `censored.py`: `fit_censored` for right-, left-, and interval-censored
  data, and for truncated observation windows.
- `diagnostics.py`: Pearson, Cox-Snell, and randomised quantile residuals,
  plus an approximate `cooks_distance` case-deletion influence measure.
- `model_selection.py`: `auto_select_distribution` to rank all candidate
  base distributions by AIC/BIC/log-likelihood.
- `py.typed` marker for type-checker support.
- Restructured to a `src/` layout with a corrected `pyproject.toml`
  (`[tool.setuptools.packages.find]`), fixing a package-discovery failure
  that made `pip install -e .` fail entirely on the previous layout.
- Expanded CI matrix (Linux/macOS/Windows, Python 3.9-3.12) plus a build
  and `twine check` job.

### Fixed
- `np.trapz` (removed in NumPy 2.0) replaced with `np.trapezoid` in tests.
- `extras.summary()` previously logged via `logger.info()` with no handler
  configured, so nothing was printed and `None` was returned, contradicting
  the documented behaviour. It now prints by default and returns the
  formatted report string (`verbose=False` to suppress printing).
- `extras.fit_moments("weibull_min", ...)` called `scipy.stats.gamma(...)`
  (a *distribution* object) instead of the gamma *function*; now uses
  `scipy.special.gamma`.
- `extras.profile_likelihood_ci` referenced `scipy.optimize.brentq` without
  importing the `scipy` top-level namespace (`NameError` on any call that
  reached that branch); `brentq` is now imported directly.
- **`extras.profile_likelihood_ci` produced statistically invalid confidence
  intervals** due to two compounding bugs:
  1. `chdtri(1 - alpha, 1)` passed arguments in the wrong order (the
     correct signature is `chdtri(df, probability)`), collapsing the
     chi-square threshold to ~0.
  2. Nuisance parameters (loc, scale) were fixed at their unconstrained MLE
     values for every candidate value of the parameter of interest instead
     of being re-optimised (true profile likelihood requires maximising
     over nuisance parameters at each point). Together these caused the
     reported 95% CI to collapse to a single point at the MLE and to
     frequently exclude the true parameter value. Both are fixed; CIs are
     now noticeably slower to compute (nested optimisation) but correct.
- `UnitTransformMixture`'s weighted component refit allowed an unbounded
  `loc`, occasionally producing simulated values outside (0,1) from `rvs()`;
  `loc` is now bounded to keep the transformed support non-negative, with a
  defensive clip as a second safety net.

## [0.2.0] - earlier internal version
- Initial `UnitTransformDist` core, registry of 50+ base distributions,
  MLE fitting, AIC/BIC, KS/CvM/AD goodness-of-fit tests, hazard functions,
  bootstrap confidence intervals, method-of-moments, and plotting utilities.
