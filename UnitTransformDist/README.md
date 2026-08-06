# unitTransformDist

**Generalized (0,1)-interval distributions from *any* base distribution on
(0, ∞), via the odds-inverting transform `Y = 1/(1+X)`.**

Unlike registry-restricted designs, the base distribution is supplied
directly as density/CDF/quantile/random-generation functions — so this
works with built-in R distributions (`dgamma`/`pgamma`/...), or with a
fully custom research distribution that has never been registered
anywhere, even one with no closed-form quantile function.

Version 0.1.0 also ships a corrected, distribution-agnostic
profile-likelihood confidence interval, Monte Carlo utilities for
studying estimator bias and coverage, and several bias-correction
methods — see [Key contribution](#key-contribution) below for the
empirical numbers behind that fix.

---

## Table of contents

- [Why this package](#why-this-package)
- [Key contribution](#key-contribution)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Package architecture](#package-architecture)
- [Design principles](#design-principles)
- [Documentation and papers](#documentation-and-papers)
- [Testing](#testing)
- [Status and known limitations](#status-and-known-limitations)
- [Related R packages](#related-r-packages)
- [Citing this package](#citing-this-package)
- [Contributing](#contributing)
- [Authors](#authors)
- [License](#license)

---

## Why this package

Most unit-interval distribution packages in R — [`betareg`](https://CRAN.R-project.org/package=betareg),
[`simplexreg`](https://CRAN.R-project.org/package=simplexreg),
[`gkwdist`](https://CRAN.R-project.org/package=gkwdist),
[`logitnorm`](https://CRAN.R-project.org/package=logitnorm) — are built
around **one specific parametric family**, fixed in advance by the
package authors. `unitTransformDist` instead takes the base
distribution as a plain argument:

```r
d <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5))
```

Swap in `dweibull`/`pweibull`, `dlnorm`/`plnorm`, or a density you wrote
yourself for a distribution that has never been implemented in R before
— nothing else in the package changes. `qdist` and `rdist` are optional;
if omitted, quantiles are obtained by numerically inverting `pdist`, and
random draws follow by inverse-CDF sampling. This is demonstrated in the
package vignette with the XGamma distribution of Sen, Maiti & Chandra
(2016), which has no closed-form CDF or quantile function and is not
implemented in any other R package.

## Key contribution

Profile-likelihood confidence intervals require **re-optimizing nuisance
parameters at every candidate value** of the parameter of interest — not
freezing them at their unconstrained MLE. Freezing them is an easy
mistake to make (we made it ourselves in an earlier prototype of this
idea) and produces intervals with no guaranteed coverage.

`profile_likelihood_ci()` does the re-optimization correctly.
`safe_profile_likelihood()` is the same fix factored out as a
standalone, distribution-agnostic routine that works on **any**
negative log-likelihood, not just objects from this package.
`compare_profile_likelihood_coverage()` lets you check the difference
yourself:

```r
true_dist <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5, rate = 1))
compare_profile_likelihood_coverage(
  true_dist, c(shape = 2.5, rate = 1), dgamma,
  start = c(shape = 1, rate = 1), param_name = "shape",
  n = 60, n_sim = 150, lower = c(shape = 1e-6, rate = 1e-6), seed = 1)
#>                  method coverage_95 n_valid
#> 1  naive_fixed_nuisance       0.587     150
#> 2 corrected_reoptimized       0.940     150
```

The naive (broken) construction's actual 95% coverage is **~59%**; the
corrected construction's is **~94%**. We repeated this across sample
sizes n = 30, 60, 100, 200 — the naive method's shortfall does **not**
shrink as n grows, and its severity depends on how correlated the
profiled and nuisance parameters are (small for Weibull and log-normal
parameterizations we tested, large for gamma's shape/rate). Full
details, proofs, and reproducible code are in the package's JSS and R
Journal papers (see [Documentation and papers](#documentation-and-papers)).

## Installation

Not yet on CRAN (see [Status](#status-and-known-limitations)). Install
the development version from source:

```r
# from a local clone or extracted tarball
install.packages("path/to/unitTransformDist", repos = NULL, type = "source")

# or, once a GitHub remote is public:
# remotes::install_github("yamenetoo/unitTransformDist")
```

Only base R (`stats`, `graphics`) is required. `ggplot2` is an optional
`Suggests` dependency used only if you want `ggplot2`-based plots
instead of the base-R fallbacks.

## Quick start

```r
library(unitTransformDist)

# Any base distribution works, built-in or custom:
d <- utd(dgamma, pgamma, qgamma, rgamma, params = list(shape = 2.5))
dutd(0.3, d)
y <- rutd(500, d)

fit <- fit_utd_mle(y, dgamma, start = c(shape = 1), lower = c(shape = 1e-6))
summary(fit)

# Corrected profile-likelihood CI: re-optimizes nuisance parameters at
# every grid point (see NEWS.md for the historical bug this fixes)
profile_likelihood_ci(fit, "shape")

# Monte Carlo study of the estimator's bias / RMSE
mc <- run_monte_carlo(d, true_params = c(shape = 2.5), ddist = dgamma,
                       start = c(shape = 1), lower = c(shape = 1e-6),
                       n = 100, n_sim = 1000, seed = 1)
summary(mc)

# Bias correction
bias_correct(fit, y, method = "jackknife")

# Diagnostics and goodness of fit
qq_plot(fit, pgamma)
ks_test(fit, pgamma)
```

A fully custom distribution with no closed-form quantile function works
identically, with no package code to modify:

```r
dcustom <- function(x, rate, log = FALSE) dexp(x, rate = rate, log = log)
pcustom <- function(q, rate, lower.tail = TRUE) pexp(q, rate, lower.tail = lower.tail)
d2 <- utd(dcustom, pcustom, params = list(rate = 1.5))
qutd(0.5, d2)  # obtained by numerically inverting pcustom
```

See `vignette("getting-started", package = "unitTransformDist")` for a
complete walkthrough, including the XGamma case study.

## Package architecture

| File | Classes | Key exported functions |
|---|---|---|
| `core.R` | `utd`, `utd_fit` | `utd()`, `dutd()`, `putd()`, `qutd()`, `rutd()`, `fit_utd_mle()`, `profile_likelihood_ci()` |
| `registry.R` *(optional convenience, nothing else depends on it)* | — | `utd_register()`, `utd_get()`, `utd_list_registered()`, `utd_preset_common()` |
| `montecarlo.R` | `utd_mc_sim` | `run_monte_carlo()`, `compare_estimators_mc()`, `mc_sample_size_study()`, `plot.utd_mc_sim()` |
| `bias_correction.R` | — (uses `utd_fit`) | `bias_correct_jackknife()`, `bias_correct_bootstrap()`, `bias_correct_median()`, `bias_correct_firth()`, `bias_correct_double_bootstrap()`, `bias_correct()` |
| `contribution.R` *(flagship — see [Key contribution](#key-contribution))* | — | `safe_profile_likelihood()`, `compare_profile_likelihood_coverage()`, `bias_correct_analytic()` **[experimental, not verified]** |
| `diagnostics.R` | — (methods on `utd_fit`) | `pearson_residuals()`, `quantile_residuals()`, `cox_snell_residuals()`, `cooks_distance()` |
| `extras.R` | — (S3 generics) | `AIC.utd_fit()`, `BIC.utd_fit()`, `wald_ci()`, `bootstrap_ci()`, `ks_test()`, `hazard()` |
| `plot.R` | — (S3 plot methods) | `plot_pdf()`, `plot_cdf()`, `qq_plot()` |

Every function above is exported, documented with runnable `@examples`,
and covered by at least one `testthat` test. A clustered diagram of this
same table is reproduced as a figure in both accompanying papers.

## Design principles

These are stated in full at the top of `R/core.R` and followed
throughout the codebase:

1. **Mathematical contract before code** — every transform is derived on
   paper first; boundary behaviour at `y → 0⁺` and `y → 1⁻` is checked
   explicitly.
2. **Input validation at the exported boundary** — user-facing functions
   validate domain and type and raise informative errors; internal
   helpers (prefixed `.`) skip redundant checks.
3. **Numerical stability** — likelihoods are computed on the log scale
   throughout; densities return `-Inf`/`0` rather than `NaN` outside the
   support.
4. **Reproducibility** — no function seeds randomness globally; callers
   control seeds explicitly via a `seed` argument.
5. **Honest failure** — every `stats::optim()` call inspects
   `$convergence` and raises a `warning()` rather than silently
   returning a misleading fit.
6. **Base-R idiom** — `d`/`p`/`q`/`r` naming mirrors base R's `stats`
   distributions; S3 generics (`print`, `summary`, `coef`, `logLik`,
   `vcov`, `AIC`, `BIC`) are implemented instead of ad hoc accessors,
   and `summary()` methods return a `summary.*` object rather than
   printing directly (matching `summary.lm()`'s convention).

## Documentation and papers

- **Vignette**: `vignette("getting-started", package = "unitTransformDist")`
  — full walkthrough including the corrected profile-likelihood CI and
  the XGamma case study.
- **Function reference**: `?utd`, `?fit_utd_mle`, `?profile_likelihood_ci`,
  or browse `help(package = "unitTransformDist")`.
- **Papers** (in `paper/`): a 20-page *Journal of Statistical Software*
  submission (`paper/jss/`) with formal propositions and proofs for the
  transform's density/CDF/quantile/hazard relationships and a
  score-invariance result underlying the whole package, and a 20-page
  *R Journal* submission (`paper/rjournal/`) focused on the software
  design and the profile-likelihood coverage study. Both compile
  cleanly with `pdflatex`; every number in both papers is reproducible
  from package code with the seeds given in the text.

## Testing

```r
# from an R session in the package root
devtools::test()

# or, from the shell, after building
R CMD check --as-cran unitTransformDist_0.1.0.tar.gz
```

`tests/testthat/test-contribution.R` specifically regression-tests the
profile-likelihood fix: it asserts that the corrected construction
achieves higher coverage than the naive one and exceeds 85% in a fresh
Monte Carlo run every time the test suite executes.

## Status and known limitations

Version 0.1.0, pre-CRAN. Stated plainly rather than hidden:

- `R CMD check --as-cran` has been exercised by direct execution of
  every example and test rather than the full automated check in a
  network-connected environment; run it yourself before relying on this
  for a CRAN submission.
- `NAMESPACE` was written by hand to mirror what `roxygen2::roxygenise()`
  would produce; regenerate it with roxygen2 if you have network access
  to CRAN.
- `bias_correct_analytic()` is explicitly experimental: its
  finite-difference third-derivative implementation produced an
  implausibly small correction in initial testing and has not been
  validated against a closed-form case. Use `bias_correct_jackknife()`
  or `bias_correct_bootstrap()` instead.
- No regression, mixture, or censored-data extension of the transform
  yet exists (sketched as proposed APIs in the papers' "Future work"
  sections, not implemented).
- None of the Monte Carlo/bootstrap loops are parallelized internally.

## Related R packages

If your data fits one of these fixed families, they will generally be
faster (several have compiled backends) and better-documented for that
specific case:

- [`betareg`](https://CRAN.R-project.org/package=betareg) — beta regression
- [`simplexreg`](https://CRAN.R-project.org/package=simplexreg) — simplex distribution regression
- [`gkwdist`](https://CRAN.R-project.org/package=gkwdist) — generalized Kumaraswamy family, C++ backend
- [`logitnorm`](https://CRAN.R-project.org/package=logitnorm) — logit-normal distribution
- [`Newdistns`](https://CRAN.R-project.org/package=Newdistns) — the closest existing precedent for
  "arbitrary base distribution" generality, using a different (generator-family) transform design

See the papers in `paper/` for a detailed comparison, including exactly
what `unitTransformDist` does and does not add relative to each.

## Citing this package

A `CITATION` file and formal DOI will be added once the package is
released on CRAN. Until then, please cite the accompanying JSS or R
Journal manuscript in `paper/`.

## Contributing

Issues and pull requests are welcome once the repository has a public
remote. Please include a minimal reproducible example and, for anything
touching `core.R` or `contribution.R`, a `testthat` test demonstrating
the behaviour you're fixing or adding.

## Authors

**Mohamad Yamen Al Mohamad** — author, maintainer
Department of Statistics, Faculty of Mathematics, Statistics, and Computer Sciences, University of Tabriz, Tabriz, Iran
Email: [yamenmohamad@tabrizu.ac.ir](mailto:yamenmohamad@tabrizu.ac.ir)

**Ali Akbar Heydari** — author
Department of Statistics, Faculty of Mathematics, Statistics, and Computer Sciences, University of Tabriz, Tabriz, Iran
Email: [heydari@tabrizu.ac.ir](mailto:heydari@tabrizu.ac.ir)

**Hossein Jabbari Khamnei** — author
Department of Statistics, Faculty of Mathematics, Statistics, and Computer Sciences, University of Tabriz, Tabriz, Iran
Email: [h_jabbari@tabrizu.ac.ir](mailto:h_jabbari@tabrizu.ac.ir)

## License

MIT — see [`LICENSE.md`](LICENSE.md).
