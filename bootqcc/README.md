<div align="center">

<img src="man/figures/logo.png" alt="bootqcc logo" height="220">

# bootqcc

**Parametric bootstrap control charts for any distribution**

[![R package](https://img.shields.io/badge/R-package-276DC3?logo=r&logoColor=white)](#installation)
[![Version](https://img.shields.io/badge/version-0.1.0-5aa9e6)](DESCRIPTION)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![R CMD check](https://img.shields.io/badge/R%20CMD%20check-passing-brightgreen)](#testing)
[![Tests](https://img.shields.io/badge/tests-8%20passing-brightgreen)](tests/testthat)

</div>

---

## Contents

1. [Overview](#overview)
2. [Why bootqcc](#why-bootqcc)
3. [How it works](#how-it-works)
4. [Installation](#installation)
5. [Quick start](#quick-start)
6. [Using your own distribution](#using-your-own-distribution)
7. [Reproducing the Lindley-geometric chart](#reproducing-the-lindley-geometric-chart)
8. [Out-of-control performance and Phase II](#out-of-control-performance-and-phase-ii)
9. [Calibrating the false-alarm rate](#calibrating-the-false-alarm-rate)
10. [Other statistics and estimators](#other-statistics-and-estimators)
11. [Simulation studies over a grid](#simulation-studies-over-a-grid)
12. [Function reference](#function-reference)
13. [Argument reference for `bootqcc()`](#argument-reference-for-bootqcc)
14. [What the result contains](#what-the-result-contains)
15. [Tips for speed and reproducibility](#tips-for-speed-and-reproducibility)
16. [Troubleshooting](#troubleshooting)
17. [Limitations](#limitations)
18. [Roadmap](#roadmap)
19. [Testing](#testing)
20. [Relation to the original code](#relation-to-the-original-code)
21. [Citation](#citation)
22. [Author and license](#author-and-license)

---

## Overview

`bootqcc` builds **Phase I control limits** for a quantile (or any statistic
of your choice) of a parametric distribution using the **parametric
bootstrap** and **maximum likelihood estimation**. It then evaluates the chart
by Monte Carlo estimates of the **average run length (ARL)**, and can monitor
new samples in **Phase II**.

The idea comes from control charts for lifetime and reliability data, where
the quantity of interest is a *percentile* of a skewed distribution, not a
mean, and where a classical Shewhart chart is a poor fit. The original code in
this repository solved the problem for one family, the Lindley-geometric
distribution. `bootqcc` generalizes it: **any family with density, quantile and
random-generation functions works**, whether it is built into R or defined by
you.

## Why bootqcc

| You have | Without bootqcc | With bootqcc |
|---|---|---|
| A new lifetime distribution | Rewrite the whole chart code for it | Pass its `d`, `q`, `r` functions |
| A built-in family (Weibull, gamma, log-normal, ...) | Write and debug a bootstrap loop | `bootqcc("weibull", ...)` |
| A statistic other than a quantile | Modify the internals | `statistic = function(est) ...` |
| A special estimator | Modify the internals | `fit_fn = function(x) ...` |
| A simulation study for a paper | Hand-written nested loops | `bootqcc_grid()` |
| Samples to monitor | Manual bookkeeping | `bootqcc_monitor()` and `plot()` |

It is also defensive: failed fits are redrawn (within a bound, with a clear
error), run lengths are capped, and capped runs are *reported* so an unreliable
ARL is never silent.

## How it works

Let `T(θ̂)` be the monitored statistic computed from fitted parameters `θ̂`.
By default `T` is the `u`-quantile of the fitted distribution, `Q(u; θ̂)`.

For each Monte Carlo iteration `m = 1, ..., M`:

1. Draw `B` samples of size `n` from the distribution with the in-control
   parameters `θ`.
2. Fit the distribution to each sample by maximum likelihood, giving
   `θ̂₁, ..., θ̂_B`.
3. Compute `T₁, ..., T_B` with `T_b = T(θ̂_b)`.
4. Record the empirical `α/2` and `1 − α/2` quantiles of the `T_b` as one lower
   and one upper limit, `LCL_m` and `UCL_m`.

The reported control limits are the averages over the `M` iterations,

```
LCL = (1/M) Σ LCL_m        UCL = (1/M) Σ UCL_m
```

and their standard deviations over `m` tell you how stable the limits are.
Averaging over `M` iterations reduces the Monte Carlo noise in the limits.

**Phase II / ARL.** Samples of size `n` are drawn repeatedly, the distribution
is fitted to each, and `T(θ̂)` is compared with `[LCL, UCL]`. The run length is
the number of samples up to and including the first signal. Repeating this
`arl_rep` times gives the ARL and its standard deviation. Drawing from shifted
parameters gives the out-of-control ARL.

## Installation

From GitHub (the package lives in a subfolder of this repository):

```r
# install.packages("remotes")
remotes::install_github("yamenetoo/Statistical-Repository", subdir = "bootqcc")
```

From a downloaded source archive:

```r
install.packages("bootqcc_0.1.0.tar.gz", repos = NULL, type = "source")
```

The package depends only on base R (`stats` and `graphics`). No compiled code
is needed.

## Quick start

```r
library(bootqcc)

ch <- bootqcc("weibull",
              params = list(shape = 2, scale = 1),
              u = 0.1,         # monitor the 10% quantile
              n = 5,           # subgroup size
              alpha = 0.0027,  # false-alarm probability used for the limits
              B = 500, M = 10, arl_rep = 100, seed = 1)
ch
```

Output (real run):

```
Parametric bootstrap control chart
Distribution : weibull (shape = 2, scale = 1)
Statistic    : 0.1-quantile, subgroup size n = 5
Settings     : alpha = 0.0027, B = 500, M = 10
LCL = 0.06954 (SD 0.0156)
UCL = 1.118 (SD 0.0897)
In-control ARL = 175.8 (SD 181.2) from 100 runs
```

Any family works as long as `d<name>`, `q<name>` and `r<name>` exist in your
session: `"weibull"`, `"gamma"`, `"lnorm"`, `"exp"`, `"logis"`, `"norm"`, and
the families of any package you have loaded.

To look at the bootstrap distribution of the statistic with the limits:

```r
plot(ch)
```

## Using your own distribution

Pass a list with the functions `d`, `q` and `r` (and an optional `name`). Each
function takes the data, the sample size, or the probabilities as its first
argument, followed by the parameters **by name**:

```r
kw <- list(
  name = "kumaraswamy",
  d = function(x, a, b) a * b * x^(a - 1) * (1 - x^a)^(b - 1),
  q = function(p, a, b) (1 - (1 - p)^(1 / b))^(1 / a),
  r = function(n, a, b) (1 - (1 - runif(n))^(1 / b))^(1 / a)
)

ch <- bootqcc(kw, params = list(a = 2, b = 3), u = 0.5, n = 15, seed = 1)
```

Alternatively, define functions named `d<name>`, `q<name>` and `r<name>` in your
session and refer to the family by `name`, exactly like the built-in ones:

```r
dmydist <- function(x, lambda, k) { ... }
qmydist <- function(p, lambda, k) { ... }
rmydist <- function(n, lambda, k) { ... }

bootqcc("mydist", list(lambda = 1, k = 2), u = 0.99, n = 5)
```

Requirements for the functions:

- `d` must return a positive density at points inside the support and may
  return `0`, `NaN` or an error outside it (these are treated as invalid
  parameter values by the optimizer);
- `q` and `r` must accept the same named parameters as `d`;
- parameter names in `params` (and `start`) must match the argument names.

## Reproducing the Lindley-geometric chart

The original code in `LGBootStrapControlChart` used the
[`LindleyPowerSeries`](https://cran.r-project.org/package=LindleyPowerSeries)
package. The same chart with `bootqcc`:

```r
library(LindleyPowerSeries)

lg <- list(
  name = "lindleygeometric",
  d = function(x, theta, p) dlindleygeometric(x, theta, p),
  q = function(u, theta, p) qlindleygeometric(u, theta, p),
  r = function(n, theta, p) rlindleygeometric(n, theta, p)
)

ch <- bootqcc(lg, params = list(theta = 2, p = 0.8),
              u = 0.99, alpha = 0.07, n = 5,
              B = 10000, M = 20, arl_rep = 1000,
              start = list(theta = 0.5, p = 0.01))
```

> This snippet mirrors the call in the original `example.R`. It was not run
> while preparing this README because `LindleyPowerSeries` was not available
> in the build environment; verify the numbers against the original script
> before using them in a publication.

## Out-of-control performance and Phase II

**ARL after a shift.** Simulate from shifted parameters to estimate the
out-of-control ARL:

```r
for (sc in c(1, 1.2, 1.5, 2)) {
  a <- arl(ch, params = list(shape = 2, scale = sc), rep = 100)
  cat(sprintf("scale = %.1f  ARL = %.1f\n", sc, a$mean))
}
```

Output (real run, with the chart above):

```
scale = 1.0  ARL = 178.9
scale = 1.2  ARL = 78.3
scale = 1.5  ARL = 20.3
scale = 2.0  ARL = 4.5
```

**Monitoring new samples.** One sample per row (a matrix or data frame), or a
list of numeric vectors:

```r
new <- rbind(matrix(rweibull(15, 2, 1), ncol = 5),   # in control
             matrix(rweibull(10, 2, 2), ncol = 5))   # scale doubled
bootqcc_monitor(ch, new)
plot(ch, data = new)       # draws the Phase II control chart
```

```
  sample statistic        LCL      UCL         status
1      1 0.6799767 0.06953615 1.118293     in control
2      2 0.4940795 0.06953615 1.118293     in control
3      3 0.4854396 0.06953615 1.118293     in control
4      4 0.8747797 0.06953615 1.118293     in control
5      5 1.2695785 0.06953615 1.118293 out of control
```

Samples for which no fit could be obtained are labelled `"fit failed"`
instead of being silently dropped.

## Calibrating the false-alarm rate

The nominal `alpha` and the realized in-control ARL are **not** tied by the
simple formula `ARL0 = 1/alpha`. In the Quick start example `alpha = 0.0027`
would suggest an ARL of about 370, but the simulated in-control ARL was about
176. The reason is that with a small subgroup size (`n = 5`) the fitted
parameters are variable, and the fitted statistic is compared with limits that
were built from the *true* parameters.

So treat `alpha` as a tuning parameter: choose it so that the simulated
in-control ARL matches your target.

```r
for (a in c(0.0027, 0.0015, 0.001, 0.0005)) {
  ch_a <- bootqcc("weibull", list(shape = 2, scale = 1), u = 0.1, n = 5,
                  alpha = a, B = 500, M = 10, arl_rep = 300, seed = 1)
  cat(sprintf("alpha = %.4f  ARL0 = %.0f\n", a, ch_a$arl$mean))
}
```

Increase `arl_rep` for a stable estimate: the ARL has a large standard
deviation (close to its mean for a geometric-like run length), so a few
hundred runs are needed for a useful standard error.

## Other statistics and estimators

**Monitor a different statistic.** Supply a function of the fitted parameter
list. Here, the Weibull mean:

```r
bootqcc("weibull", list(shape = 2, scale = 1),
        statistic = function(est) est$scale * gamma(1 + 1 / est$shape),
        n = 8, alpha = 0.01)
```

Other ideas: reliability at a mission time `1 - p<name>(t, ...)`, hazard
rate at a time point, or a tail probability.

**Use your own estimator.** `fit_fn` takes a sample and returns a named list
of estimates (or `NULL` on failure). Closed-form estimators are much faster
than numerical optimization:

```r
bootqcc("exp", list(rate = 1), u = 0.5, n = 10,
        fit_fn = function(x) list(rate = 1 / mean(x)))
```

**Help the optimizer.** For awkward families supply starting values and
bounds:

```r
bootqcc("gamma", list(shape = 3, scale = 2), u = 0.9, n = 10,
        start = list(shape = 2, scale = 1),
        optim_method = "L-BFGS-B",
        lower = c(1e-6, 1e-6), upper = c(Inf, Inf))
```

## Simulation studies over a grid

`bootqcc_grid()` runs `bootqcc()` for every combination of sample size,
quantile level and parameter vector, and returns one table. This replaces the
nested loops of a typical Monte Carlo study:

```r
g <- expand.grid(shape = c(1.5, 2), scale = 1)
res <- bootqcc_grid("weibull", g,
                    u_values = c(0.1, 0.9), n_values = 5,
                    alpha = 0.05, B = 200, M = 5, arl_rep = 50, seed = 4)
res
```

```
  n   u shape scale   MLCL  MUCL   SDLCL  SDUCL MARL SDARL censored
1 5 0.1   1.5     1 0.0654 0.812 0.00808 0.0389 15.2  16.9        0
2 5 0.1   2.0     1 0.1292 0.855 0.01189 0.0308 15.2  16.9        0
3 5 0.9   1.5     1 0.7817 2.680 0.03242 0.1104 14.1  13.0        0
4 5 0.9   2.0     1 0.8312 2.095 0.02591 0.0645 14.1  13.0        0
```

The columns are the mean and standard deviation of the lower and upper limits
(`MLCL`, `MUCL`, `SDLCL`, `SDUCL`), the mean and standard deviation of the ARL
(`MARL`, `SDARL`), and the number of censored runs. Because `seed` is passed to
every call, the run lengths can coincide across rows; omit `seed` to get
independent simulations per row.

Save the result for your paper:

```r
write.csv(res, "res.csv", row.names = FALSE)
```

## Function reference

| Function | Purpose |
|---|---|
| `bootqcc()` | Build a chart: bootstrap limits plus in-control ARL |
| `arl()` | Re-estimate the ARL of an existing chart, in control or under shifted parameters |
| `bootqcc_monitor()` | Apply a chart to new samples (Phase II) |
| `bootqcc_grid()` | Simulation study over sample sizes, quantiles and parameters |
| `print()` | Summary of a chart |
| `plot()` | Bootstrap distribution with limits, or the Phase II chart when `data` is given |

Help pages: `?bootqcc`, `?arl`, `?bootqcc_monitor`, `?bootqcc_grid`.

## Argument reference for `bootqcc()`

| Argument | Default | Meaning |
|---|---|---|
| `distr` | | Family name (`"weibull"`) or list `list(d, q, r, name)` |
| `params` | | Named list of the in-control parameters |
| `u` | `0.5` | Probability level of the monitored quantile |
| `n` | `5` | Subgroup (sample) size |
| `alpha` | `0.0027` | Total false-alarm probability for the limits |
| `B` | `1000` | Bootstrap replicates per Monte Carlo iteration |
| `M` | `30` | Monte Carlo iterations for the limits |
| `arl_rep` | `100` | Simulated run lengths for the ARL |
| `statistic` | `NULL` | Function of the fitted parameter list, replaces the quantile |
| `start` | `params` | Starting values for the maximum likelihood fit |
| `fit_fn` | `NULL` | Custom estimator, replaces the built-in maximum likelihood fit |
| `optim_method` | `"Nelder-Mead"` | Method passed to `stats::optim()` |
| `lower`, `upper` | `-Inf`, `Inf` | Bounds for `"L-BFGS-B"` or `"Brent"` |
| `max_run_length` | `1e5` | Cap on a single run length |
| `max_fit_tries` | `20` | Failed-fit budget: at most `max_fit_tries * B` attempts per iteration |
| `seed` | `NULL` | Random seed |
| `verbose` | `FALSE` | Print progress |
| `env` | caller | Environment where `d<name>`, `q<name>`, `r<name>` are looked up |

## What the result contains

`bootqcc()` returns an object of class `"bootqcc"`, a list with:

| Element | Content |
|---|---|
| `limits` | Named vector with `LCL` and `UCL` (means over the `M` iterations) |
| `sd_limits` | Standard deviation of the `M` lower and upper limits |
| `lcl_runs`, `ucl_runs` | The `M` individual limits |
| `arl` | List with `mean`, `sd`, `runs` and `censored` |
| `boot_stats` | Bootstrap statistics of the last Monte Carlo iteration |
| `params`, `u`, `n`, `alpha`, `B`, `M` | Settings used |
| `distr`, `distr_name` | The resolved distribution functions and its name |

```r
ch$limits
ch$arl$mean
hist(ch$arl$runs)
```

## Tips for speed and reproducibility

- **Cost.** One chart needs about `M * B` fits for the limits plus the fits of
  the ARL runs. Reduce `B` and `M` while exploring (`B = 200`, `M = 5`), then
  increase them for the final numbers.
- **Closed-form estimators** passed through `fit_fn` can be orders of
  magnitude faster than numerical maximum likelihood.
- **Reproducibility.** Set `seed`. For parallel studies, run separate
  `bootqcc_grid()` calls in separate R sessions with different seeds.
- **Stable limits.** Check `sd_limits`: if the limits vary a lot between Monte
  Carlo iterations, increase `B`.
- **Stable ARL.** Check `arl$censored` and the spread of `arl$runs`; increase
  `arl_rep` for a tighter standard error.

## Troubleshooting

| Message or symptom | Likely cause and fix |
|---|---|
| `Could not find function 'dxxx'` | The `d`, `q` or `r` function for that name is not defined or not loaded. Define it, load its package, or pass a list. |
| `'params' must be a fully named list` | Use `list(shape = 2, scale = 1)`, not `c(2, 1)`. |
| `Too many failed fits while building the limits` | Parameters or starting values do not suit the family, or `n` is too small for the model. Provide `start`, try `L-BFGS-B` with bounds, or increase `n`. |
| ARL reports censored runs | The limits are too wide or the family almost never signals; check `max_run_length`, `alpha` and the parameters. |
| Very different results across seeds | Increase `B`, `M` and `arl_rep`. |
| Argument names do not match | The names in `params` must equal the argument names of your `d`, `q`, `r` functions. |

## Limitations

- **Known in-control parameters.** Phase I draws samples from the parameters
  you give. A function that estimates the parameters from real Phase I data
  and then builds the limits is not included yet (see the roadmap).
- **Single statistic per chart.** One quantile or one custom statistic is
  monitored at a time.
- **Univariate, independent data** with subgroups of fixed size `n`.
- **Speed.** The loops are plain R; large studies (for example `B = 10000` with
  `arl_rep = 1000`) take a long time.
- **Estimator.** The default is Nelder-Mead maximum likelihood, which can be
  fragile for many-parameter or heavily constrained families.
- **Validation.** The package is tested on Weibull, exponential and
  Kumaraswamy models. Results for new families should be validated against
  your own simulations.

## Roadmap

- Estimate the parameters from real Phase I data and build the chart from them
- Parallel execution of the bootstrap and ARL loops
- Parameter transformations for stable fitting of constrained families
- Vignette with a full worked case study
- Additional chart types (individual observations, one-sided charts)
- Preparation for CRAN submission

Suggestions are welcome: open an issue in the repository.

## Testing

```r
devtools::test()      # the testthat suite (8 test blocks)
devtools::check()     # R CMD check
```

The tests cover built-in and user-defined distributions, custom statistics and
estimators, the effect of a shift on the ARL, Phase II monitoring, the grid
function, and input validation.

## Relation to the original code

This package generalizes the Lindley-geometric chart in
[`LGBootStrapControlChart`](../LGBootStrapControlChart) of this repository,
whose README cites Prof. Dr. Ali Akbar Heydari (University of Tabriz,
Department of Statistics). Compared with the original script, `bootqcc`:

- works for any distribution instead of one;
- fixes the swapped assignment of the two fitted parameters;
- always returns the mean ARL (it was undefined when `verbos = FALSE`);
- computes the limits with empirical quantiles instead of an integer index
  that could be zero;
- caps the run length so a chart that never signals cannot loop forever.

## Citation

If you use `bootqcc` in your work, please cite it as:

```
Al-Mohamad, M. Y., Hydari, A. A., and Jabbari Khamnei, H. (2026). bootqcc:
Parametric Bootstrap Control Charts for Any Distribution. R package version 0.1.0.
https://github.com/yamenetoo/Statistical-Repository/tree/main/bootqcc
```

## Author and license

**Mohamad Yamen Al-Mohamad** (corresponding author and maintainer),
Department of Statistics, University of Tabriz
· [ORCID](https://orcid.org/0009-0009-1665-7203)
· [Google Scholar](https://scholar.google.com/citations?user=un_6SQwAAAAJ)
· [yamenmohamad@tabrizu.ac.ir](mailto:yamenmohamad@tabrizu.ac.ir)

**Ali Akbar Hydari**, Department of Statistics, University of Tabriz
· [heydari@tabrizu.ac.ir](mailto:heydari@tabrizu.ac.ir)

**Hossein Jabbari Khamnei**, Department of Statistics, University of Tabriz
· [h_jabbari@tabrizu.ac.ir](mailto:h_jabbari@tabrizu.ac.ir)

Released under the MIT license, see [LICENSE](LICENSE).
