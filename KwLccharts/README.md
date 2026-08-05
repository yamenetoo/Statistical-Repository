# KwLccharts

Parametric bootstrap percentile control charts for the **Kumaraswamy–Lindley
(KwL) distribution** — an R package accompanying the manuscript *"A
Parametric Bootstrap Control Chart for Kumaraswamy–Lindley Percentiles"*.

The methodology follows Al-Lami, Jabbari Khamnei & Heydari (2025, *PLOS ONE*
20(2): e0316449), *"A parametric bootstrap control chart for Lindley
Geometric percentiles,"* adapted from the two-parameter Lindley-Geometric
distribution to the three-parameter Kumaraswamy–Lindley distribution
obtained via the Kumaraswamy-G construction (Cordeiro & de Castro, 2011).

## Installation

```r
install.packages("KwLccharts_0.1.0.tar.gz", repos = NULL, type = "source")
```

## What's included

| Function | Purpose |
|---|---|
| `dKwL`, `pKwL`, `qKwL`, `rKwL` | Density, cdf, closed-form quantile (via the Lambert W function), and random generation for the KwL distribution |
| `gLindley`, `GLindley`, `QLindley` | Density, cdf and quantile of the Lindley baseline |
| `lambertW0`, `lambertWm1` | Self-contained (no external dependency) Lambert W function, both real branches |
| `fit_KwL` | Maximum-likelihood estimation of (θ, a, b), multi-start L-BFGS-B |
| `KwLBootCC` | Parametric bootstrap percentile control chart (Phase-I control limits) |
| `KwL_ARL` | Phase-II average-run-length (ARL) simulation, in-control or shifted |

## Quick example

```r
library(KwLccharts)
set.seed(1)

## KwL(theta = 1, a = 1.5, b = 1.5), 5th-percentile bootstrap control chart
lim <- KwLBootCC(u = 0.05, alpha = 0.0027, theta = 1, a = 1.5, b = 1.5,
                  n = 5, B = 300, MC = 10)
lim$MLCL; lim$MUCL

## Phase-II average run length, in control
KwL_ARL(u = 0.05, LCL = lim$MLCL, UCL = lim$MUCL,
        theta1 = 1, a1 = 1.5, b1 = 1.5, n_runs = 20, max_run = 700)
```

## Reproducing the manuscript results

The `scripts/` directory contains the exact R scripts used to produce every
number reported in the manuscript (control-limit table, ARL table, and the
real-data example on the gastric-cancer survival-time data set of Bekker,
Roux & Mosteit, 2000). Because three-parameter MLE refitting inside a
bootstrap loop is computationally demanding, the settings used in the
manuscript (B = 300 bootstrap iterations, MC = 10 outer repetitions, 20
Phase-II runs per ARL entry) are modest; increase `B`, `MC` and `n_runs` for
production-quality precision.

## Validation performed

- The KwL density is confirmed to integrate to 1 (numerical integration).
- `qKwL` is confirmed to invert `pKwL` to within 1e-6 across the support.
- `fit_KwL` is confirmed to be consistent: recovered parameters converge to
  the true generating values as the sample size grows (verified at
  n = 1,000 / 5,000 / 20,000).
- `tests/test-all.R` runs these checks automatically (`R CMD check` clean,
  no external test-framework dependency required).

## License

MIT (see `LICENSE`).
