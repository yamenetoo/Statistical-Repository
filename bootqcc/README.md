# bootqcc

Parametric bootstrap control charts for **any** distribution.

`bootqcc` builds Phase I control limits for a quantile (or any statistic) of a
parametric distribution by parametric bootstrap with maximum likelihood
estimation, then evaluates the chart through the average run length (ARL). It
generalizes the Lindley-geometric code in
[`Statistical-Repository/LGBootStrapControlChart`](https://github.com/yamenetoo/Statistical-Repository/tree/main/LGBootStrapControlChart)
to any family that has density, quantile and random-generation functions.

## Installation

```r
# install.packages("remotes")
remotes::install_github("yamenetoo/Statistical-Repository", subdir = "bootqcc")
```

## Quick start: a built-in distribution

```r
library(bootqcc)

ch <- bootqcc("weibull", params = list(shape = 2, scale = 1),
              u = 0.1,        # monitor the 10% quantile
              n = 5,          # subgroup size
              alpha = 0.0027, # false-alarm probability for the limits
              B = 1000, M = 30, arl_rep = 100, seed = 1)
ch
```

Any name works as long as `d<name>`, `q<name>` and `r<name>` exist: `"gamma"`,
`"lnorm"`, `"exp"`, `"weibull"`, `"logis"`, and so on, including those from
other packages once they are loaded.

## Your own distribution

Pass a list with the functions `d`, `q`, `r` (and an optional `name`). Each
function takes the data, the sample size or the probabilities first, then the
parameters **by name**.

```r
kw <- list(
  name = "kumaraswamy",
  d = function(x, a, b) a * b * x^(a - 1) * (1 - x^a)^(b - 1),
  q = function(p, a, b) (1 - (1 - p)^(1 / b))^(1 / a),
  r = function(n, a, b) (1 - (1 - runif(n))^(1 / b))^(1 / a)
)

ch <- bootqcc(kw, params = list(a = 2, b = 3), u = 0.5, n = 10, seed = 1)
```

For the Lindley-geometric case of the original code, load the
`LindleyPowerSeries` package and wrap its `dlindleygeometric`,
`qlindleygeometric` and `rlindleygeometric` in such a list (or define them with
those names in your session and call `bootqcc("lindleygeometric", ...)`).

## Out-of-control behavior and Phase II

```r
# ARL after a scale shift
arl(ch, params = list(shape = 2, scale = 1.5), rep = 200)

# Monitor new samples (one sample per row) and plot
new <- matrix(rweibull(50, 2, 1.6), ncol = 5)
bootqcc_monitor(ch, new)
plot(ch, data = new)
```

## Other statistics and estimators

```r
# Monitor the mean instead of a quantile
bootqcc("weibull", list(shape = 2, scale = 1),
        statistic = function(est) est$scale * gamma(1 + 1 / est$shape))

# Use your own estimator (must return a named list, or NULL on failure)
bootqcc("exp", list(rate = 1), u = 0.5,
        fit_fn = function(x) list(rate = 1 / mean(x)))
```

## Simulation study over a grid

```r
g <- expand.grid(shape = c(1.5, 2), scale = 1)
res <- bootqcc_grid("weibull", g, u_values = c(0.01, 0.1, 0.9, 0.99),
                    n_values = c(4, 5, 6), alpha = 0.07,
                    B = 1000, M = 10, arl_rep = 500)
write.csv(res, "res.csv", row.names = FALSE)
```

## Main functions

| Function | Purpose |
|---|---|
| `bootqcc()` | Build a chart: bootstrap limits plus in-control ARL |
| `arl()` | Re-estimate the ARL, in control or under shifted parameters |
| `bootqcc_monitor()` | Apply the chart to new samples |
| `bootqcc_grid()` | Simulation study over sample sizes, quantiles and parameters |
| `print()`, `plot()` | Summaries and control chart plots |

## Notes

- Fits that fail to converge are discarded and redrawn; the number of attempts
  is bounded (`max_fit_tries`) and an informative error is raised otherwise.
- Run lengths are capped at `max_run_length`; capped runs are reported as
  `censored` so an unreliable ARL is never silent.
- Limits use the empirical `alpha/2` and `1 - alpha/2` quantiles of the
  bootstrap statistics, averaged over `M` Monte Carlo iterations.
- Starting values default to the true parameters. Supply `start` (and `lower`,
  `upper` with `optim_method = "L-BFGS-B"`) for hard-to-fit families.

## License

MIT
