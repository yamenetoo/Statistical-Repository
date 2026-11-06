# BKwG

**The beta Kumaraswamy-G family of distributions in R — distribution functions for any baseline, six nested sub-models, six estimation methods, likelihood-ratio / Wald / goodness-of-fit tests, and the structural properties of the family.**

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/yamenetoo/Statistical-Repository/blob/main/BKwG/notebooks/BKwG_examples.ipynb)
![version](https://img.shields.io/badge/version-0.1.0-blue)
![license](https://img.shields.io/badge/license-GPL--3-green)
![R](https://img.shields.io/badge/R-%E2%89%A5%203.5-276DC3)

---

## Contents

1. [Overview](#overview)
2. [The BKw-G family](#the-bkw-g-family)
3. [Installation](#installation)
4. [Quick start](#quick-start)
5. [Baseline distributions](#baseline-distributions)
6. [Sub-models](#sub-models)
7. [Distribution functions](#distribution-functions)
8. [Structural properties](#structural-properties)
9. [Parameter estimation](#parameter-estimation)
10. [Hypothesis tests and goodness of fit](#hypothesis-tests-and-goodness-of-fit)
11. [Model comparison](#model-comparison)
12. [User-defined baselines](#user-defined-baselines)
13. [Numerical details](#numerical-details)
14. [Function reference](#function-reference)
15. [Examples notebook](#examples-notebook)
16. [Testing](#testing)
17. [Citation](#citation)
18. [Authors and license](#authors-and-license)

---

## Overview

The beta Kumaraswamy-G (BKw-G) family of Handique, Chakraborty and Ali (2017)
puts the Kumaraswamy-G distribution of Cordeiro and de Castro (2011) inside the
beta-generated construction of Eugene, Lee and Famoye (2002). Starting from
**any** continuous baseline `G`, it adds four shape parameters that control
skewness, tail weight and the shape of the hazard rate (increasing, decreasing,
bathtub, upside-down bathtub).

`BKwG` turns the family into a working tool:

- `d/p/q/r` functions plus survival, hazard and cumulative hazard, accurate in
  the tails, for **12 built-in baselines or any baseline you define**;
- the **full model and six nested sub-models** (Kw-G, Beta-G, EKw-G, BE-G, E-G, G);
- **six estimators**: maximum likelihood, maximum product of spacings,
  minimum Anderson–Darling, minimum Cramér–von Mises, least squares and
  weighted least squares, with any parameter held fixed;
- **inference**: standard errors and confidence intervals, likelihood-ratio and
  Wald tests for the sub-models, W\*, A\* and Kolmogorov–Smirnov statistics with
  parametric-bootstrap p-values, and AIC/CAIC/BIC/HQIC comparison tables;
- **properties**: moments, moment generating function, skewness and kurtosis
  (moment and quantile based), Shannon and Rényi entropies, mean residual life,
  order statistics, the analytical score vector and the infinite-mixture
  representation of the density.

---

## The BKw-G family

Let `G(x)` and `g(x)` be the cdf and pdf of a baseline distribution and

$$K(x) = 1 - \left[1 - G(x)^{a}\right]^{b}$$

the Kumaraswamy-G cdf. The BKw-G distribution with parameters
$a, b, \alpha, \beta > 0$ has

$$F(x) = I_{K(x)}(\alpha, \beta) = \frac{1}{B(\alpha,\beta)} \int_0^{K(x)} t^{\alpha-1} (1-t)^{\beta-1}\, dt ,$$

$$f(x) = \frac{a\, b\, g(x)\, G(x)^{a-1} \left[1 - G(x)^a\right]^{b\beta - 1} K(x)^{\alpha-1}}{B(\alpha, \beta)} ,$$

and a closed-form quantile function

$$Q(u) = G^{-1}\!\left( \left\{ 1 - \left[1 - z_u\right]^{1/b} \right\}^{1/a} \right), \qquad z_u = \text{qbeta}(u; \alpha, \beta).$$

The hazard is $h(x) = f(x) / [1 - F(x)]$. In the notation of Handique et al.
(2017), `alpha` $= m$ and `beta` $= n$.

**Mixture representation.** Expanding $K^{\alpha-1}$ gives

$$f(x) = \sum_{j=0}^{\infty} w_j \, f_{\text{KwG}}\big(x;\, a,\, b(\beta + j)\big), \qquad w_j = \frac{(-1)^j \binom{\alpha-1}{j}}{B(\alpha,\beta)\,(\beta + j)},$$

so properties of BKw-G follow from those of Kw-G (`bkwg_mixture()`,
`dbkwg_series()`).

---

## Installation

```r
# from GitHub
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("yamenetoo/Statistical-Repository", subdir = "BKwG")

# or from a source tarball
install.packages("BKwG_0.1.0.tar.gz", repos = NULL, type = "source")
```

The package depends only on base R (`stats`, `graphics`). `numDeriv` is used
for the observed information matrix when installed (otherwise `optimHess`),
and `testthat` runs the tests.

---

## Quick start

```r
library(BKwG)

## simulate from a BKw-Weibull distribution
set.seed(123)
x <- rbkwg(300, a = 2, b = 0.8, alpha = 0.7, beta = 1.5,
           baseline = "weibull", par = c(shape = 1.5, scale = 2))

## fit the full model by maximum likelihood
fit <- fit_bkwg(x, baseline = "weibull")
summary(fit)
confint(fit, type = "log")
plot(fit)                        # histogram + density, cdf, P-P and Q-Q plots

## is the extra flexibility needed?
lrt_bkwg(fit, "KwG")             # H0: alpha = beta = 1
lrt_bkwg(fit, "G")               # H0: a plain Weibull is enough

## goodness of fit and model comparison
gof_bkwg(fit, B = 200, seed = 1)
compare_bkwg(x, baseline = c("weibull", "gamma", "lnorm"))
```

---

## Baseline distributions

Every function takes `baseline =` (a name or a `bkwg_baseline()` object) and
`par =` (the baseline parameters, named or in the order below). Missing
parameters take their defaults.

| `baseline` | Parameters (`par`) | Support | Baseline cdf `G(x)` |
|---|---|---|---|
| `"exp"` | `rate` | $(0, \infty)$ | $1 - e^{-\lambda x}$ |
| `"weibull"` | `shape`, `scale` | $(0, \infty)$ | $1 - e^{-(x/s)^k}$ |
| `"norm"` | `mean`, `sd` | $\mathbb{R}$ | $\Phi((x-\mu)/\sigma)$ |
| `"lnorm"` | `meanlog`, `sdlog` | $(0, \infty)$ | $\Phi((\log x-\mu)/\sigma)$ |
| `"gamma"` | `shape`, `rate` | $(0, \infty)$ | regularised gamma |
| `"logis"` | `location`, `scale` | $\mathbb{R}$ | $1/(1+e^{-(x-\mu)/s})$ |
| `"gumbel"` | `location`, `scale` | $\mathbb{R}$ | $\exp\{-e^{-(x-\mu)/s}\}$ |
| `"lomax"` | `shape`, `scale` | $(0, \infty)$ | $1-(1+x/\lambda)^{-\alpha}$ |
| `"llogis"` | `shape`, `scale` | $(0, \infty)$ | $1/(1+(x/s)^{-k})$ |
| `"frechet"` | `shape`, `scale` | $(0, \infty)$ | $\exp\{-(x/s)^{-k}\}$ |
| `"burrxii"` | `c`, `k`, `scale` | $(0, \infty)$ | $1-[1+(x/s)^c]^{-k}$ |
| `"unif"` | — | $(0, 1)$ | $x$ (gives the beta-Kumaraswamy law on the unit interval) |

```r
bkwg_baselines()                     # the same table inside R
dbkwg(1:3, 2, 1, 1.5, 2, baseline = "lomax", par = c(shape = 3, scale = 2))
```

---

## Sub-models

| `submodel` | Name | Restriction | Free shapes |
|---|---|---|---|
| `"BKwG"` | Beta Kumaraswamy-G | none | 4 |
| `"KwG"` | Kumaraswamy-G (Cordeiro & de Castro, 2011) | $\alpha = \beta = 1$ | 2 |
| `"BetaG"` | Beta-G (Eugene et al., 2002) | $a = b = 1$ | 2 |
| `"EKwG"` | Exponentiated Kumaraswamy-G | $\beta = 1$ | 3 |
| `"BEG"` | Beta exponentiated-G | $b = 1$ | 3 |
| `"EG"` | Exponentiated-G (Lehmann type I) | $b = \alpha = \beta = 1$ | 1 |
| `"G"` | The baseline itself | $a = b = \alpha = \beta = 1$ | 0 |

All restrictions lie in the interior of the parameter space, so likelihood-ratio
statistics have the usual $\chi^2$ limit.

---

## Distribution functions

```r
th <- list(a = 2, b = 0.7, alpha = 1.6, beta = 2.3)
x  <- c(0.5, 1, 2)

dbkwg(x, 2, 0.7, 1.6, 2.3, par = c(1.5, 2))                 # density
pbkwg(x, 2, 0.7, 1.6, 2.3, par = c(1.5, 2))                 # cdf
pbkwg(x, 2, 0.7, 1.6, 2.3, par = c(1.5, 2), lower.tail = FALSE, log.p = TRUE)
qbkwg(c(.1, .5, .9), 2, 0.7, 1.6, 2.3, par = c(1.5, 2))     # quantiles
rbkwg(5, 2, 0.7, 1.6, 2.3, par = c(1.5, 2))                 # random numbers
sbkwg(x, 2, 0.7, 1.6, 2.3, par = c(1.5, 2))                 # survival
hbkwg(x, 2, 0.7, 1.6, 2.3, par = c(1.5, 2))                 # hazard
chbkwg(x, 2, 0.7, 1.6, 2.3, par = c(1.5, 2))                # cumulative hazard

## shapes of the density and the hazard
plot_bkwg(2, 1, 1.5, 2, par = c(2, 1), to = 3)
plot_bkwg(0.5, 3, 2, 0.5, par = c(2, 1), from = 0, to = 3, add = TRUE, col = 2)
plot_bkwg(2, 0.4, 0.5, 3, par = c(1, 1), what = "hazard", to = 4)
```

Shape parameters are scalars; `x`, `q` and `p` are vectors. Random numbers are
generated by inverting a beta variate, so they are exact.

---

## Structural properties

All quantities are computed by numerical integration over the quantile
function, so they work for every baseline.

```r
mbkwg(1:4, 2, 1, 1.5, 2, par = c(2, 1))          # raw moments E[X^r]
mgf_bkwg(c(-1, 0.5), 2, 1, 1.5, 2, par = c(2, 1)) # moment generating function
bkwg_stats(2, 1, 1.5, 2, par = c(2, 1))          # mean, variance, sd, skewness,
                                                 # kurtosis, median, Bowley, Moors
entropy_bkwg(2, 1, 1.5, 2, par = c(2, 1))                         # Shannon
entropy_bkwg(2, 1, 1.5, 2, par = c(2, 1), type = "renyi", delta = 0.5)
mrl_bkwg(c(0.5, 1), 2, 1, 1.5, 2, par = c(2, 1)) # mean residual life

dos_bkwg(1, i = 10, N = 10, 2, 1, 1.5, 2, par = c(2, 1))  # density of the maximum
pos_bkwg(1, i = 1,  N = 10, 2, 1, 1.5, 2, par = c(2, 1))  # cdf of the minimum

bkwg_mixture(alpha = 2.5, beta = 1.5, J = 20)              # mixture weights
score_bkwg(x = rbkwg(100, 2, 1, 1.5, 2, par = c(2, 1)),    # analytical score
           2, 1, 1.5, 2, par = c(2, 1))
```

| Quantity | Formula used |
|---|---|
| Raw moments | $E[X^r] = \int_0^1 Q(u)^r\,du$ |
| MGF | $M(t) = \int_0^1 e^{tQ(u)}\,du$ |
| Shannon entropy | $-\int_0^1 \log f(Q(u))\,du$ |
| Rényi entropy | $(1-\delta)^{-1}\log\int_0^1 f(Q(u))^{\delta-1}\,du$ |
| Mean residual life | $\int_{F(t)}^1 Q(u)\,du / [1-F(t)] - t$ |
| Bowley skewness | $[Q(3/4)+Q(1/4)-2Q(1/2)]/[Q(3/4)-Q(1/4)]$ |
| Moors kurtosis | $[Q(7/8)-Q(5/8)+Q(3/8)-Q(1/8)]/[Q(6/8)-Q(2/8)]$ |
| Order statistics | $f_{i:N}=\frac{N!}{(i-1)!(N-i)!} f F^{i-1}(1-F)^{N-i}$, $F_{i:N}=I_{F}(i,N-i+1)$ |

When a moment does not exist (heavy tails), the function returns `NA` with a
warning instead of a wrong number.

---

## Parameter estimation

```r
fit_bkwg(x, baseline = "weibull", submodel = "BKwG", fixed = NULL,
         start = NULL, method = "mle", nstart = 3, hessian = TRUE)
```

| `method` | Estimator | Objective minimised |
|---|---|---|
| `"mle"` | Maximum likelihood | $-\sum \log f(x_i)$ |
| `"mps"` | Maximum product of spacings | $-\sum \log[F(x_{(i)}) - F(x_{(i-1)})]$ |
| `"ad"` | Minimum Anderson–Darling distance | $-n - \frac1n\sum (2i-1)[\log F(x_{(i)}) + \log S(x_{(n+1-i)})]$ |
| `"cvm"` | Minimum Cramér–von Mises distance | $\frac{1}{12n} + \sum [F(x_{(i)}) - \frac{2i-1}{2n}]^2$ |
| `"lse"` | Ordinary least squares | $\sum [F(x_{(i)}) - \frac{i}{n+1}]^2$ |
| `"wlse"` | Weighted least squares | $\sum \frac{(n+1)^2(n+2)}{i(n-i+1)} [F(x_{(i)}) - \frac{i}{n+1}]^2$ |

- **Sub-models and fixed parameters**: `submodel = "KwG"`, or
  `fixed = c(scale = 1, beta = 1)` for any combination of shape and baseline
  parameters.
- **Optimisation**: on an unconstrained scale (log for positive parameters,
  logit for bounded ones), Nelder–Mead followed by BFGS, from the best
  `nstart` points of a grid of shape values combined with moment/quantile
  starting values for the baseline. Supply `start =` to override.
- **Standard errors** come from the inverse observed information matrix
  (methods `"mle"` and `"mps"`).

Methods for fitted objects: `print`, `summary` (Wald z-tests and information
criteria), `coef` (`complete = TRUE` adds the fixed parameters), `vcov`,
`logLik`, `AIC`, `BIC`, `nobs`, `confint` (`type = "wald"` or `"log"`),
`residuals` (`"quantile"`, `"pit"`, `"cox-snell"`), `quantile` and `plot`.

```r
fit_mps <- fit_bkwg(x, baseline = "weibull", method = "mps")
fit_kw  <- fit_bkwg(x, baseline = "weibull", submodel = "KwG")
fit_fix <- fit_bkwg(x, baseline = "weibull", fixed = c(scale = 2))
coef(fit_kw, complete = TRUE)
qqnorm(residuals(fit)); abline(0, 1)
```

> **Identifiability.** When $\alpha = 1$, the parameters $b$ and $\beta$ enter
> the model only through the product $b\beta$, and scale-type baseline
> parameters can be partly confounded with the shape parameters. If the
> standard errors are very large, fix a scale parameter or use a sub-model
> selected by `lrt_bkwg()` / `compare_bkwg()`.

---

## Hypothesis tests and goodness of fit

```r
## likelihood ratio test of a nested sub-model (refits it automatically)
lrt_bkwg(fit, "KwG")
lrt_bkwg(fit, fit_kw)               # or pass an already fitted sub-model

## joint Wald test
wald_bkwg(fit, c(alpha = 1, beta = 1))

## goodness of fit
gof_bkwg(fit)                       # W*, A*, KS (asymptotic KS p-value)
gof_bkwg(fit, B = 500, seed = 1)    # parametric bootstrap p-values
```

`gof_bkwg` reports the Cramér–von Mises $W^*$ and Anderson–Darling $A^*$
statistics of Chen and Balakrishnan (1995) and the Kolmogorov–Smirnov
statistic. With `B > 0`, samples are drawn from the fitted model and the model
is re-estimated on each one, so the p-values account for parameter estimation.
`lrt_bkwg` and `wald_bkwg` return standard `"htest"` objects.

---

## Model comparison

```r
tab <- compare_bkwg(x, baseline = c("weibull", "gamma", "lnorm"),
                    submodels = c("BKwG", "KwG", "BetaG", "EKwG", "G"))
tab                                  # sorted by AIC
attr(tab, "fits")[["KwG-gamma"]]     # any fitted model in the table
```

Columns: number of free parameters `k`, `logLik`, `AIC`, `CAIC`, `BIC`, `HQIC`,
`W*`, `A*`, `KS` and its p-value — the table usually reported in papers that
introduce a new distribution.

---

## User-defined baselines

Any continuous distribution can be used. Only the cdf and pdf are required; a
quantile function is optional (otherwise it is found numerically) and
`log.p`/`log` arguments improve tail accuracy.

```r
rayleigh <- bkwg_baseline("rayleigh",
  pG = function(q, par) ifelse(q > 0, 1 - exp(-q^2 / (2 * par[["sigma"]]^2)), 0),
  dG = function(x, par) ifelse(x > 0, x / par[["sigma"]]^2 *
                               exp(-x^2 / (2 * par[["sigma"]]^2)), 0),
  qG = function(p, par) par[["sigma"]] * sqrt(-2 * log1p(-p)),
  par = c(sigma = 1), lower = 0, upper = Inf, support = c(0, Inf),
  start = function(x) c(sigma = sqrt(mean(x^2) / 2)))

y <- rbkwg(300, 2, 1, 1.5, 0.8, baseline = rayleigh, par = c(sigma = 2))
fit_bkwg(y, baseline = rayleigh)
```

---

## Numerical details

- Everything is computed on the log scale with `log1mexp` / `log1pexp`
  helpers: $\log G$, $\log(1-G^a)$, $\log K$ and $\log(1-K)$. The cdf switches
  between `pbeta(K, α, β)` and `pbeta(1-K, β, α, lower.tail = FALSE)` so that
  both tails keep full precision.
- The quantile function uses `qbeta` with the complementary tail and
  `log1mexp`, and random generation inverts `rbeta(n, β, α)`.
- The density, cdf, quantile function, analytical score and mixture
  representation were checked numerically against direct integration and
  finite differences for seven baselines.

---

## Function reference

| Function | Purpose |
|---|---|
| `dbkwg`, `pbkwg`, `qbkwg`, `rbkwg` | density, cdf, quantile function, random generation |
| `sbkwg`, `hbkwg`, `chbkwg` | survival, hazard, cumulative hazard |
| `score_bkwg` | analytical score of $(a, b, \alpha, \beta)$ |
| `bkwg_baselines`, `bkwg_baseline` | list built-in baselines; define a new one |
| `mbkwg`, `mgf_bkwg`, `bkwg_stats` | moments, mgf, summary measures |
| `entropy_bkwg`, `mrl_bkwg` | Shannon / Rényi entropy, mean residual life |
| `dos_bkwg`, `pos_bkwg` | order statistics |
| `bkwg_mixture`, `dbkwg_series` | mixture weights and series density |
| `plot_bkwg` | pdf / cdf / survival / hazard plots |
| `fit_bkwg`, `bkwg_submodels` | estimation; table of sub-models |
| `lrt_bkwg`, `wald_bkwg` | likelihood ratio and Wald tests |
| `gof_bkwg` | W\*, A\*, KS with bootstrap p-values |
| `compare_bkwg` | comparison table over sub-models and baselines |

Every function has a help page with examples: `help(package = "BKwG")`.

---

## Examples notebook

[`notebooks/BKwG_examples.ipynb`](notebooks/BKwG_examples.ipynb) walks through
every feature: shapes of the density and hazard, properties, simulation,
estimation with all six methods, tests, model comparison on a real data set, a
user-defined baseline and a small Monte Carlo study. Open it in Colab with the
badge at the top (choose *Runtime → Change runtime type → R*).

---

## Testing

```r
devtools::test("BKwG")           # testthat suite
devtools::check("BKwG")          # full R CMD check
```

The tests cover: the density integrating to one for every baseline, cdf versus
integrated density, quantile/cdf inversion, closed forms of all sub-models,
random generation (KS test), the analytical score against a numerical
gradient, the mixture representation, moments, entropies and mean residual
life against known exponential/normal/Weibull results, order statistics,
nesting of the fitted log-likelihoods, all estimation methods, and every
method for fitted objects.

---

## Citation

```r
citation("BKwG")
```

```bibtex
@Manual{BKwG,
  title  = {{BKwG}: The Beta Kumaraswamy-G Family of Distributions},
  author = {Mohamad Yamen {AL Mohamad} and Ali Akbar Heydari and Hossein {Jabbari Khamnei}},
  year   = {2026},
  note   = {R package version 0.1.0},
  url    = {https://github.com/yamenetoo/Statistical-Repository/tree/main/BKwG}
}
```

Please also cite the paper that introduced the family:

> Handique, L., Chakraborty, S. and Ali, M. M. (2017). Beta generated
> Kumaraswamy-G family of distributions. *Pakistan Journal of Statistics*,
> 33(6), 467–490. Preprint: [arXiv:1603.00634](https://arxiv.org/abs/1603.00634).

Other references:

- Cordeiro, G. M. and de Castro, M. (2011). A new family of generalized
  distributions. *Journal of Statistical Computation and Simulation*, 81(7),
  883–898.
- Eugene, N., Lee, C. and Famoye, F. (2002). Beta-normal distribution and its
  applications. *Communications in Statistics – Theory and Methods*, 31(4),
  497–512.
- Chen, G. and Balakrishnan, N. (1995). A general purpose approximate
  goodness-of-fit test. *Journal of Quality Technology*, 27, 154–161.

---

## Authors and license

| Author | ORCID | Email |
|---|---|---|
| **Mohamad Yamen AL Mohamad** (maintainer) | [0000-0003-4658-9095](https://orcid.org/0000-0003-4658-9095) | yamenmohamad@tabrizu.ac.ir |
| **Ali Akbar Heydari** | [0000-0002-1236-3189](https://orcid.org/0000-0002-1236-3189) | heydari@tabrizu.ac.ir |
| **Hossein Jabbari Khamnei** | [0000-0003-0902-0967](https://orcid.org/0000-0003-0902-0967) | h_jabbari@tabrizu.ac.ir |

Department of Statistics, Faculty of Mathematics, Statistics and Computer
Science, University of Tabriz, Tabriz, Iran.

License: GPL (≥ 3).
