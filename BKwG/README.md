# BKwG: The Beta Kumaraswamy-G Family of Distributions

`BKwG` is an R package for the beta Kumaraswamy-G (BKw-G) family of
distributions (Handique, Chakraborty and Ali, 2017). For any baseline
distribution with cdf G and pdf g, the BKw-G distribution has

```
K(x) = 1 - [1 - G(x)^a]^b                       (Kumaraswamy-G cdf)
F(x) = I_{K(x)}(alpha, beta)                     (regularised incomplete beta)
f(x) = a b g(x) G(x)^(a-1) [1 - G(x)^a]^(b*beta - 1) K(x)^(alpha - 1) / B(alpha, beta)
Q(u) = G^{-1}( {1 - [1 - z_u]^(1/b)}^(1/a) ),   z_u = qbeta(u, alpha, beta)
```

with `a, b, alpha, beta > 0` (`alpha = m`, `beta = n` in Handique et al.).

## Installation

```r
# from the source tarball
install.packages("BKwG_0.1.0.tar.gz", repos = NULL, type = "source")

# or from GitHub
# remotes::install_github("yamenetoo/Statistical-Repository", subdir = "BKwG")
```

## What is in the package

| Area | Functions |
|---|---|
| Distribution | `dbkwg`, `pbkwg`, `qbkwg`, `rbkwg`, `sbkwg`, `hbkwg`, `chbkwg` |
| Baselines | 12 built-in (`exp`, `weibull`, `norm`, `lnorm`, `gamma`, `logis`, `gumbel`, `lomax`, `llogis`, `frechet`, `burrxii`, `unif`) + any user baseline via `bkwg_baseline()` |
| Score | `score_bkwg` (analytical score of a, b, alpha, beta) |
| Properties | `mbkwg` (moments), `mgf_bkwg`, `bkwg_stats` (mean, variance, skewness, kurtosis, Bowley, Moors), `entropy_bkwg` (Shannon, Renyi), `mrl_bkwg`, `dos_bkwg` / `pos_bkwg` (order statistics), `bkwg_mixture` / `dbkwg_series` (mixture of Kw-G densities), `plot_bkwg` |
| Sub-models | BKw-G, Kw-G, Beta-G, exponentiated Kw-G, beta exponentiated-G, exponentiated-G, G (`bkwg_submodels()`) |
| Estimation | `fit_bkwg` with MLE, maximum product of spacings, Anderson-Darling, Cramer-von Mises, least squares and weighted least squares; fixed parameters; standard errors |
| Inference | `lrt_bkwg` (likelihood ratio), `wald_bkwg`, `gof_bkwg` (W*, A*, KS, parametric bootstrap p-values), `compare_bkwg` (AIC, CAIC, BIC, HQIC table) |
| Methods | `print`, `summary`, `coef`, `vcov`, `logLik`, `AIC`, `BIC`, `nobs`, `confint`, `residuals`, `quantile`, `plot` |

## Example

```r
library(BKwG)

set.seed(123)
x <- rbkwg(300, a = 2, b = 0.8, alpha = 0.7, beta = 1.5,
           baseline = "weibull", par = c(shape = 1.5, scale = 2))

fit <- fit_bkwg(x, baseline = "weibull")
summary(fit)
confint(fit, type = "log")
plot(fit)

# nested sub-models
lrt_bkwg(fit, "KwG")          # H0: alpha = beta = 1
lrt_bkwg(fit, "G")            # H0: plain Weibull
wald_bkwg(fit, c(alpha = 1, beta = 1))

# goodness of fit and model comparison
gof_bkwg(fit, B = 200, seed = 1)
compare_bkwg(x, baseline = c("weibull", "gamma"))

# properties
bkwg_stats(a = 2, b = 0.8, alpha = 0.7, beta = 1.5, par = c(1.5, 2))
entropy_bkwg(a = 2, b = 0.8, alpha = 0.7, beta = 1.5, par = c(1.5, 2))
plot_bkwg(a = 2, b = 0.8, alpha = 0.7, beta = 1.5, par = c(1.5, 2),
          what = "hazard")
```

A user-defined baseline:

```r
rayleigh <- bkwg_baseline("rayleigh",
  pG = function(q, par) ifelse(q > 0, 1 - exp(-q^2 / (2 * par[["sigma"]]^2)), 0),
  dG = function(x, par) ifelse(x > 0, x / par[["sigma"]]^2 *
                               exp(-x^2 / (2 * par[["sigma"]]^2)), 0),
  par = c(sigma = 1), lower = 0, upper = Inf, support = c(0, Inf))

fit_bkwg(rbkwg(200, 2, 1, 1.5, 0.8, baseline = rayleigh, par = c(sigma = 2)),
         baseline = rayleigh)
```

## Reference

Handique, L., Chakraborty, S. and Ali, M. M. (2017). Beta generated
Kumaraswamy-G family of distributions. *Pakistan Journal of Statistics*,
33(6), 467-490. Preprint: arXiv:1603.00634.
