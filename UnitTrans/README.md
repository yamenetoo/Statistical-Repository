# UnitTrans

**Distributions on the unit interval via transformations of baseline laws**

`UnitTrans` builds new probability distributions on (0, 1) by applying a monotone
transformation `V = g(Y)` — rational, exponential, logit, arcsecant-hyperbolic,
probit or complementary log-log — to a baseline random variable `Y` with any
available `d/p/q/r` functions.

It provides:

- density, distribution, quantile and random generation (`dtrans`, `ptrans`, `qtrans`, `rtrans`)
- maximum-likelihood fitting with information criteria (`fit_trans`)
- automatic model search over baseline/transformation pairs (`fit_all`, `best_model`)
- goodness-of-fit tests with parametric bootstrap (`gof_trans`)
- moments, survival and hazard functions (`moments_trans`, `survival_trans`, `hazard_trans`)
- simulation studies for validating the estimator (`simulation_study`)

## Installation

```r
# install.packages("remotes")
remotes::install_github("yamenetoo/Statistical-Repository", subdir = "UnitTrans", build_vignettes = TRUE)
```

## Quick start

```r
library(UnitTrans)

trans_list()                                        # available transformations
detect_transformations("gamma", c(shape = 3, rate = 1))

set.seed(1)
x <- rtrans(300, "gamma", c(2, 1), "rational")      # simulate data

fit <- fit_trans(x, "gamma", "rational")            # fit one model
summary(fit)

res  <- fit_all(x, c("gamma", "lnorm", "weibull"), criterion = "AIC")
best <- best_model(res)                             # rank all candidates
plot(best)
gof_trans(best, B = 50, seed = 1)
confint(best)
```

See the vignette for the full workflow: `vignette("workflow", package = "UnitTrans")`.

## Running the tests

```sh
R CMD build UnitTrans
R CMD check UnitTrans_0.1.0.tar.gz
```

or inside R: `testthat::test_local()`.

## Author

Mohamad Yamen Al-Mohamad — Department of Statistics, University of Tabriz

## License

GPL-3
