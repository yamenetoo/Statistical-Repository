# unitTransformDist

Generalized (0,1)-interval distributions from *any* base distribution on
(0, Inf), via the odds-inverting transform `Y = 1/(1+X)`.

Unlike registry-restricted designs, the base distribution is supplied
directly as density/CDF/quantile/random-generation functions -- so this
works with built-in R distributions (`dgamma`/`pgamma`/...), or with a
fully custom research distribution that has never been registered
anywhere, even one with no closed-form quantile function.

## Design principles

See the header comment in `R/core.R` for the full list this package was
built against: mathematical contract before code, input validation at the
exported boundary, log-scale numerical stability, no hidden global
`set.seed()`, honest convergence reporting, and base-R S3 idiom
(`print`, `summary`, `coef`, `logLik`, `vcov`) throughout.

## Quick start

```r
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
```

## Status

Version 0.1.0. `R CMD check` has not yet been run end-to-end in a
network-connected environment (this development environment could not
reach CRAN); all mathematical contracts and examples above have been
verified by direct execution. Run `R CMD check --as-cran` before any
CRAN submission.

## License

MIT
