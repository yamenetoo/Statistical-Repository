# StressSystems

Stress-strength reliability for **any pair of continuous distributions**, for
single components and for multicomponent systems, with a full set of
estimators, confidence intervals and a ready Monte Carlo engine.

## Install

```r
# install.packages("remotes")
remotes::install_github("yamenetoo/Statistical-Repository", subdir = "StressSystems")
```

## Distributions

Any family that follows R's `d/p/q/r` naming works by name; the parameters are
read from the density's arguments:

```r
library(StressSystems)
X <- ss_dist("toppleone")                       # built in
Y <- ss_dist("genrayleigh")                     # built in (shape, rate)
W <- ss_dist("weibull", par = c(shape = 2, scale = 1))
G <- ss_dist("gamma")                           # shape, rate
```

User-defined families only need a density and a distribution function:

```r
D <- ss_dist(d = mydens, p = mycdf, lower = 0, upper = Inf,
             constraints = list(a = "positive", b = "unit"))
```

## Reliability

| Function | Quantity |
|---|---|
| `ss_R()` | R = P(cY < X) |
| `ss_multistress()` | P(cY_j < X for all j), shared or independent components |
| `ss_R3()` | P(Z < Y < X) |
| `ss_sk()` | s-out-of-k system |
| `ss_cascade()` | n-stage cascade with stress factor c^(i-1) |

R is computed as the integral of F_Y(Q_X(u)/c) over u in (0, 1), which stays
bounded even when a density is unbounded at the edge of its support.

## Estimation

`ss_fit(x, dist, method)` with `method` one of `"mle"`, `"mps"`, `"lse"`,
`"wlse"`, `"cvm"`, `"ad"`, `"jackknife"`, `"bootstrap"`, `"bayes"`.
Optimisation runs on an unconstrained scale from screened random starts, so
no starting values are needed.

```r
x <- rtoppleone(30, 0.7); y <- rgenrayleigh(30, 1, 2)
ss_estimate_R(x, y, X, Y, method = "mle")                 # delta-method CI
ss_estimate_R(x, y, X, Y, method = "mle", ci = "boot")    # bootstrap CI
ss_estimate_R(x, y, X, Y, method = "bayes")               # posterior interval
```

## Simulation

```r
sim <- ss_simulate(X, Y, par_x = c(shape = 0.7), par_y = c(shape = 1, rate = 2),
                   n = c(10, 25, 50), MC = 1000,
                   methods = c("mle", "mps", "jackknife", "bayes"))
sim$summary      # bias, MSE, coverage, interval length, failures
sim$parameters   # bias and MSE of every parameter
```

`ss_simulate_grid()` runs the same study over grids of parameter values. See
`inst/scripts/toppleone_genrayleigh_study.R` for a complete study.

## Authors

Mohamad Yamen AL Mohamad, Ali Akbar Heydari, Hossein Jabbari Khamnei —
Department of Statistics, University of Tabriz.
