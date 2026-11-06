# StressSystems 0.1.0

* First release.
* Distributions: built-in Topp-Leone and generalized Rayleigh families, any
  R family found by name, and user-supplied density/distribution functions.
* Reliability: `ss_R()`, `ss_multistress()`, `ss_R3()`, `ss_sk()` and
  `ss_cascade()`.
* Estimation: maximum likelihood, maximum product of spacings, least squares,
  weighted least squares, Cramer-von Mises, Anderson-Darling, jackknife and
  bootstrap bias correction, and Bayesian Metropolis-Hastings.
* Inference for R: delta-method, parametric bootstrap and posterior intervals
  (`ss_estimate_R()`).
* Monte Carlo studies: `ss_simulate()` and `ss_simulate_grid()`.
