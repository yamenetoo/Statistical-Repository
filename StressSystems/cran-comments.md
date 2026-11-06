## Submission

This is the first submission of StressSystems.

StressSystems computes and estimates stress-strength reliability
R = P(Y < X) for arbitrary pairs of continuous distributions, for single
components and for multi-stress, three-variable, s-out-of-k and cascade
systems, with several estimation methods, interval estimation for R and
Monte Carlo tools.

## Test environments

<!-- Update this section and the results below with the actual output
     of the checks before submitting. -->

* local: R release
* win-builder: R-devel and R-release
* R-hub: Linux, macOS, Windows

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Notes for reviewers

* Examples that run bootstrap, MCMC or larger Monte Carlo studies are wrapped
  in \donttest{} because they take more than 5 seconds.
* `ss_simulate()` sets the random seed only when the user passes `seed`
  (default `NULL`), so results can be reproduced on request.
* Parallel execution uses at most the number of cores given by the user
  (default 1).
