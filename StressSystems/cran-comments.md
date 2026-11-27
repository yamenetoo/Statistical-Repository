## Submission

This is the first submission of StressSystems.

StressSystems computes and estimates stress-strength reliability
R = P(Y < X) for arbitrary pairs of continuous distributions, for single
components and for multi-stress, three-variable, s-out-of-k and cascade
systems, with several estimation methods, interval estimation for R and
Monte Carlo tools.

## Test environments

* local: Ubuntu 24.04, R 4.3.3
* win-builder: R-release (R 4.6.2, Windows Server 2022)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.
* Possibly misspelled words in DESCRIPTION (Amin, Bhattacharyya, Cheng,
  Cramer, Kotz, Lumelskii, Mises, Pensky, spacings, von) are author
  names and standard statistical terms.

## Notes for reviewers

* Examples that run bootstrap, MCMC or larger Monte Carlo studies are wrapped
  in \donttest{} because they take more than 5 seconds.
* `ss_simulate()` sets the random seed only when the user passes `seed`
  (default `NULL`), so results can be reproduced on request.
* Parallel execution uses at most the number of cores given by the user
  (default 1).
