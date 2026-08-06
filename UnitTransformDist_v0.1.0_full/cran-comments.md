## Submission notes

This is the first submission of unitTransformDist to CRAN.

## Test environments

* Local development environment (R 4.3.3, Ubuntu), no network access to
  CRAN mirrors during development -- `R CMD check --as-cran` has NOT yet
  been run end-to-end in a network-connected environment and must be run
  before an actual submission. All mathematical contracts and every
  `@examples` block have been verified by direct `Rscript` execution
  instead (see `tests/testthat/` and the development log in NEWS.md).

## R CMD check results

Not yet run in a connected environment. TODO before submission:
* Run `devtools::check(cran = TRUE)` locally.
* Run `rhub::check_for_cran()` across platforms.
* Regenerate `man/*.Rd` and `NAMESPACE` with `roxygen2::roxygenise()`
  (the current NAMESPACE was written by hand to mirror what roxygen2
  would produce, because this development environment could not reach
  CRAN to install roxygen2).

## Downstream dependencies

None (first submission).
