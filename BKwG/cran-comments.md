## Submission

This is the first submission of BKwG.

BKwG implements the beta Kumaraswamy-G family of distributions
(Handique, Chakraborty and Ali, 2017, Pakistan Journal of Statistics 33(6),
467-490; <doi:10.48550/arXiv.1603.00634>) for arbitrary baseline
distributions, with estimation, hypothesis tests and goodness-of-fit tools.

## Test environments

* local: R release
* win-builder: R-devel and R-release
* R-hub / GitHub Actions: Linux, macOS, Windows

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Notes for reviewers

* Examples that fit the full six-parameter model or run a parametric
  bootstrap are wrapped in \donttest{} because they take more than 5 seconds.
* Slower model-fitting tests are skipped on CRAN with testthat::skip_on_cran().
