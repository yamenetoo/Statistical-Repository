## Submission

This is a new submission (resubmission after the first review).

Changes in response to the reviewer:

* References describing the methods were added to the Description field of
  DESCRIPTION in the form authors (year) <doi:...>; the arXiv e-print is cited
  through its arXiv DOI.
* No example, test or vignette uses more than one core. The only script that
  uses parallel code (`inst/validation/exponential-validation.R`, not run during
  checks) is limited to at most 2 cores.

## Test environments

* local Ubuntu 24.04, R 4.3.3
* win-builder, R-release (4.6.2) and R-devel (2026-11-03 r90835)

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Mohamad Yamen Al-Mohamad <yamenmohamad@tabrizu.ac.ir>'
  New submission.

  Possibly misspelled words in DESCRIPTION: ARL, Kumaraswamy, Lindley.
  These are correct: ARL is the standard abbreviation of "average run
  length", Kumaraswamy and Lindley are the names of the two distributions
  mentioned, and the other flagged words (Efron, Gandy, Heydari, Jabbari,
  Khamnei, Kvaloy, Lami, Padgett, Spurrier, Tibshirani) are author names in
  the references.

## Downstream dependencies

There are no downstream dependencies.
