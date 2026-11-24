# FuzzyMFD

**Membership functions, fuzzification and defuzzification in R.**

FuzzyMFD implements the methods of

> Ross, T. J. (2000). Membership functions, fuzzification and defuzzification. In P. S. Szczepaniak, P. J. G. Lisboa and J. Kacprzyk (eds.), *Fuzzy Systems in Medicine*, 49–78. Springer-Verlag, Berlin Heidelberg.

Fuzzy sets on the real line are stored as **exact piecewise-linear membership functions** (knots and membership values, zero outside the knots). Triangles, trapezoids, crisp intervals and all their unions and intersections are therefore represented without grid error, and the centroid is computed in closed form.

## Installation

```r
# install.packages("remotes")
remotes::install_github("yamenetoo/Statistical-Repository", subdir = "FuzzyMFD")
```

Worked examples: [`notebooks/FuzzyMFD_examples.ipynb`](notebooks/FuzzyMFD_examples.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/yamenetoo/Statistical-Repository/blob/main/FuzzyMFD/notebooks/FuzzyMFD_examples.ipynb), and the script `inst/examples/ross2000_examples.R`, which reproduces every numerical example of the chapter.

## What is in the package

| Section of the chapter | Functions |
|---|---|
| 1. Set membership | `fuzzy_set`, `fuzzy_set_fun`, `fs_triangular`, `fs_trapezoidal`, `fs_rectangular`, `fs_singleton`, `fs_gaussian`, `membership`, `mf_indicator` (Eq. 1), `mf_triangular`, `mf_trapezoidal`, `mf_gaussian`, `mf_bell`, `mf_sigmoid` |
| 2. Features of the membership function | `fuzzy_core`, `fuzzy_support`, `fuzzy_boundary`, `fuzzy_height`, `is_normal`, `is_subnormal`, `is_convex` (Eq. 3), `is_fuzzy_number`, `crossover_points`, `fuzzy_prototype`, `summary()` |
| Set operations | `fuzzy_union` (max), `fuzzy_intersection` (min), `fuzzy_complement` — exact for piecewise-linear sets |
| 3. Fuzzification | `fuzzify` (crisp readings with absolute or relative imprecision), `agreement` (Figure 1.7) |
| 4.2 Rank ordering | `rank_ordering` with `print`/`plot` methods (Table 1.1, Figure 1.9) |
| 4.3 Inductive reasoning | `entropy_threshold` (Eqs. 4–13), `entropy_partition` (PRI, SEC1–2, TER1–4, …), `partition_sets`, `inductive_mf` (Figures 1.11–1.14) |
| 5. Defuzzification | `defuzzify`, `defuzz_max` (Eq. 15), `defuzz_centroid` (Eq. 16), `defuzz_weighted` (Eq. 17), `defuzz_meanmax` (Eq. 18) |
| Examples and plots | `ross_colors`, `ross_ellipse`, `ross2000_sets`, `plot()`, `lines()`, `plot_fuzzy_sets` |

## Quick start

```r
library(FuzzyMFD)

## features of a membership function
A <- fs_trapezoidal(1, 3, 5, 8, name = "A")
summary(A)

## rank ordering (Table 1.1)
data(ross_colors)
rank_ordering(ross_colors)

## inductive reasoning: thresholds 0.35 (SEC1), 1.5 (PRI), 27.5 (SEC2)
data(ross_ellipse)
res <- inductive_mf(ross_ellipse$x, ross_ellipse$class, levels = 2, range = c(0, 100))
res$partition
plot_fuzzy_sets(res$sets, xlab = "x = a/b")

## defuzzification of the mean-arterial-pressure example (Section 5.2)
B <- ross2000_sets("map")
defuzzify(B, "centroid")  # 4.953  (the chapter reports 4.9)
defuzzify(B, "weighted")  # 5.417  (5.41)
defuzzify(B, "meanmax")   # 6.5
```

## Notes on the inductive-reasoning example

* Thresholds are chosen among candidate points, by default the mid-points between adjacent distinct values; the chapter evaluates only some of them. Pass `candidates` to reproduce Tables 1.3–1.5 exactly; with all mid-points the primary and negative-side thresholds are the same (1.5 and 0.35).
* The positive side of the example contains a single class, so the entropy carries no class information there. With `pure = "middle"` (the default) such a region is split at its middle candidate, which gives the chapter's 27.5; `pure = "entropy"` uses the minimum of `S` and `pure = "none"` does not split.

## Author

Mohamad Yamen AL Mohamad, Department of Statistics, University of Tabriz — yamenmohamad@tabrizu.ac.ir — ORCID [0009-0009-1665-7203](https://orcid.org/0009-0009-1665-7203)

License: GPL (>= 3)
