## Worked examples of Ross (2000), "Membership Functions, Fuzzification and
## Defuzzification", reproduced with FuzzyMFD.
##   source(system.file("examples", "ross2000_examples.R", package = "FuzzyMFD"))

library(FuzzyMFD)

## Section 1, Figure 1.1: crisp set A = [5, 7] and fuzzy set H "near 6 feet"
hs <- ross2000_sets("height")
membership(hs$A, c(6, 4.99))          # chi_A(x1) = 1, chi_A(x2) = 0
membership(hs$H, c(5.917, 6, 7))
plot_fuzzy_sets(hs, xlab = "height (ft)", main = "Figure 1.1")

## Section 2: features of membership functions (Figures 1.2-1.5)
Tz <- fs_trapezoidal(1, 3, 5, 8, name = "A")
summary(Tz)
convex <- fuzzy_set(c(0, 2, 5), c(0, 1, 0))
nonconvex <- fuzzy_set(c(0, 2, 3, 4, 6), c(0, 1, 0.7, 0.9, 0))
c(is_convex(convex), is_convex(nonconvex))
AB <- fuzzy_intersection(fs_triangular(0, 2, 5), fs_trapezoidal(2, 4, 6, 8))
is_convex(AB)                          # Figure 1.5

## Section 3: fuzzification of a voltage reading with +/- 1 percent error
R <- fuzzify(10, spread = 0.01, relative = TRUE, name = "reading")
fuzzy_support(R)
v <- ross2000_sets("voltage")
agreement(v$low, 5.8)                  # crisp reading (Figure 1.7a)
agreement(v$medium, fuzzify(7.2, 1))   # fuzzy reading (Figure 1.7b)

## Section 4.2: rank ordering (Table 1.1, Figure 1.9)
data(ross_colors)
ro <- rank_ordering(ross_colors)
ro
plot(ro, main = "Figure 1.9: best colour")

## Section 4.3: inductive reasoning (Tables 1.2-1.5, Figures 1.13-1.14)
data(ross_ellipse)
entropy_threshold(ross_ellipse$x, ross_ellipse$class,
                  candidates = c(0.7, 1.0, 1.5, 3.45))           # Table 1.3
neg <- ross_ellipse$x < 1.5
entropy_threshold(ross_ellipse$x[neg], ross_ellipse$class[neg],
                  candidates = c(0.175, 0.35, 0.7))               # Table 1.4
entropy_threshold(ross_ellipse$x[!neg], ross_ellipse$class[!neg],
                  candidates = 27.5)                              # Table 1.5
res <- inductive_mf(ross_ellipse$x, ross_ellipse$class, levels = 2,
                    range = c(0, 100))
res$partition
plot_fuzzy_sets(partition_sets(1.5, range = c(0, 100)), main = "Figure 1.13")
plot_fuzzy_sets(res$sets, main = "Figure 1.14", xlab = "x = a/b")

## Section 5: defuzzification
out <- ross2000_sets("output")
plot_fuzzy_sets(list(C = fuzzy_union(out)), main = "Figure 1.15c")

B <- ross2000_sets("map")              # Figures 1.20-1.22
c(centroid = defuzzify(B, "centroid"), # 4.95 (reported as 4.9 in the chapter)
  weighted = defuzzify(B, "weighted"), # 5.41
  meanmax = defuzzify(B, "meanmax"))   # 6.5
plot_fuzzy_sets(list(union = fuzzy_union(B)), legend = FALSE,
                xlab = "MAP (mm Hg)", main = "Figures 1.23-1.25")
abline(v = c(defuzzify(B, "centroid"), defuzzify(B, "weighted"),
             defuzzify(B, "meanmax")), lty = 2:4)
