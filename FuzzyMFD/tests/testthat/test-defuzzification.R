test_that("defuzzification reproduces Section 5.2", {
  B <- ross2000_sets("map")
  ## exact centroid of the aggregated output (the chapter rounds to 4.9)
  U <- fuzzy_union(B)
  num <- 0.3 * (1/3) + 0.3 * (3.6^2 - 1) / 2 +
    ((4^3 - 3.6^3) / 3 - 3 * (4^2 - 3.6^2) / 2) / 2 +
    0.5 * (5.5^2 - 4^2) / 2 +
    ((6^3 - 5.5^3) / 3 - 5 * (6^2 - 5.5^2) / 2) + (7^2 - 6^2) / 2 +
    (8 * (8^2 - 7^2) / 2 - (8^3 - 7^3) / 3)
  den <- 0.3 / 2 + 0.3 * 2.6 + ((4^2 - 3.6^2) / 2 - 3 * 0.4) / 2 + 0.5 * 1.5 +
    ((6^2 - 5.5^2) / 2 - 5 * 0.5) + 1 + (8 - (8^2 - 7^2) / 2)
  expect_equal(defuzzify(B, "centroid"), num / den)
  expect_equal(round(defuzzify(B, "centroid"), 3), 4.953)
  expect_equal(defuzzify(B, "weighted"), (0.3 * 2.5 + 0.5 * 5 + 1 * 6.5) / 1.8)
  expect_equal(round(defuzzify(B, "weighted"), 2), 5.42)
  expect_equal(defuzzify(B, "meanmax"), 6.5)
  expect_equal(defuzzify(B, "max"), 6)
  expect_equal(defuzzify(B, "max", which = "last"), 7)
  expect_equal(defuzz_centroid(fs_triangular(0, 1, 2)), 1)
  expect_equal(defuzz_centroid(fs_rectangular(2, 4)), 3)
  expect_equal(defuzz_meanmax(fs_trapezoidal(0, 2, 6, 8)), 4)
})
