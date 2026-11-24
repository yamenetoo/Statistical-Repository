test_that("rank ordering reproduces Table 1.1", {
  data(ross_colors, package = "FuzzyMFD", envir = environment())
  ro <- rank_ordering(ross_colors)
  expect_equal(ro$responses, 10000)
  expect_equal(ro$table$total, c(2248, 2377, 1782, 2087, 1506))
  expect_equal(round(ro$table$percent, 1), c(22.5, 23.8, 17.8, 20.9, 15.1))
  expect_equal(ro$table$rank, c(2, 1, 4, 3, 5))
  expect_equal(names(sort(ro$membership)),
               c("blue", "yellow", "green", "red", "orange"))
  expect_equal(max(ro$membership), 1)
})

test_that("entropy thresholds reproduce Tables 1.3-1.5", {
  data(ross_ellipse, package = "FuzzyMFD", envir = environment())
  x <- ross_ellipse$x
  cl <- ross_ellipse$class
  et <- entropy_threshold(x, cl, candidates = c(0.7, 1.0, 1.5, 3.45))
  tab <- et$table
  expect_equal(tab$p1, c(6/7, 7/8, 8/9, 8/10))
  expect_equal(tab$p2, c(2/7, 2/8, 2/9, 3/10))
  expect_equal(tab$q1, c(3/7, 2/6, 1/5, 1/4))
  expect_equal(tab$q2, c(5/7, 5/6, 1, 1))
  expect_equal(round(tab$Sp, 3), c(0.490, 0.463, 0.439, 0.540))
  expect_equal(round(tab$Sq, 2), c(0.60, 0.52, 0.32, 0.35))
  expect_equal(round(tab$S, 2), c(0.55, 0.49, 0.40, 0.49))
  expect_equal(et$threshold, 1.5)
  ## all mid-points give the same primary threshold
  expect_equal(entropy_threshold(x, cl)$threshold, 1.5)

  ## negative side (Table 1.4)
  neg <- x < 1.5
  e2 <- entropy_threshold(x[neg], cl[neg], candidates = c(0.175, 0.35, 0.7))
  expect_equal(round(e2$table$S, 3), c(0.454, 0.398, 0.459))
  expect_equal(e2$threshold, 0.35)

  ## positive side (Table 1.5): a single class
  e3 <- entropy_threshold(x[!neg], cl[!neg])
  expect_true(e3$pure)
  expect_equal(e3$threshold, 27.5)
  expect_equal(round(e3$table$S[e3$table$x == 27.5], 3), 0.366)

  part <- entropy_partition(x, cl, levels = 2)
  expect_equal(part$thresholds$threshold, c(0.35, 1.5, 27.5))
  expect_equal(part$thresholds$label, c("SEC1", "PRI", "SEC2"))
})

test_that("partition sets (Figure 1.14)", {
  sets <- partition_sets(c(0.35, 1.5, 27.5), range = c(0, 100))
  expect_named(sets, c("NG", "ZE", "PO"))
  expect_equal(membership(sets$NG, c(0, 0.35, 1.5)), c(1, 1, 0))
  expect_equal(membership(sets$ZE, c(0.35, 1.5, 27.5)), c(0, 1, 0))
  expect_equal(membership(sets$PO, c(1.5, 27.5, 100)), c(0, 1, 1))
  two <- partition_sets(1.5, range = c(0, 100))
  expect_equal(membership(two$NG, 1.5), 0.5)
  expect_equal(membership(two$PO, 1.5), 0.5)
  seven <- partition_sets(c(1, 2, 3, 4, 5, 6, 7), range = c(0, 8))
  expect_named(seven, c("NB", "NM", "NS", "ZE", "PS", "PM", "PB"))
  res <- inductive_mf(c(0, 0.1, 0.15, 0.2, 0.2, 0.5, 0.9, 1.1, 1.9, 5, 50, 100),
                      c(1, 1, 1, 1, 1, 2, 1, 1, 2, 2, 2, 2), levels = 2, range = c(0, 100))
  expect_length(res$sets, 3)
})
