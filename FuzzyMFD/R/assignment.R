## Membership value assignment (Ross 2000, Section 4) -----------------------

## 4.2 Rank ordering ---------------------------------------------------------

rank_ordering <- function(P) {
  ## P[i, j] = number of respondents preferring item i to item j
  P <- as.matrix(P)
  if (nrow(P) != ncol(P) || nrow(P) < 2L) stop("'P' must be a square matrix")
  if (is.null(rownames(P))) rownames(P) <- colnames(P)
  if (is.null(rownames(P))) rownames(P) <- paste0("item", seq_len(nrow(P)))
  colnames(P) <- rownames(P)
  diag(P) <- NA
  if (any(P < 0, na.rm = TRUE)) stop("preference counts must be non-negative")
  pairs <- P + t(P)
  n_resp <- unique(round(pairs[upper.tri(pairs)], 8))
  if (length(n_resp) > 1L)
    warning("P[i, j] + P[j, i] is not constant: the matrix is not a ",
            "complete reciprocal (antisymmetric) comparison table")
  total <- rowSums(P, na.rm = TRUE)
  percent <- 100 * total / sum(total)
  tab <- data.frame(item = rownames(P), total = total,
                    percent = percent,
                    rank = rank(-total, ties.method = "min"),
                    membership = percent / max(percent),
                    row.names = NULL, stringsAsFactors = FALSE)
  out <- list(table = tab, preferences = P,
              responses = sum(total),
              membership = stats::setNames(tab$membership, tab$item))
  class(out) <- "rank_ordering"
  out
}

print.rank_ordering <- function(x, digits = 3, ...) {
  cat("Rank ordering from pairwise preferences (", x$responses,
      " responses)\n\n", sep = "")
  P <- x$preferences
  tab <- matrix(format(P), nrow(P), dimnames = dimnames(P))
  tab[is.na(P)] <- "-"
  tab <- cbind(tab, Total = x$table$total,
               Percent = format(round(x$table$percent, 1)), Rank = x$table$rank)
  print(tab, quote = FALSE, right = TRUE)
  cat("\nMembership in 'best' (percent / max percent), ascending:\n")
  m <- sort(x$membership)
  print(round(m, digits))
  invisible(x)
}

plot.rank_ordering <- function(x, ylab = "membership", xlab = "",
                               main = "Rank-ordering membership", ...) {
  m <- sort(x$membership)
  graphics::plot(seq_along(m), m, type = "h", lwd = 6, xaxt = "n",
                 ylim = c(0, 1), xlab = xlab, ylab = ylab, main = main, ...)
  graphics::axis(1, at = seq_along(m), labels = names(m))
  invisible(m)
}

## 4.3 Inductive reasoning: entropy minimisation ---------------------------

entropy_threshold <- function(x, class, candidates = NULL,
                              pure = c("middle", "entropy", "none"),
                              k = 1) {
  ## Equations (7)-(13): entropy S(x) = p(x) Sp(x) + q(x) Sq(x) of a split of
  ## the sample region at threshold x; samples with value <= x fall in the
  ## p region, the others in the q region.
  pure <- match.arg(pure)
  x <- as.numeric(x)
  if (length(x) != length(class)) stop("'x' and 'class' must have the same length")
  if (anyNA(x) || anyNA(class)) stop("'x' and 'class' must be free of NA")
  cls <- if (is.factor(class)) levels(droplevels(class)) else sort(unique(class))
  class <- as.character(class)
  cls <- as.character(cls)
  if (is.null(candidates)) {
    u <- sort(unique(x))
    candidates <- if (length(u) >= 2L) (u[-1L] + u[-length(u)]) / 2 else numeric(0)
  }
  candidates <- sort(as.numeric(candidates))
  n <- length(x)
  rows <- lapply(candidates, function(t) {
    left <- x <= t
    nl <- sum(left)
    nr <- n - nl
    pk <- vapply(cls, function(cl) (sum(left & class == cl) + 1) / (nl + 1), numeric(1))
    qk <- vapply(cls, function(cl) (sum(!left & class == cl) + 1) / (nr + 1), numeric(1))
    p <- nl / n
    Sp <- -k * sum(pk * log(pk))
    Sq <- -k * sum(qk * log(qk))
    c(x = t, stats::setNames(pk, paste0("p", seq_along(cls))),
      stats::setNames(qk, paste0("q", seq_along(cls))),
      p = p, q = 1 - p, Sp = Sp, Sq = Sq, S = p * Sp + (1 - p) * Sq)
  })
  tab <- if (length(rows)) as.data.frame(do.call(rbind, rows)) else NULL
  is_pure <- length(unique(class)) == 1L
  threshold <- NA_real_
  if (!is.null(tab) && nrow(tab)) {
    if (is_pure && pure == "middle") {
      threshold <- candidates[ceiling(length(candidates) / 2)]
    } else if (!(is_pure && pure == "none")) {
      threshold <- tab$x[which.min(tab$S)]
    }
  }
  out <- list(table = tab, threshold = threshold,
              S = if (is.na(threshold)) NA_real_ else tab$S[tab$x == threshold][1L],
              classes = cls, n = n, pure = is_pure, rule = pure)
  class(out) <- "entropy_threshold"
  out
}

print.entropy_threshold <- function(x, digits = 3, ...) {
  cat("Entropy-minimisation threshold (", x$n, " samples, classes: ",
      paste(x$classes, collapse = ", "), ")\n", sep = "")
  if (is.null(x$table)) {
    cat("no candidate thresholds\n")
    return(invisible(x))
  }
  print(round(x$table, digits), row.names = FALSE)
  cat("\nthreshold =", signif(x$threshold, digits + 2),
      "  S =", round(x$S, digits))
  if (x$pure) cat("  (pure region, rule = '", x$rule, "')", sep = "")
  cat("\n")
  invisible(x)
}

entropy_partition <- function(x, class, levels = 2, candidates = NULL,
                              pure = c("middle", "entropy", "none"), k = 1) {
  ## repeated partitioning (Section 4.3.1): the primary threshold (PRI) splits
  ## the region in two, secondary thresholds (SEC1, SEC2) split each side,
  ## tertiary thresholds (TER1..TER4) split each of the four parts, ...
  pure <- match.arg(pure)
  x <- as.numeric(x)
  class <- as.character(class)
  levels <- as.integer(levels)
  if (levels < 1L) stop("'levels' must be >= 1")
  prefix <- c("PRI", "SEC", "TER", "QUA", "QUI")
  regions <- list(rep(TRUE, length(x)))
  res <- list()
  steps <- list()
  for (lev in seq_len(levels)) {
    new_regions <- list()
    j <- 0L
    for (r in regions) {
      j <- j + 1L
      if (sum(r) < 2L) {
        new_regions <- c(new_regions, list(r & FALSE, r & FALSE))
        next
      }
      xr <- x[r]
      cand <- if (is.null(candidates)) NULL else
        candidates[candidates > min(xr) & candidates < max(xr)]
      et <- entropy_threshold(xr, class[r], candidates = cand, pure = pure, k = k)
      lab <- if (lev == 1L) "PRI" else
        paste0(if (lev <= length(prefix)) prefix[lev] else paste0("L", lev), j)
      if (!is.na(et$threshold)) {
        res[[length(res) + 1L]] <- data.frame(level = lev, label = lab,
                                              threshold = et$threshold, S = et$S,
                                              n = sum(r), stringsAsFactors = FALSE)
        steps[[lab]] <- et
        new_regions <- c(new_regions, list(r & x <= et$threshold, r & x > et$threshold))
      } else {
        new_regions <- c(new_regions, list(r & FALSE, r & FALSE))
      }
    }
    regions <- new_regions
  }
  tab <- if (length(res)) do.call(rbind, res) else
    data.frame(level = integer(0), label = character(0), threshold = numeric(0),
               S = numeric(0), n = integer(0))
  tab <- tab[order(tab$threshold), , drop = FALSE]
  rownames(tab) <- NULL
  out <- list(thresholds = tab, steps = steps)
  class(out) <- "entropy_partition"
  out
}

print.entropy_partition <- function(x, digits = 4, ...) {
  cat("Entropy-minimisation partition thresholds\n")
  tab <- x$thresholds
  tab$S <- round(tab$S, 3)
  tab$threshold <- signif(tab$threshold, digits)
  print(tab, row.names = FALSE)
  invisible(x)
}

partition_sets <- function(thresholds, range, labels = NULL, overlap = 0.5) {
  ## Membership functions whose toes are the thresholds (Figures 1.11-1.14):
  ## with k >= 2 sorted thresholds t1 < ... < tk on [lo, hi] the first set is
  ## 1 on [lo, t1] falling to 0 at t2, the middle sets are triangles with
  ## peak ti and feet t(i-1), t(i+1), and the last set rises from t(k-1) to 1
  ## at tk and stays 1 up to hi. A single threshold gives two sets crossing at
  ## membership 0.5 over the threshold.
  if (inherits(thresholds, "entropy_partition")) thresholds <- thresholds$thresholds$threshold
  t <- sort(unique(as.numeric(thresholds)))
  k <- length(t)
  if (!k) stop("no thresholds")
  lo <- range[1L]
  hi <- range[2L]
  if (!(lo <= t[1L] && t[k] <= hi && lo < hi)) stop("thresholds must lie inside 'range'")
  nset <- if (k == 1L) 2L else k
  if (is.null(labels)) labels <- .default_labels(nset)
  if (length(labels) != nset) stop("need ", nset, " labels")
  if (k == 1L) {
    if (!(overlap > 0 && overlap < 1)) stop("'overlap' must be in (0, 1)")
    w <- overlap * min(t - lo, hi - t)
    if (w <= 0) stop("the threshold must lie strictly inside 'range'")
    return(stats::setNames(list(
      fuzzy_set(unique(c(lo, t - w, t + w)), c(if (t - w > lo) 1, 1, 0), name = labels[1]),
      fuzzy_set(unique(c(t - w, t + w, hi)), c(0, 1, if (hi > t + w) 1), name = labels[2])),
      labels))
  }
  sets <- vector("list", k)
  sets[[1L]] <- if (t[1L] > lo) fuzzy_set(c(lo, t[1L], t[2L]), c(1, 1, 0), name = labels[1L]) else
    fuzzy_set(c(t[1L], t[2L]), c(1, 0), name = labels[1L])
  if (k > 2L) {
    for (i in 2:(k - 1L))
      sets[[i]] <- fs_triangular(t[i - 1L], t[i], t[i + 1L], name = labels[i])
  }
  sets[[k]] <- if (t[k] < hi) fuzzy_set(c(t[k - 1L], t[k], hi), c(0, 1, 1), name = labels[k]) else
    fuzzy_set(c(t[k - 1L], t[k]), c(0, 1), name = labels[k])
  stats::setNames(sets, labels)
}

inductive_mf <- function(x, class, levels = 2, range = base::range(x),
                         labels = NULL, candidates = NULL,
                         pure = c("middle", "entropy", "none"), overlap = 0.5) {
  pure <- match.arg(pure)
  part <- entropy_partition(x, class, levels = levels, candidates = candidates,
                            pure = pure)
  sets <- partition_sets(part, range = range, labels = labels, overlap = overlap)
  list(partition = part, sets = sets)
}

.default_labels <- function(n) {
  switch(as.character(n),
         "2" = c("NG", "PO"),
         "3" = c("NG", "ZE", "PO"),
         "5" = c("NB", "NS", "ZE", "PS", "PB"),
         "7" = c("NB", "NM", "NS", "ZE", "PS", "PM", "PB"),
         paste0("S", seq_len(n)))
}
