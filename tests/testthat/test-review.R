# Counterexamples from the independent technical review (12 Nov 2026)

test_that("betweenness shares tied best paths and is invariant to relabelling", {
  m <- matrix(0, 4, 4); m[1, 2:3] <- 0.8; m[2:3, 4] <- 0.8
  x <- fsn(M = m, labels = c("A", "B", "C", "D"))
  b <- ijar_betweenness(x, normalize = "none")
  expect_equal(b$paths_through, c(0, 0.5, 0.5, 0))
  b2 <- ijar_betweenness(x[c("A", "C", "B", "D")], normalize = "none")
  expect_equal(b2$paths_through[match(b$actor, b2$actor)], b$paths_through)
  # longer optimal routes and random relabelling
  set.seed(11)
  y <- simulate_fsn(12, 0.35, seed = 11)
  for (s in c(2, 3)) {
    b <- ijar_betweenness(y, max_steps = s, normalize = "none")
    p <- sample(n_actors(y)); yp <- y[p]
    bp <- ijar_betweenness(yp, max_steps = s, normalize = "none")
    expect_equal(bp$paths_through[match(b$actor, bp$actor)], b$paths_through)
  }
  expect_true(all(is.na(ijar_betweenness(fsn(M = matrix(c(0, .5, .5, 0), 2)), normalize = "pairs")$betweenness)))
})

test_that("split betweenness agrees with brute-force enumeration of optimal paths", {
  set.seed(5)
  for (rep in 1:5) {
    x <- simulate_fsn(7, 0.4, seed = rep)
    n <- n_actors(x); sc <- (x$L + x$M + x$R) / 3
    for (s in 2:3) {
      bf <- numeric(n)
      for (i in 1:n) for (j in 1:n) {
        if (i == j) next
        paths <- list(c(i, j))
        if (s >= 2) for (h in setdiff(1:n, c(i, j))) paths[[length(paths) + 1]] <- c(i, h, j)
        if (s >= 3) for (h in setdiff(1:n, c(i, j))) for (k in setdiff(1:n, c(i, j, h)))
          paths[[length(paths) + 1]] <- c(i, h, k, j)
        v <- sapply(paths, function(p) min(sc[cbind(p[-length(p)], p[-1])]))
        if (max(v) <= 0) next
        opt <- paths[v == max(v)]; len <- sapply(opt, length); opt <- opt[len == min(len)]
        for (p in opt) { mid <- p[-c(1, length(p))]; bf[mid] <- bf[mid] + 1 / length(opt) }
      }
      expect_equal(ijar_betweenness(x, max_steps = s, normalize = "none")$paths_through, bf)
    }
  }
})

test_that("invalid parameters are rejected", {
  l <- md <- rr <- matrix(0, 3, 3); md[1:2, 3] <- .5; rr[1:2, 3] <- 1
  expect_error(fuzzy_centrality(fsn(l, md, rr), "in_owa", w = c(-1, 2)), "non-negative")
  expect_error(centrality_owa(md, c(-1, 2)), "non-negative")
  expect_error(centrality_katz(fsn_example()$M, beta = -0.1), "non-negative")
  expect_error(tfn_score(tfn(.1, .2, .9), "liou_wang", lambda = 2), "lambda")
  expect_error(best_path_strength(fsn_example()$M, max_steps = 0), "positive integer")
  expect_error(ijar_best_paths(fsn_example(), max_steps = 0), "positive integer")
  expect_error(tfn(NA, 0.2, 0.3), "finite")
  m <- matrix(0, 5, 5); m[1:3, 1:3] <- .8; diag(m) <- 0
  expect_error(ijar_best_paths(fsn(M = m), max_steps = 4, tie_break = "indirect"), "max_steps <= 2")
})

test_that("undirected networks: comparisons are bounds and Monte Carlo draws are symmetric", {
  l <- matrix(c(0, .1, .1, 0), 2); md <- matrix(c(0, .5, .5, 0), 2); rr <- matrix(c(0, .9, .9, 0), 2)
  fc <- fuzzy_centrality(fsn(l, md, rr, directed = FALSE), "in_strength")
  expect_false(possibility_matrix(fc)$exact)
  W <- matrix(c(0, .3, .3, 0), 2); u <- fsn(W * .5, W, W * 2, directed = FALSE)
  mc <- mc_centrality(u, function(A) rep(max(abs(A - t(A))), 2), B = 20, keep_draws = TRUE)
  expect_true(all(mc$draws == 0))
})

test_that("unknown ties give interval cores and are not dropped as isolated", {
  fu <- fuzzy_centrality(set_unknown(fsn_simulated(), "nonrespondents"), "in_strength")
  expect_true(fu$interval_core)
  expect_output(print(fu), "cores are intervals")
  m <- matrix(0, 3, 3); m[1, 2] <- .5; U <- matrix(FALSE, 3, 3); U[3, 1] <- TRUE
  expect_equal(drop_isolated(fsn(M = m, unknown = U))$labels, c("1", "2", "3"))
  expect_warning(x <- fsn(M = matrix(c(0, NA, 1, 0), 2)), "treated as null")
  x <- fsn(M = matrix(c(0, NA, 1, 0), 2), na = "unknown")
  expect_true(x$unknown[2, 1])
  expect_error(fsn(M = matrix(c(0, NA, 1, 0), 2), na = "error"))
})

test_that("certified score is documented as heuristic; alpha = 0 certification agrees with CoG", {
  fc <- fuzzy_centrality(fsn_simulated(), "in_closeness_prod")
  C <- certified_pairs(fc, 0); cg <- tfn_score(fc$triangle, "cog")
  expect_true(all(outer(cg, cg, ">")[C]))
})

test_that("workflow with a single ranking method and participation records", {
  res <- fsna_analyze(fsn_example(), output_dir = NULL, indices = "in_strength", rank_methods = "core")
  expect_true(is.na(res$certification_summary$kendall_cog_necessity))
  d <- data.frame(rater = c("A", "B"), rated = c("B", ""), low = c(2, NA), best = c(4, NA), high = c(6, NA))
  x <- fsn_from_ratings(fuzzify_responses(d, "three_point", scale = c(0, 10)))
  expect_equal(sort(x$labels), c("A", "B"))
  expect_true(all(x$respondents))
  d2 <- data.frame(rater = c("A", "B"), rated = c("B", NA), label = c("strong", NA))
  expect_equal(nrow(fuzzify_responses(d2, "linguistic")), 2)
})

test_that("fsna_demo writes the report and one PDF page per figure", {
  d <- file.path(tempdir(), "demo_test")
  res <- fsna_demo(n = 6, open = FALSE, output_dir = d)
  expect_true(file.exists(res$files$report))
  expect_true(file.exists(file.path(d, "figures.pdf")))
  expect_true(any(grepl("Fig01b_glyphs", res$files$figures)))
})
