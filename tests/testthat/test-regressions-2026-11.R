# Regression tests for the defects found in the review of 20 November 2026

test_that("selection and permutation of actors keep the whole-network context", {
  f <- fuzzy_centrality(fsn_example(), "in_closeness_prod")
  p <- possibility_matrix(f)$possibility
  ids <- c(2, 1, 3:8)
  expect_equal(unname(possibility_matrix(f[ids])$possibility), unname(p[ids, ids]))
  expect_equal(unname(possibility_matrix(f[1:3])$possibility), unname(p[1:3, 1:3]))
  z <- fuzzy_centrality(fsn_example(), "in_closeness_prod")        # no cache yet
  expect_equal(unname(possibility_matrix(z[1:3])$possibility), unname(p[1:3, 1:3]))
  z <- fuzzy_centrality(fsn_example(), "in_closeness_prod")
  expect_equal(unname(possibility_matrix(z[ids])$possibility), unname(p[ids, ids]))
  z <- fuzzy_centrality(fsn_example(), "in_closeness_prod")
  zz <- z[c("C", "A", "F")][2:3]                                    # nested, by name
  expect_equal(unname(possibility_matrix(zz)$possibility), unname(p[c("A", "F"), c("A", "F")]))
  expect_equal(possibility_matrix(f)$possibility, p)                # parent unchanged
  expect_equal(unname(possibility_matrix(f[rep(TRUE, 8)])$possibility), unname(p))
  expect_error(f["nobody"], "unknown actor")
  expect_error(f[c(TRUE, FALSE)], "one value per actor")
  expect_error(f[20], "out of range")
})

test_that("as_igraph exports one edge per unordered pair in undirected networks", {
  skip_if_not_installed("igraph")
  x <- fsn(M = matrix(c(0, .5, .5, 0), 2), directed = FALSE)
  g <- as_igraph(x)
  expect_equal(igraph::ecount(g), 1)
  expect_equal(unname(igraph::strength(g)), c(.5, .5))
  expect_false(any(igraph::which_multiple(g)))
  y <- fsn_example(); g2 <- as_igraph(y)
  expect_equal(igraph::ecount(g2), nrow(fsn_ties(y)))               # directed unchanged
})

test_that("plot_fuzzy_scatter matches the axes by actor", {
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off())
  x <- tfn(c(.1, .4), c(.2, .5), c(.3, .6), names = c("A", "B"))
  d <- plot_fuzzy_scatter(x, x[2:1])
  expect_equal(d$x, d$y)
  expect_error(plot_fuzzy_scatter(x, tfn(.1, .2, .3, names = "A")), "same actors")
  expect_error(plot_fuzzy_scatter(x, x, colour = c(B = 1)), "does not cover")
})

test_that("numerical brackets of possibility and necessity contain the exact values", {
  l <- m <- r <- matrix(0, 3, 3)
  m[3, 1] <- .5; r[3, 1] <- 1
  l[3, 2] <- m[3, 2] <- r[3, 2] <- 1 - .5e-5
  z <- fuzzy_centrality(fsn(l, m, r), "in_owa", w = c(1, 0))
  res <- possibility_matrix(z)
  # exact Pi(1 >= 2) = 1e-5, N(2 > 1) = 0.99999
  expect_lte(res$possibility_lower[1, 2], 1e-5)
  expect_gte(res$possibility_upper[1, 2], 1e-5)
  expect_lte(res$necessity_lower[2, 1], 0.99999)
  expect_gte(res$necessity_upper[2, 1], 0.99999)
  expect_lt(res$necessity_lower[2, 1], 1)                            # not a guaranteed N = 1
  expect_error(possibility_matrix(z, tol = 0), "positive")
  expect_error(possibility_matrix(z, tol = Inf), "positive")
  expect_silent(possibility_matrix(z, tol = 1e-300))                # terminates
  # certification at alpha = 0 implies a guaranteed N = 1 (not the converse)
  fc <- fuzzy_centrality(fsn_simulated(n = 15, respondents = 8, seed = 3), "in_closeness_min")
  pm <- possibility_matrix(fc); C0 <- certified_pairs(fc, 0)
  off <- row(C0) != col(C0)
  expect_true(all(pm$necessity_lower[off][C0[off]] == 1))
})

test_that("guaranteed N = 1 with touching supports, and exact comparisons with one", {
  # supports [0, 0.5] and [0.5, 1] touch at 0.5 with membership 0: exact N = 1,
  # cuts separated at every alpha > 0, closed supports not strictly separated
  l <- m <- r <- matrix(0, 3, 3)
  l[3, 1] <- 0; m[3, 1] <- .2; r[3, 1] <- .5
  l[3, 2] <- .5; m[3, 2] <- .8; r[3, 2] <- 1
  f <- fuzzy_centrality(fsn(l, m, r), "in_strength")
  expect_equal(possibility_matrix(f)$necessity_lower[2, 1], 1)
  expect_false(certified_pairs(f, 0)[2, 1])
  # two such pairs: actors 2 > 1 (supports [0.25, 0.5] and [0, 0.25]) and
  # 1 > 3 (support [0, 0.25] and the crisp value 0 of actor 3)
  expect_equal(certification(f)$N1_not_separated_guaranteed, 2)
  # true Pi = 1e-13: N is not 1 and must not be counted as guaranteed
  l <- m <- r <- matrix(0, 3, 3); m[3, 1] <- .5; r[3, 1] <- 1
  l[3, 2] <- m[3, 2] <- r[3, 2] <- 1 - .5e-13
  pm <- possibility_matrix(fuzzy_centrality(fsn(l, m, r), "in_strength"))
  expect_lt(pm$necessity_lower[2, 1], 1)
  rf <- rank_fuzzy(fuzzy_centrality(fsn(l, m, r), "in_strength"), "necessity")
  expect_true(is.data.frame(rf))
})

test_that("interpolation is refused for nonlinear indices", {
  L <- M <- R <- matrix(0, 4, 4)
  L[3, 1] <- 0; M[3, 1] <- .1; R[3, 1] <- .9
  L[4, 1] <- .5; M[4, 1] <- .6; R[4, 1] <- .7
  L[3, 2] <- M[3, 2] <- R[3, 2] <- .65
  f <- fuzzy_centrality(fsn(L, M, R), "in_owa", w = c(1, 0, 0), alpha = c(0, 1))
  expect_error(possibility_matrix(f, "interpolate"), "affine")
  b <- possibility_matrix(f, "bisection", tol = 1e-10)
  expect_lte(b$possibility_lower[1, 2], 0.5); expect_gte(b$possibility_upper[1, 2], 0.5)
  expect_equal(b$possibility[1, 2], 0.5, tolerance = 1e-9)
  expect_equal(possibility_matrix(f)$rule, "bisection")
})

test_that("scatter of two vectors without identifiers", {
  grDevices::pdf(file = tempfile(fileext = ".pdf")); on.exit(grDevices::dev.off())
  x <- tfn(c(.1, .2), c(.2, .3), c(.4, .5))
  expect_warning(plot_fuzzy_scatter(x, x), "by position")
  d <- suppressWarnings(plot_fuzzy_scatter(x, x))
  expect_equal(d$actor, c("1", "2")); expect_equal(d$x, d$y)
})

test_that("simulation_study computes the indicators at alpha = 0.5", {
  sim <- simulation_study(n = 5, p = .5, kappa = 2, nrep = 1, indices = "in_strength",
                          B = 2, alpha = c(0, 1), seed = 4, progress = FALSE)
  set.seed(4 * 100003L + 1)
  x <- simulate_fsn(5, .5, 2, model = "logit", sigma = c(.05, .20))
  z <- fuzzy_centrality(x, "in_strength", alpha = c(0, .5, 1))
  expect_equal(sim$cert05, sum(certified_pairs(z, .5)) / choose(5, 2))
  expect_equal(sim$alpha_eval, 0.5)
})

test_that("input validation and edge cases", {
  m <- matrix(0, 2, 2); u <- matrix(FALSE, 2, 2); u[1, 2] <- TRUE
  expect_error(fsn(M = m, directed = FALSE, unknown = u), "symmetric")
  expect_error(fsn(M = m, unknown = u, unknown_range = c(0.5, 0.2)), "unknown_range")
  expect_s3_class(fsn_simulated(n = 10, respondents = 5, sigma = c(.01, .02)), "fsn")
  z <- simulate_ratings(n = 8, rate_prob = 0, seed = 1)
  expect_true(is.data.frame(z)); expect_equal(nrow(z), 0)
  expect_true(all(c("rater", "rated", "l", "m", "r") %in% names(z)))
  ci <- check_isotone(fsn(M = matrix(0, 3, 3)), "in_strength")
  expect_true(ci); expect_equal(attr(ci, "n_checked"), 0L)
  s <- simulate_fsn(6, seed = 1)
  expect_equal(dim(s[1:3]$latent), c(3L, 3L))
})

test_that("fsna_analyze handles max_steps = Inf, empty networks and closes devices", {
  before <- grDevices::dev.list()
  expect_type(fsna_analyze(fsn_example(), indices = "in_strength", max_steps = Inf,
                           output_dir = tempfile("fsna_inf_")), "list")
  suppressWarnings(r <- fsna_analyze(fsn(M = matrix(0, 3, 3)), indices = "in_strength",
                                     output_dir = tempfile("fsna_zero_")))
  expect_type(r, "list")
  expect_identical(grDevices::dev.list(), before)
})

test_that("the questionnaire cache is invalidated when the file changes", {
  skip_if_not(fsna_data_available(), "questionnaire data not available")
  path <- tempfile(fileext = ".csv")
  file.copy(fsna_data_path(), path)
  first <- fsna_questionnaire(path)
  raw <- utils::read.csv(path, sep = ";", check.names = FALSE)
  raw <- raw[raw$dipartimento == 0, , drop = FALSE]
  utils::write.table(raw, path, sep = ";", row.names = FALSE, qmethod = "double")
  second <- fsna_questionnaire(path)
  expect_true(all(c("A", "B") %in% first$department))
  expect_true(nrow(second) > 0 && all(second$department == "B"))
})
