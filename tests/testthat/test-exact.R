test_that("exact cuts are nested and contain Monte Carlo draws", {
  x <- simulate_fsn(15, 0.3, kappa = 2, seed = 4)
  for (ix in c("in_strength", "in_closeness_min", "in_closeness_prod", "in_closeness_luk", "in_katz")) {
    fc <- fuzzy_centrality(x, ix, alpha = seq(0, 1, 0.1))
    mc <- mc_centrality(x, ix, B = 30, beta = fc$index$params$beta)
    expect_true(all(mc$lower >= fc$lower[, 1] - 1e-10 & mc$upper <= fc$upper[, 1] + 1e-10))
    expect_equal(as.numeric(fc$lower[, 1]), unclass(fc$triangle)$l)
    expect_equal(as.numeric(fc$upper[, 11]), unclass(fc$triangle)$m)
  }
})

test_that("strength is triangular and agrees with the IJAR degree", {
  A <- fsn_simulated()
  fc <- fuzzy_centrality(A, "in_strength", normalize = "n-1")
  expect_lt(max(triangular_approx(fc)$max_error), 1e-10)
  d <- ijar_degree(A, "in")
  expect_equal(unname(fc$lower[, 1]), d$l)
})

test_that("possibility, necessity, certification and top set (department A)", {
  skip_if_not(fsna_data_available(), "questionnaire data not available")
  fc <- fuzzy_centrality(fsn_department("A"), "in_strength")
  cmp <- possibility_matrix(fc)
  expect_true(cmp$exact)
  off <- row(cmp$necessity) != col(cmp$necessity)
  expect_true(all(cmp$necessity[off] <= cmp$possibility[off] + 1e-12))
  C <- certified_pairs(fc, 0)
  expect_true(all(cmp$necessity[C] >= 1 - 1e-12))
  ce <- certification(fc)
  expect_equal(ce$N1_not_separated, 0); expect_equal(ce$N1_not_separated_grid, 31)
  expect_equal(max(top_set(fc)), 1)
  ri <- rank_intervals(fc, 0); expect_true(all(ri$best <= ri$worst))
  for (m in c("cog", "possibility", "necessity", "certified", "top")) expect_equal(nrow(rank_fuzzy(fc, m)), 51)
})

test_that("non-isotone indices are detected", {
  expect_false(check_isotone(fsn_example(), function(W) -colSums(W), n_check = 10))
  expect_error(fuzzy_centrality(fsn_example(), function(W) -colSums(W)))
})

test_that("Table 4 of Porreca and Maturo (2026) is reproduced", {
  skip_if_not(fsna_data_available(), "questionnaire data not available")
  skip_on_cran()
  tab <- list(A = list(n = 51, resp = 24, in_strength = c(0.419, 0.675, 31, 3),
                       in_closeness_min = c(0.030, 0.068, 4, 43), in_closeness_prod = c(0.032, 0.256, 10, 24),
                       in_katz = c(0.158, 0.463, 20, 12)),
              B = list(n = 55, resp = 28, in_strength = c(0.680, 0.844, 31, 1),
                       in_closeness_min = c(0.155, 0.304, 24, 32), in_closeness_prod = c(0.064, 0.308, 18, 27),
                       in_katz = c(0.426, 0.678, 38, 8)))
  for (d in c("A", "B")) {
    x <- fsn_department(d)
    expect_equal(n_actors(x), tab[[d]]$n); expect_equal(sum(x$respondents), tab[[d]]$resp)
    for (ix in c("in_strength", "in_closeness_min", "in_closeness_prod", "in_katz")) {
      fc <- fuzzy_centrality(x, ix); np <- n_actors(x) * (n_actors(x) - 1) / 2
      got <- c(round(sum(certified_pairs(fc, 0)) / np, 3), round(sum(certified_pairs(fc, 0.5)) / np, 3),
               certification(fc)$N1_not_separated_grid, length(top_set(fc, alpha = 0.5)))
      expect_equal(certification(fc)$N1_not_separated, 0, info = paste(d, ix))
      expect_equal(got, tab[[d]][[ix]], info = paste(d, ix))
    }
  }
  mu <- top_set(fuzzy_centrality(fsn_department("A"), "in_strength"))
  expect_equal(round(unname(sort(mu, TRUE)[1:3]), 2), c(1, 0.87, 0.85))
})

test_that("bisection locates the suprema between grid levels", {
  x <- fsn_simulated(n = 20, respondents = 12, p = 0.4, seed = 5)
  fc <- fuzzy_centrality(x, "in_closeness_prod")
  Pb <- possibility_matrix(fc)$possibility
  Pu <- possibility_matrix(fc, "upper")$possibility
  Pg <- possibility_matrix(fc, "grid")$possibility
  expect_equal(possibility_matrix(fc)$rule, "bisection")
  expect_true(all(Pg <= Pb + 1e-12 & Pb <= Pu + 1e-12))
  # condition holds at the returned level and fails just above it (within tol)
  i <- which(Pb > 0 & Pb < 1, arr.ind = TRUE)[1, ]
  a <- Pb[i[1], i[2]]
  ct <- function(al) { w <- fsn_cut(x, al); list(lo = fc$index$fun(w$lower), up = fc$index$fun(w$upper)) }
  expect_gte(ct(a)$up[i[1]] - ct(a)$lo[i[2]], 0)
  expect_lt(ct(min(1, a + 2e-4))$up[i[1]] - ct(min(1, a + 2e-4))$lo[i[2]], 0)
  # strength: interpolation
  expect_equal(possibility_matrix(fuzzy_centrality(x, "in_strength"))$rule, "interpolate")
})
