test_that("questionnaire data of the departments (private, run only when available)", {
  skip_if_not(fsna_data_available(), "questionnaire data not available")
  d <- fsna_questionnaire()
  expect_equal(sort(unique(d$department)), c("A", "B"))
  expect_equal(n_actors(fsn_department("A")), 51); expect_equal(n_actors(fsn_department("B")), 55)
  expect_equal(n_actors(fsn_department("B", members = "all")), 65)
})

test_that("missing questionnaire data give an informative error", {
  expect_error(fsna_questionnaire(file = tempfile()), "not distributed")
})

test_that("importers and simulated survey network", {
  x <- fsn_simulated()
  expect_equal(n_actors(x), 50); expect_equal(sum(x$respondents), 25)
  expect_true(all(x$M[!x$respondents, ] == 0))
  f <- system.file("extdata", "three_point_example.csv", package = "FSNA")
  d <- fuzzify_responses(read_ratings(f), "three_point", scale = c(0, 10))
  expect_true(all(d$l <= d$m & d$m <= d$r))
  x <- fsn_from_ratings(d)
  expect_s3_class(x, "fsn")
  s <- simulate_ratings(20, seed = 1)
  expect_s3_class(fsn_from_ratings(s, from = "rater", to = "rated"), "fsn")
  h <- make_collector(c("A", "B"), file = tempfile(fileext = ".html"))
  expect_true(any(grepl('"members":\\["A","B"\\]', readLines(h))))
})

test_that("raw questionnaire is parsed as in the original processing", {
  skip_if_not_installed("jsonlite")
  skip_if_not(fsna_data_available(), "questionnaire data not available")
  raw <- read_fuzzy_questionnaire(fsna_data_path())
  expect_equal(nrow(raw), 223 + 255)
  expect_equal(mean(round(raw$m_raw) == raw$button), 1, tolerance = 0.05)
})
