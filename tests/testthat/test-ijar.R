test_that("best paths prefer strong indirect routes", {
  x <- fsn_example()
  bp <- ijar_best_paths(x, 3)
  expect_equal(best_path_route(bp, "A", "D"), c("A", "B", "D"))
  expect_true(all(bp$score >= 0))
})
