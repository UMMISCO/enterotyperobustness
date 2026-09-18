X <- simulate_counts(n = 80)
classes <- stats::setNames(c("A", "B")[attr(X, "group") %% 2 + 1], rownames(X))

test_that("erpm draws stratified subsets per fraction and is reproducible", {
  r <- erpm(X, classes, fractions = c(0.5, 0.8), n_iter = 2, marker_features = markers,
            verbose = FALSE)
  expect_s3_class(r, "erpm_result")
  expect_named(r$results, c("Q0.50", "Q0.80"))
  f <- r$results[["Q0.50"]]
  expect_length(f$index[[1]], 40)
  expect_equal(as.vector(table(classes[f$index[[1]]])), c(20, 20))
  expect_equal(sum(!is.na(f$assignments[, 1])), 40)
  expect_identical(f$index[[1]], stratified_subsample(classes, 0.5, seed = 1))
  expect_null(f$fits)
  r2 <- erpm(X, classes, fractions = 0.5, n_iter = 2, marker_features = markers, verbose = FALSE)
  expect_identical(f$assignments, r2$results[[1]]$assignments)
})

test_that("fits can be spilled to disk", {
  d <- withr::local_tempdir()
  r <- erpm(X, classes, fractions = 0.5, n_iter = 1, marker_features = markers,
            fits_dir = d, verbose = FALSE)
  expect_true(file.exists(file.path(d, "fits_Q0.50.rds")))
  expect_null(r$results[[1]]$fits)
})
