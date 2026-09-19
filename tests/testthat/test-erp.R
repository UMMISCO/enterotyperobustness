X <- simulate_counts(n = 80)
div <- attr(X, "diversity")
anchors <- rownames(X)[1:40]
pool <- rownames(X)[41:80]

test_that("erp reproduces itself and structures its result", {
  r1 <- erp(X, anchors, pool, n_complementary = 12, n_iter = 3, marker_features = markers,
            diversity = div, verbose = FALSE)
  r2 <- erp(X, anchors, pool, n_complementary = 12, n_iter = 3, marker_features = markers,
            diversity = div, verbose = FALSE)
  expect_s3_class(r1, "erp_result")
  expect_identical(r1$assignments, r2$assignments)
  expect_identical(r1$index, r2$index)
  expect_equal(dim(r1$assignments), c(80, 3))
  expect_false(anyNA(r1$assignments[anchors, ]))
  expect_equal(sum(!is.na(r1$assignments[pool, 1])), 12)
  expect_identical(colnames(r1$posteriors[[1]]), enterotype_names())
  expect_equal(nrow(r1$posteriors[[1]]), 52)
  expect_equal(dim(r1$label_map), c(3, 4))
  expect_true(all(apply(r1$label_map, 1, function(m) setequal(m, 1:4))))
  expect_equal(nrow(r1$flags), 3); expect_equal(nrow(r1$gap), 3)
  expect_length(r1$fits, 3)
  # the draws are set.seed(i); sample.int(pool)
  set.seed(2); expect_identical(r1$index[[2]], sample.int(40, 12))
})

test_that("labels are the MAP component renamed, consistent with the posteriors", {
  r <- erp(X, anchors, pool, n_complementary = 12, n_iter = 2, marker_features = markers,
           verbose = FALSE)
  for (i in 1:2) {
    p <- r$posteriors[[i]]
    expect_identical(unname(r$assignments[rownames(p), i]),
                     colnames(p)[apply(p, 1, which.max)])
  }
})

test_that("n_cores does not change the result", {
  skip_on_os("windows")
  a <- erp(X, anchors, pool, 12, n_iter = 4, marker_features = markers, n_cores = 1, verbose = FALSE, keep_fits = FALSE)
  b <- erp(X, anchors, pool, 12, n_iter = 4, marker_features = markers, n_cores = 2, verbose = FALSE, keep_fits = FALSE)
  expect_identical(a$assignments, b$assignments)
  expect_equal(a$posteriors, b$posteriors)
})

test_that("with-replacement draws are handled and recorded", {
  r <- erp(X, anchors, pool[1:5], n_complementary = 12, n_iter = 2, marker_features = markers,
           verbose = FALSE, keep_fits = FALSE)
  expect_true(r$params$replace)
  expect_equal(nrow(r$posteriors[[1]]), 40 + length(unique(pool[1:5][r$index[[1]]])))
  expect_error(erp(X, anchors, pool[1:5], 12, n_iter = 2, marker_features = markers,
                   replace = "never", verbose = FALSE), "without replacement")
})

test_that("p_complementary + cohort_size gives the paper's counts", {
  r <- erp(X, anchors, pool, p_complementary = 0.3, cohort_size = 285, n_iter = 1,
           marker_features = markers, verbose = FALSE, keep_fits = FALSE, replace = "always")
  expect_equal(r$params$n_complementary, 86)
  expect_error(erp(X, anchors, pool, n_iter = 1, marker_features = markers), "n_complementary")
  expect_error(erp(X, anchors, c(pool, anchors[1]), 5, n_iter = 1, marker_features = markers), "both")
})

test_that("a diversity vector that misses a sample is an error, not a silent flag", {
  expect_error(erp(X, anchors, pool, 12, n_iter = 2, marker_features = markers,
                   diversity = div[-1], verbose = FALSE, keep_fits = FALSE), "no `diversity` value")
  fit <- fit_dmm(X[1:40, ], K = 4, seed = 1)
  expect_error(annotate_clusters(fit, X[1:40, ], markers, diversity = div[-1]), "no `diversity` value")
})

test_that("erp does not disturb the caller's random stream", {
  set.seed(77); ref <- stats::runif(2)
  set.seed(77)
  invisible(erp(X, anchors, pool, 12, n_iter = 2, marker_features = markers, verbose = FALSE, keep_fits = FALSE))
  expect_identical(stats::runif(2), ref)
})

test_that("seeds must be finite integers", {
  expect_error(erp(X, anchors, pool, 12, n_iter = 2, iteration_seeds = c(1, NA), marker_features = markers, verbose = FALSE), "iteration_seeds")
  expect_error(erp(X, anchors, pool, 12, n_iter = 2, dmm_seed = 1.5, marker_features = markers, verbose = FALSE), "dmm_seed")
})
