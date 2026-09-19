X <- simulate_counts()
div <- attr(X, "diversity")
fit <- fit_dmm(X, K = 4, seed = 314)

test_that("fit_dmm returns a DMN with K components and named accessors", {
  expect_s4_class(fit, "DMN")
  expect_equal(dim(dmm_alpha(fit)), c(ncol(X), 4))
  expect_equal(dim(dmm_posterior(fit)), c(nrow(X), 4))
  expect_equal(colSums(dmm_centroids(fit)), rep(1, 4), ignore_attr = TRUE)
  expect_length(dmm_assign(fit), nrow(X))
})

test_that("rank annotation recovers the simulated structure and flags nothing", {
  ann <- annotate_clusters(fit, X, markers, diversity = div, method = "rank",
                           extra_checks = list(FS = c("Faecalibacterium", "Subdoligranulum")))
  expect_s3_class(ann, "enterotype_annotation")
  expect_setequal(names(ann$map), enterotype_names())
  expect_setequal(ann$map, 1:4)
  expect_identical(ann$component_names[ann$map], names(ann$map))
  # simulated group 1 is Bact2, 3 is Prevotella, 4 is Ruminococcus
  grp <- attr(X, "group"); comp <- dmm_assign(fit)
  majority <- function(g) as.integer(names(which.max(table(comp[grp == g]))))
  expect_equal(unname(ann$map["Bact2"]), majority(1))
  expect_equal(unname(ann$map["Prevotella"]), majority(3))
  expect_equal(unname(ann$map["Ruminococcus"]), majority(4))
  expect_false(ann$flags$bact2_not_lowest_diversity)
  expect_false(ann$flags$rum_not_highest_diversity)
  expect_false(ann$flags$rum_not_highest_FS)
  expect_gt(ann$gap$bact_gap, 0)
  expect_false(ann$gap$prevotella_collision)
})

test_that("match annotation agrees with rank on a clean fit and reports margins", {
  ref <- reference_centroids(X, markers, div, fit = fit)
  expect_identical(colnames(ref), enterotype_names())
  ann_m <- annotate_clusters(fit, X, markers, div, method = "match", reference = ref)
  ann_r <- annotate_clusters(fit, X, markers, div, method = "rank")
  expect_identical(ann_m$map, ann_r$map)
  expect_equal(ann_m$gap$matching_cost, 0, tolerance = 1e-12)
  expect_gt(ann_m$gap$matching_margin_min, 0)
  # a different fit (other seed) still maps onto the reference names
  fit2 <- fit_dmm(X, K = 4, seed = 99)
  ann2 <- annotate_clusters(fit2, X, markers, div, method = "match", reference = ref)
  expect_setequal(ann2$map, 1:4)
})

test_that("annotation refuses wrong inputs", {
  expect_error(annotate_clusters(fit, X, c(bacteroides = "Bacteroides")), "prevotella")
  expect_error(annotate_clusters(fit, X, c(bacteroides = "Nope", prevotella = "Prevotella")), "not in")
  expect_error(annotate_clusters(fit, X, markers, method = "match"), "reference")
})

test_that("annotate_clusters refuses a count table that is not the fitted table in the fitted order", {
  expect_error(annotate_clusters(fit, X[rev(seq_len(nrow(X))), ], markers), "fitted order")
  expect_error(annotate_clusters(fit, X[, rev(seq_len(ncol(X)))], markers), "features of `fit`")
  expect_s3_class(annotate_clusters(fit, X, markers), "enterotype_annotation")
})
