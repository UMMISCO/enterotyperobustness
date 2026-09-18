root <- Sys.getenv("ENTEROTYPES_MIRROR", "/Users/jdz/PROGRAMMATION/R/ENTEROTYPES")
exp_dir <- file.path(root, "analyses/2.Estimating_DMM/Article/Experiments")

test_that("legacy ERP reader rebuilds the stored Healthy-Healthy labels from the fits", {
  skip_if_not(dir.exists(exp_dir), "ENTEROTYPES mirror not available")
  r <- read_legacy_erp(exp_dir, "Healthy-Healthy", keep_fits = FALSE)
  expect_s3_class(r, "erp_result")
  expect_length(r$samples$anchors, 142)
  expect_equal(ncol(r$assignments), 500)
  expect_equal(r$legacy_agreement, 1)
  expect_identical(colnames(r$posteriors[[1]]), enterotype_names())
  expect_equal(sum(r$flags$bact2_not_lowest_diversity), 0)
  expect_equal(sum(r$flags$rum_not_highest_diversity), 3)   # Supplementary Table S1 (stored 2022 flags)
})
