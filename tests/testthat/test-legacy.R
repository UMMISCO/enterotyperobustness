# The 2022 result files are private (MetaCardis); these tests run only where
# ENTEROTYPES_MIRROR points to the archive, and are skipped everywhere else.
root <- Sys.getenv("ENTEROTYPES_MIRROR", "")
exp_dir <- if (nzchar(root)) file.path(root, "analyses/2.Estimating_DMM/Article/Experiments") else ""
conf_dir <- if (nzchar(root)) file.path(root, "analyses/2.Estimating_DMM/Article/Confidence") else ""

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

test_that("legacy ERPM reader rebuilds one stored fraction with its diversity flags", {
  skip_if_not(dir.exists(conf_dir), "ENTEROTYPES mirror not available")
  f <- file.path(conf_dir, "Conf0.55_TotalEnvironment.RData")
  skip_if_not(file.exists(f), "Conf0.55 file not available")
  r <- read_legacy_erpm(f)
  expect_s3_class(r, "erpm_fraction")
  # the stored environments do not all carry `Percen`; script 00 maps files to fractions
  expect_true(is.na(r$fraction) || (r$fraction > 0 && r$fraction < 1))
  expect_equal(ncol(r$assignments), 500)
  expect_equal(length(r$index[[1]]) / length(r$samples$cohort), 0.5, tolerance = 0.01)   # 50 % draws, per-class rounding
  expect_equal(r$params$seed, 314L)
  expect_identical(colnames(r$posteriors[[1]]), enterotype_names())
  rb <- robustness(r)
  stored <- r$legacy_prob_enterotype
  expect_true(is.list(stored) || is.data.frame(stored) || is.matrix(stored))
  expect_output(print(r), "ERPM fraction")
})
