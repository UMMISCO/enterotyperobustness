test_that("plot_entropy and plot_alluvial build", {
  skip_if_not_installed("ggplot2"); skip_if_not_installed("ggalluvial")
  d <- expand.grid(sample = sprintf("s%02d", 1:10), condition = c("A", "B"), stringsAsFactors = FALSE)
  d$entropy <- stats::runif(nrow(d)); d$stratum <- rep(c("Bact1", "Bact2"), 10)
  d$enterotype <- d$stratum
  pw <- data.frame(stratum = c("Bact1", "Bact2"), cond1 = "A", cond2 = "B", p = c(0.01, 0.5), p_adj = c(0.02, 0.6))
  p <- plot_entropy(d, stratum = "stratum", pairwise = pw)
  expect_s3_class(p, "ggplot")
  expect_s3_class(plot_alluvial(d), "ggplot")
  expect_silent(invisible(ggplot2::ggplot_build(p)))
})
