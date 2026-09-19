set.seed(3)
n <- 30; conds <- c("A", "B", "C")
base <- stats::runif(n, 0.1, 0.6)
long <- do.call(rbind, lapply(conds, function(cc) data.frame(
  sample = sprintf("s%02d", 1:n), condition = cc,
  entropy = pmin(1, base + (cc == "C") * 0.2 + stats::rnorm(n, 0, 0.03)),
  stratum = rep(c("Bact1", "Bact2"), length.out = n), stringsAsFactors = FALSE)))

test_that("paired_entropy_test detects a within-sample shift and keeps the pairing", {
  skip_if_not_installed("PMCMRplus")
  r <- paired_entropy_test(long, stratum = "stratum")
  expect_equal(nrow(r$omnibus), 2)
  expect_true(all(r$omnibus$n == 15))
  expect_true(all(r$omnibus$p < 0.01))
  expect_equal(nrow(r$pairwise), 6)
  pc <- r$pairwise[r$pairwise$cond2 == "C" | r$pairwise$cond1 == "C", ]
  expect_true(all(pc$p_adj < 0.05))
  # a sample missing one condition is dropped from the block design
  r2 <- paired_entropy_test(long[-1, ], stratum = "stratum")
  expect_equal(r2$omnibus$n[r2$omnibus$stratum == "Bact1"], 14)
})

test_that("cochran_q_test is the binary Friedman test with McNemar post hoc", {
  lab <- cbind(A = rep(c("Bact1", "Bact2"), c(20, 10)), B = rep(c("Bact1", "Bact2"), c(20, 10)),
               C = rep(c("Bact1", "Bact2"), c(5, 25)))
  rownames(lab) <- sprintf("s%02d", 1:30)
  q <- cochran_q_test(lab, "Bact2")
  expect_lt(q$test$p, 0.001)
  expect_equal(unname(q$prevalence), c(1/3, 1/3, 25/30))
  expect_equal(q$pairwise$p[q$pairwise$cond1 == "A" & q$pairwise$cond2 == "B"], 1)
})

test_that("anchor_permutation_test keeps the row structure and finds a real effect", {
  lab <- cbind(A = rep(c("Bact1", "Bact2"), c(20, 10)), B = rep(c("Bact1", "Bact2"), c(20, 10)),
               C = rep(c("Bact1", "Bact2"), c(5, 25)))
  rownames(lab) <- sprintf("s%02d", 1:30)
  r <- anchor_permutation_test(lab, n_perm = 499)
  expect_equal(dim(r$table), c(3, 2))
  expect_equal(unname(r$table["C", "Bact2"]), 25)
  expect_lt(r$p, 0.05)
  expect_equal(nrow(r$cells), 6)
  # no effect: each sample identical across conditions -> statistic 0, p = 1
  same <- cbind(A = lab[, "A"], B = lab[, "A"], C = lab[, "A"])
  r0 <- anchor_permutation_test(same, n_perm = 99)
  expect_equal(r0$statistic, 0); expect_equal(r0$p, 1)
})

test_that("bootstrap_ci and normal_ci behave", {
  skip_if_not_installed("boot")
  x <- c(rep(1, 40), stats::runif(10, 0.5, 1))
  b <- bootstrap_ci(x, R = 500)
  expect_true(b[["lower"]] <= b[["mean"]] && b[["mean"]] <= b[["upper"]])
  expect_lte(b[["upper"]], 1)
  expect_equal(unname(bootstrap_ci(rep(1, 5))[1:3]), c(1, 1, 1))
  pc <- bootstrap_ci(x, R = 500, type = "percentile")
  expect_named(pc, c("mean", "lower", "upper", "n"))
  expect_true(pc[["lower"]] <= pc[["mean"]] && pc[["mean"]] <= pc[["upper"]])
  nc <- normal_ci(x); nn <- normal_ci(x, "n")
  expect_lt(nn[["upper"]] - nn[["lower"]], nc[["upper"]] - nc[["lower"]])
})

test_that("legacy_chisq matches chisq.test on the stacked table", {
  lab <- cbind(A = rep(c("Bact1", "Bact2"), c(20, 10)), C = rep(c("Bact1", "Bact2"), c(5, 25)))
  rownames(lab) <- sprintf("s%02d", 1:30)
  r <- legacy_chisq(lab)
  ref <- suppressWarnings(stats::chisq.test(table(rep(c("A", "C"), each = 30), as.vector(lab))))
  expect_equal(r$statistic, unname(ref$statistic))
  expect_equal(nrow(r$cells), 4)
})

test_that("paired_entropy_test refuses a stratum that varies within sample", {
  skip_if_not_installed("PMCMRplus")
  bad <- long
  bad$stratum[bad$condition == "C" & bad$sample %in% sprintf("s%02d", 1:10)] <- "Bact2"   # 10 samples change stratum in C
  expect_error(paired_entropy_test(bad, stratum = "stratum"), "varies within sample")
})

test_that("Nemenyi p-values are family-wise already and are not adjusted twice", {
  skip_if_not_installed("PMCMRplus")
  r <- paired_entropy_test(long, post_hoc = "nemenyi")
  expect_identical(r$pairwise$p_adj, r$pairwise$p)
  rc <- paired_entropy_test(long, post_hoc = "conover")
  expect_true(all(rc$pairwise$p_adj >= rc$pairwise$p))
})

test_that("cochran_q_test drops incomplete samples and reports the n it used", {
  lab <- cbind(A = rep(c("Bact1", "Bact2"), 10), B = rep(c("Bact1", "Bact2", "Bact2", "Bact1"), 5))
  rownames(lab) <- sprintf("s%02d", 1:20)
  lab2 <- lab; lab2[1, "B"] <- NA
  r <- cochran_q_test(lab2, "Bact2")
  expect_equal(r$test$n, 19)
  expect_equal(r$test$statistic, cochran_q_test(lab[-1, ], "Bact2")$test$statistic)
})

test_that("seeded functions leave the caller's random stream untouched", {
  skip_if_not_installed("boot")
  lab <- cbind(A = rep(c("Bact1", "Bact2", "Prevotella"), 10), B = rep(c("Bact2", "Bact1", "Prevotella"), 10))
  rownames(lab) <- sprintf("s%02d", 1:30)
  x <- stats::runif(20)
  set.seed(2024); ref <- stats::runif(3)
  set.seed(2024)
  invisible(bootstrap_ci(x, R = 200))
  invisible(anchor_permutation_test(lab, n_perm = 20))
  invisible(draw_complementary(letters, 5, seed = 1))
  invisible(stratified_subsample(stats::setNames(rep(c("a", "b"), 10), sprintf("s%02d", 1:20)), 0.5, seed = 1))
  expect_identical(stats::runif(3), ref)
})

test_that("degenerate inputs give defined results instead of NaN or opaque errors", {
  skip_if_not_installed("PMCMRplus")
  lab <- cbind(A = rep(c("Bact1", "Bact2"), 10), B = rep(c("Bact1", "Bact2"), 10)); rownames(lab) <- sprintf("s%02d", 1:20)
  r <- cochran_q_test(lab, "Prevotella")                    # absent everywhere
  expect_equal(r$test$statistic, 0); expect_equal(r$test$p, 1); expect_true(all(r$pairwise$p == 1))
  const <- long; const$entropy <- 0.3
  expect_message(rc <- paired_entropy_test(const), "identical")
  expect_equal(rc$omnibus$p, 1); expect_equal(nrow(rc$pairwise), 0)
})

test_that("anchor_permutation_test uses complete cases and reports n", {
  lab <- cbind(A = rep(c("Bact1", "Bact2", "Prevotella"), 10), B = rep(c("Bact2", "Bact1", "Prevotella"), 10), C = rep(c("Bact1", "Bact2", "Prevotella"), 10))
  rownames(lab) <- sprintf("s%02d", 1:30)
  lab2 <- lab; lab2[1, "C"] <- NA
  r <- anchor_permutation_test(lab2, n_perm = 50)
  expect_equal(r$n, 29); expect_equal(sum(r$table), 29 * 3)
  expect_equal(r$statistic, anchor_permutation_test(lab[-1, ], n_perm = 50)$statistic)
})
