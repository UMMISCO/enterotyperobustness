set.seed(11)
n <- 40; conds <- c("Own", "B", "C")
base <- stats::runif(n, 0.05, 0.5)
lab0 <- rep(c("Bact1", "Bact2"), length.out = n)
long <- do.call(rbind, lapply(conds, function(cc) data.frame(
  sample = sprintf("s%02d", 1:n), condition = cc,
  entropy = pmin(1, base + (cc == "C") * 0.25 + stats::rnorm(n, 0, 0.03)),
  label = lab0, stringsAsFactors = FALSE)))

test_that("mixed_entropy_model recovers a within-individual shift and returns tidy tables", {
  skip_if_not_installed("lme4"); skip_if_not_installed("lmerTest"); skip_if_not_installed("emmeans")
  r <- mixed_entropy_model(long, label = "label")
  expect_s3_class(r, "mixed_entropy_model")
  expect_true(r$converged)
  expect_equal(r$n_samples, 40)
  expect_true("condition" %in% r$omnibus$term)
  expect_lt(r$omnibus$p[r$omnibus$term == "condition"], 1e-6)
  expect_equal(nrow(r$pairwise), 6)                       # 3 pairs x 2 labels
  expect_setequal(unique(r$pairwise$stratum), c("Bact1", "Bact2"))
  pc <- r$pairwise[r$pairwise$cond1 == "C" | r$pairwise$cond2 == "C", ]
  expect_true(all(pc$p_adj < 0.01))
  expect_true(all(c("stratum", "cond1", "cond2", "p", "p_adj") %in% names(r$pairwise)))
  expect_gt(r$variance$icc, 0.3)                          # base varies by individual
  expect_true(!is.null(r$omnibus_rank))
  r0 <- mixed_entropy_model(long, rank_check = FALSE)
  expect_equal(unique(r0$pairwise$stratum), "all")
  expect_null(r0$omnibus_rank)
})

test_that("condition names containing hyphens survive the contrast split", {
  skip_if_not_installed("emmeans")
  l2 <- long; l2$condition <- c(Own = "Healthy-Healthy", B = "Healthy-Obese", C = "Healthy-Diabetes")[l2$condition]
  r <- mixed_entropy_model(l2, rank_check = FALSE)
  expect_setequal(unique(c(r$pairwise$cond1, r$pairwise$cond2)),
                  c("Healthy-Healthy", "Healthy-Obese", "Healthy-Diabetes"))
})

lab <- cbind(Own = rep(c("Bact1", "Bact2", "Prevotella"), c(20, 12, 8)),
             B = rep(c("Bact1", "Bact2", "Prevotella"), c(20, 12, 8)),
             C = rep(c("Bact1", "Bact2", "Prevotella"), c(30, 4, 6)))
rownames(lab) <- sprintf("s%02d", 1:40)
set.seed(5); flip <- sample(40, 6); lab[flip, "B"] <- ifelse(lab[flip, "B"] == "Bact1", "Bact2", "Bact1")

test_that("multinomial_mixed_test fits or fails gracefully", {
  skip_if_not_installed("mclogit")
  r <- multinomial_mixed_test(lab, reference_condition = "Own")
  expect_true(is.logical(r$converged))
  if (r$converged) {
    expect_true(is.finite(r$wald$statistic))
    expect_equal(r$wald$df, 4)                             # 2 non-reference conditions x 2 non-reference enterotypes
    expect_setequal(r$per_enterotype$enterotype, c("Bact2", "Prevotella"))
    expect_equal(nrow(r$coefficients), 4)
  } else {
    # a fit that ran but did not converge keeps its tables for diagnosis and says why
    expect_true(length(r$messages) > 0)
  }
  same <- cbind(Own = lab[, "Own"], B = lab[, "Own"], C = lab[, "Own"])
  r0 <- multinomial_mixed_test(same, reference_condition = "Own")
  expect_true(is.logical(r0$converged))                    # degenerate design: no error either way
})

test_that("multinomial_permutation_test returns a valid p-value and detects the shift", {
  skip_if_not_installed("nnet")
  r <- multinomial_permutation_test(lab, reference_condition = "Own", n_perm = 199)
  expect_true(is.finite(r$statistic)); expect_equal(r$df, 4)
  expect_true(r$p > 0 && r$p <= 1)
  expect_lt(r$p, 0.05)
})

test_that("multinomial_permutation_test handles a two-label design (binary logit)", {
  skip_if_not_installed("nnet")
  set.seed(5)
  two <- cbind(Own = rep(c("Bact1", "Bact2"), 20),
               B = rep(c("Bact1", "Bact2"), 20),
               C = sample(rep(c("Bact1", "Bact2"), 20)))
  rownames(two) <- sprintf("s%02d", 1:40)
  r <- multinomial_permutation_test(two, reference_condition = "Own", n_perm = 99)
  expect_true(is.finite(r$statistic))
  expect_equal(r$df, 2)                                   # 2 non-reference conditions x 1 logit
  expect_true(r$p > 0 && r$p <= 1)
  expect_length(r$coefficients, 3)
  one <- two; one[] <- "Bact1"
  expect_error(multinomial_permutation_test(one, reference_condition = "Own", n_perm = 9), "fewer than two")
})

test_that("mixed_entropy_model uses Satterthwaite df for the contrasts, whatever emmeans' default", {
  skip_if_not_installed("lme4"); skip_if_not_installed("lmerTest"); skip_if_not_installed("emmeans")
  m <- mixed_entropy_model(long, label = "label")
  d <- data.frame(.y = long$entropy, .cond = factor(long$condition), .lab = factor(long$label), .id = factor(long$sample))
  fit <- lmerTest::lmer(.y ~ .cond * .lab + (1 | .id), data = d, REML = TRUE)
  ct <- as.data.frame(emmeans::contrast(emmeans::emmeans(fit, ~ .cond | .lab, lmer.df = "satterthwaite"), "pairwise", adjust = "none"))
  expect_equal(m$pairwise$df, ct$df, tolerance = 1e-8)
  expect_equal(m$pairwise$p, ct$p.value, tolerance = 1e-8)
})

test_that("multinomial_permutation_test drops incomplete samples and counts valid permutations", {
  skip_if_not_installed("nnet")
  lab2 <- lab; lab2[1, "B"] <- NA
  r <- multinomial_permutation_test(lab2, reference_condition = "Own", n_perm = 49)
  expect_equal(r$n, nrow(lab) - 1); expect_equal(r$n_perm, 49)
  expect_true(r$p >= 1 / 50 && r$p <= 1)
})
