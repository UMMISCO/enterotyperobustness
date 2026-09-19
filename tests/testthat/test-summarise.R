X <- simulate_counts(n = 60)
anchors <- rownames(X)[1:30]; pool <- rownames(X)[31:60]
r <- erp(X, anchors, pool, 10, n_iter = 4, marker_features = markers, verbose = FALSE, keep_fits = FALSE)

test_that("assignment_frequency uses the number of iterations drawn as denominator", {
  fr <- assignment_frequency(r)
  expect_equal(dim(fr), c(60, 4))
  expect_equal(unname(rowSums(fr[anchors, ])), rep(1, 30))
  nd <- attr(fr, "n_drawn")
  expect_true(all(nd[anchors] == 4))
  never <- names(nd)[nd == 0]
  if (length(never)) expect_true(all(is.na(fr[never, ])))
  # by hand for one anchor
  a <- r$assignments[anchors[1], ]
  expect_equal(unname(fr[anchors[1], "Bact2"]), mean(a == "Bact2"))
})

test_that("entropy is 0 for unanimous and 1 for uniform", {
  fr <- rbind(c(1, 0, 0, 0), c(0.25, 0.25, 0.25, 0.25))
  colnames(fr) <- enterotype_names()
  expect_equal(unname(assignment_entropy(fr)), c(0, 1))
})

test_that("robustness is the max frequency with ties counted", {
  rb <- robustness(r)
  fr <- assignment_frequency(r)
  expect_equal(rb$robustness[match(anchors, rb$sample)], unname(apply(fr[anchors, ], 1, max)))
  expect_true(is.logical(rb$tie))
})

test_that("mean_posterior aligns columns by name", {
  mp <- mean_posterior(r)
  expect_identical(colnames(mp), enterotype_names())
  expect_equal(unname(rowSums(mp[anchors, ])), rep(1, 30))
})

test_that("core_variable implements the published unanimity rule", {
  lab <- rbind(a = c("Bact1", "Bact1", "Bact1"), b = c("Bact1", "Bact2", "Bact1"),
               c = c("Prevotella", "Prevotella", "Ruminococcus"))
  cv <- core_variable(lab, "unanimous")
  expect_equal(cv$status[cv$sample == "a" & cv$enterotype == "Bact1"], "core")
  expect_equal(sort(cv$status[cv$sample == "b"]), c("variable", "variable"))   # variable for Bact1 and Bact2
  expect_false("a" %in% cv$sample[cv$enterotype == "Bact2"])
  cv1 <- core_variable(lab, "k_minus_1")
  expect_equal(cv1$status[cv1$sample == "b" & cv1$enterotype == "Bact1"], "core")
  freqs <- list(matrix(c(1,0,0,0, 1,0,0,0, 0,0,1,0), 3, 4, byrow = TRUE, dimnames = list(c("a","b","c"), enterotype_names())),
                matrix(c(1,0,0,0, 0.4,0.6,0,0, 0,0,0.5,0.5), 3, 4, byrow = TRUE, dimnames = list(c("a","b","c"), enterotype_names())),
                matrix(c(1,0,0,0, 0.9,0.1,0,0, 0,0,0.6,0.4), 3, 4, byrow = TRUE, dimnames = list(c("a","b","c"), enterotype_names())))
  sd_rule <- core_variable(lab, "sd", frequencies = freqs, threshold = 0.2)
  expect_equal(sd_rule$status, c("core", "variable", "variable"))
})

test_that("merge_pairs: A keeps the 4-way argmax, B re-assigns on the sums", {
  # x: argmax Bact1, and Bact1+Bact2 is also the larger sum -> A and B agree
  # y: argmax Bact1 (0.4) but Prevotella+Ruminococcus = 0.55 > 0.45 -> A and B differ
  fr <- rbind(x = c(0.5, 0.1, 0.4, 0), y = c(0.4, 0.05, 0.35, 0.2))
  colnames(fr) <- enterotype_names()
  A <- merge_pairs(fr, method = "A"); B <- merge_pairs(fr, method = "B")
  expect_equal(A$pair, c("Bact1Bact2", "Bact1Bact2"))
  expect_equal(A$value, c(0.6, 0.45))
  expect_equal(B$pair, c("Bact1Bact2", "PrevotellaRuminococcus"))
  expect_equal(B$value, c(0.6, 0.55))
  other <- list(c("Bact2", "Ruminococcus"), c("Bact1", "Prevotella"))
  A2 <- merge_pairs(fr, other, "A"); B2 <- merge_pairs(fr, other, "B")
  expect_equal(A2$pair[1], "Bact1Prevotella"); expect_equal(B2$pair[1], "Bact1Prevotella")
})

test_that("labels_by_condition returns anchors x conditions", {
  r2 <- erp(X, anchors, pool, 10, n_iter = 4, marker_features = markers, iteration_seeds = 5:8,
            verbose = FALSE, keep_fits = FALSE)
  lab <- labels_by_condition(list(c1 = r, c2 = r2))
  expect_equal(dim(lab), c(30, 2))
  expect_true(all(lab %in% enterotype_names()))
})

test_that("core_variable drops incomplete rows instead of calling them variable", {
  lab <- rbind(a = c("Bact1", "Bact1", NA), b = c("Bact1", "Bact1", "Bact1"))
  expect_message(cv <- core_variable(lab, "unanimous"), "dropped")
  expect_false("a" %in% cv$sample)
  expect_equal(cv$status[cv$sample == "b"], "core")
})

test_that("mean_posterior's sd attribute is the sample standard deviation", {
  mp <- mean_posterior(r)
  s <- attr(mp, "sd")
  a <- anchors[1]
  vals <- vapply(r$posteriors, function(p) p[a, "Bact1"], numeric(1))
  expect_equal(unname(s[a, "Bact1"]), stats::sd(vals))
  expect_true(all(is.na(s[attr(assignment_frequency(r), "n_drawn") < 2, ])))
})

test_that("assignment_entropy ignores extra columns and rejects non-frequencies", {
  fr <- assignment_frequency(r)[anchors, ]
  e <- assignment_entropy(fr)
  expect_identical(assignment_entropy(cbind(fr, n_drawn = 4)), e)
  expect_error(assignment_entropy(unname(cbind(fr, 7))), "frequencies")   # no enterotype names: every column is checked
  expect_true(all(e >= 0 & e <= 1))
})

test_that("assignment_frequency refuses a whole ERPM object with a helpful message", {
  fake <- structure(list(results = list()), class = "erpm_result")
  expect_error(assignment_frequency(fake), "one fraction")
})

test_that("robustness defaults to the seed recorded in the result", {
  expect_identical(robustness(r), robustness(r, seed = r$params$seed))
  expect_equal(enterotyperobustness:::result_seed(r), 314L)
})

test_that("input validation: entropy rows must sum to 1, pairs must be disjoint, core_variable may be empty", {
  fr <- assignment_frequency(r)[anchors, ]
  bad <- fr; bad[1, "Bact1"] <- bad[1, "Bact1"] + 0.6
  expect_error(assignment_entropy(bad), "sum to 1")
  expect_error(merge_pairs(fr, pairs = list(c("Bact1", "Bact2"), c("Bact2", "Prevotella", "Ruminococcus"))), "disjoint")
  lab <- rbind(a = c("Bact1", NA)); 
  expect_message(cv <- core_variable(lab), "dropped")
  expect_s3_class(cv, "data.frame"); expect_equal(nrow(cv), 0)
})
