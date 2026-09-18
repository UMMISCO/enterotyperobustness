test_that("integerise_counts truncates like DirichletMultinomial 1.28 did", {
  x <- matrix(c(0.2, 1.7, 3491.7, 2), 2, 2, dimnames = list(c("a", "b"), c("f1", "f2")))
  y <- integerise_counts(x, "trunc")
  expect_identical(storage.mode(y), "integer")
  expect_equal(as.vector(y), c(0L, 1L, 3491L, 2L))
  expect_equal(as.vector(integerise_counts(x, "round")), c(0L, 2L, 3492L, 2L))
  expect_error(integerise_counts(x, "none"), "non-integer")
  expect_identical(dimnames(y), dimnames(x))
})

test_that("check_abundance rejects what the pipeline cannot use", {
  m <- matrix(1, 2, 2)
  expect_error(check_abundance(m), "row names")
  dimnames(m) <- list(c("a", "a"), c("x", "y"))
  expect_error(check_abundance(m), "duplicated")
})

test_that("jsd is a distance", {
  p <- c(0.5, 0.5, 0); q <- c(0, 0.5, 0.5)
  expect_equal(jsd(p, p), 0)
  expect_equal(jsd(p, q), jsd(q, p))
  expect_gt(jsd(p, q), 0)
})

test_that("argmax_rows breaks ties reproducibly and reports them", {
  m <- rbind(c(0.5, 0.5, 0, 0), c(0.1, 0.2, 0.7, 0), c(NA, NA, NA, NA))
  a1 <- argmax_rows(m, seed = 1); a2 <- argmax_rows(m, seed = 1)
  expect_identical(a1, a2)
  expect_true(a1$tie[1]); expect_false(a1$tie[2])
  expect_equal(a1$index[2], 3L); expect_true(is.na(a1$index[3]))
})

test_that("union_counts keeps repeated rows with unique names and aligns diversity", {
  A <- matrix(1:12, 4, 3, dimnames = list(c("a", "b", "c", "d"), c("x", "y", "z")))
  div <- c(a = 1, b = 2, c = 3, d = 4)
  u <- union_counts(A, c("a", "b", "a", "d"), div)
  expect_equal(nrow(u$counts), 4)
  expect_equal(rownames(u$counts), c("a", "b", "a.1", "d"))
  expect_equal(unname(u$diversity), c(1, 2, 1, 4))
  expect_equal(u$ids_unique, c("a", "b", "d")); expect_equal(u$keep, c(TRUE, TRUE, FALSE, TRUE))
  expect_error(check_abundance(A[c("a", "a"), ]), "duplicated")
})

test_that("adjusted_rand_index handles edge cases and agrees with the definition", {
  expect_equal(adjusted_rand_index(c(1, 1, 2, 2), c("x", "x", "y", "y")), 1)
  expect_equal(adjusted_rand_index(c(1, 1, 2, 2), c("y", "y", "x", "x")), 1)
  expect_lt(adjusted_rand_index(c(1, 1, 2, 2, 3, 3), c(1, 2, 1, 2, 1, 2)), 0.1)
  expect_true(is.na(adjusted_rand_index("a", "b")))
  expect_true(is.na(adjusted_rand_index(c(NA, "a"), c("b", NA))))
  expect_error(adjusted_rand_index(1:3, 1:2), "same length")
})
