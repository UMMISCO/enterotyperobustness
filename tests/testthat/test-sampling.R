test_that("draw_complementary reproduces set.seed(i); sample()", {
  pool <- sprintf("p%02d", 1:20)
  set.seed(7); ref <- sample.int(20, 5)
  expect_identical(draw_complementary(pool, 5, seed = 7), ref)
  expect_error(draw_complementary(pool, 30, seed = 1), "without replacement")
  expect_length(draw_complementary(pool, 30, seed = 1, replace = TRUE), 30)
})

test_that("stratified_subsample preserves class proportions with the 2022 rounding", {
  classes <- stats::setNames(rep(c("1", "1 RH", "2a", "8"), c(50, 10, 30, 110)),
                             sprintf("s%03d", 1:200))
  ids <- stratified_subsample(classes, 0.3, seed = 505)
  expect_length(ids, 60)
  expect_equal(as.vector(table(classes[ids])[c("1", "1 RH", "2a", "8")]), c(15, 3, 9, 33))
  expect_identical(ids, stratified_subsample(classes, 0.3, seed = 505))
  expect_false(identical(ids, stratified_subsample(classes, 0.3, seed = 506)))
  expect_error(stratified_subsample(unname(classes), 0.3, 1), "named")
})

test_that("default class order is the string sort used in 2022", {
  classes <- stats::setNames(c("8", "2b", "1 RH", "1", "2a", "3"), letters[1:6])
  expect_identical(sort(unique(classes)), c("1", "1 RH", "2a", "2b", "3", "8"))
})

test_that("the paper's fraction grid has 55 values", {
  f <- erpm_default_fractions()
  expect_length(f, 55)
  expect_equal(range(f), c(0.05, 0.95))
})
