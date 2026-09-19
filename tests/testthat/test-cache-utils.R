# The cache helpers live in inst/scripts (they are script infrastructure, not
# package API); they are tested here so that a regression cannot silently
# reuse a stale intermediate result.
f <- system.file("scripts", "cache_utils.R", package = "enterotyperobustness")
if (!nzchar(f)) f <- test_path("..", "..", "inst", "scripts", "cache_utils.R")
skip_if_not(file.exists(f), "cache_utils.R not found")
source(f, local = TRUE)

test_that("a cache is reused only when parameters and the input matrix match", {
  withr::local_envvar(CACHE_STRICT = "", CACHE_MIGRATE = "")
  A <- matrix(1:12, 3, 4, dimnames = list(paste0("s", 1:3), paste0("f", 1:4)))
  m <- cache_manifest(list(K = 4, algorithm = "v1"), inputs = list(A = A))
  file <- withr::local_tempfile(fileext = ".rds")
  cache_write("result", file, m)
  expect_identical(suppressMessages(cache_read(file, m)), "result")
  expect_null(suppressMessages(cache_read(file, cache_manifest(list(K = 5, algorithm = "v1"), inputs = list(A = A)))))
  B <- A; B[1, 1] <- B[1, 1] + 1; B[1, 2] <- B[1, 2] - 1                 # same dims, same total, different cells
  expect_null(suppressMessages(cache_read(file, cache_manifest(list(K = 4, algorithm = "v1"), inputs = list(A = B)))))
  m_dmm <- m; m_dmm$dmm_version <- "0.0.0"
  expect_null(suppressMessages(cache_read(file, m_dmm)))                  # a different DirichletMultinomial: never reused
  m_commit <- m; m_commit$git_commit <- "deadbee"
  expect_identical(suppressMessages(cache_read(file, m_commit)), "result")   # commit differences: message only ...
  withr::local_envvar(CACHE_STRICT = "1")
  expect_null(suppressMessages(cache_read(file, m_commit)))              # ... unless the run is strict
})

test_that("a cache without manifest is ignored unless migration is requested explicitly", {
  withr::local_envvar(CACHE_STRICT = "", CACHE_MIGRATE = "")
  file <- withr::local_tempfile(fileext = ".rds"); saveRDS(list(a = 1), file)
  m <- cache_manifest(list(algorithm = "v1"), inputs = list(1:3))
  expect_null(suppressMessages(cache_read(file, m)))
  withr::local_envvar(CACHE_MIGRATE = "1")
  expect_identical(suppressMessages(cache_read(file, m)), list(a = 1))
  withr::local_envvar(CACHE_MIGRATE = "")
  expect_identical(suppressMessages(cache_read(file, m)), list(a = 1))   # the manifest was written during migration
  expect_true(!is.null(readRDS(file)$manifest))
})
