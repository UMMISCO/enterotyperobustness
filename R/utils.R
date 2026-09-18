# Internal helpers: input validation, count integerisation, distances.

#' Validate an abundance matrix
#'
#' @param abundance Numeric matrix, samples in rows, features in columns, with
#'   row and column names.
#' @return The matrix, unchanged, or an error.
#' @keywords internal
check_abundance <- function(abundance) {
  if (!is.matrix(abundance) || !is.numeric(abundance)) {
    stop("`abundance` must be a numeric matrix (samples in rows, features in columns).")
  }
  if (is.null(rownames(abundance)) || is.null(colnames(abundance))) {
    stop("`abundance` needs row names (sample IDs) and column names (feature names).")
  }
  if (anyDuplicated(rownames(abundance))) {
    stop("`abundance` has duplicated sample IDs in its row names.")
  }
  if (any(abundance < 0, na.rm = TRUE) || anyNA(abundance)) {
    stop("`abundance` must be non-negative and complete (no NA).")
  }
  abundance
}

#' Validate a vector of sample IDs against an abundance matrix
#' @keywords internal
check_samples <- function(ids, abundance, what) {
  if (!is.character(ids) || length(ids) == 0) {
    stop(sprintf("`%s` must be a non-empty character vector of sample IDs.", what))
  }
  missing <- setdiff(ids, rownames(abundance))
  if (length(missing)) {
    stop(sprintf("%d %s sample(s) are not rows of `abundance`, e.g. %s",
                 length(missing), what, paste(utils::head(missing, 3), collapse = ", ")))
  }
  if (anyDuplicated(ids)) stop(sprintf("`%s` contains duplicated IDs.", what))
  ids
}

#' Turn a count-like matrix into integer counts
#'
#' The 2022 analyses ran on a rarefied and *scaled* table whose cells were not
#' integers. `DirichletMultinomial` 1.28 accepted it and silently applied
#' `mode(x) <- "integer"`, i.e. truncation towards zero. Newer versions refuse
#' non-integer input, so the same operation is now explicit. `"trunc"` therefore
#' reproduces the published fits; `"round"` does not.
#'
#' @param x Numeric matrix.
#' @param method `"trunc"` (default, reproduces the published analyses),
#'   `"round"`, or `"none"` (error if any cell is not integer-valued).
#' @return An integer matrix with the same dimnames.
#' @export
integerise_counts <- function(x, method = c("trunc", "round", "none")) {
  method <- match.arg(method)
  y <- switch(method,
    trunc = trunc(x),
    round = round(x),
    none  = {
      if (any(x != round(x))) {
        stop("`abundance` has non-integer cells; use `integerise = \"trunc\"` or \"round\".")
      }
      x
    })
  storage.mode(y) <- "integer"
  dimnames(y) <- dimnames(x)
  y
}

#' Jensen-Shannon distance between two probability vectors
#' @keywords internal
jsd <- function(p, q, eps = 1e-12) {
  p <- p / sum(p); q <- q / sum(q)
  m <- (p + q) / 2
  kl <- function(a, b) sum(ifelse(a > 0, a * log((a + eps) / (b + eps)), 0))
  sqrt(0.5 * kl(p, m) + 0.5 * kl(q, m))
}

#' Pairwise Jensen-Shannon distances between the columns of two matrices
#' @return A `ncol(A)` x `ncol(B)` matrix.
#' @keywords internal
jsd_matrix <- function(A, B) {
  out <- matrix(NA_real_, ncol(A), ncol(B), dimnames = list(colnames(A), colnames(B)))
  for (i in seq_len(ncol(A))) for (j in seq_len(ncol(B))) out[i, j] <- jsd(A[, i], B[, j])
  out
}

#' Canonical enterotype names for the four-enterotype scheme
#' @keywords internal
enterotype_names <- function() c("Bact1", "Bact2", "Prevotella", "Ruminococcus")

#' Row-wise argmax with a seeded random tie-break
#'
#' @return A list with `index` (column index per row, NA for all-NA rows) and
#'   `tie` (logical, whether the maximum was tied).
#' @keywords internal
argmax_rows <- function(m, seed = NULL) {
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    set.seed(seed)
  }
  idx <- rep(NA_integer_, nrow(m)); tie <- rep(FALSE, nrow(m))
  for (i in seq_len(nrow(m))) {
    r <- m[i, ]
    if (all(is.na(r))) next
    w <- which(r == max(r, na.rm = TRUE))
    tie[i] <- length(w) > 1
    idx[i] <- if (length(w) > 1) w[sample.int(length(w), 1)] else w
  }
  list(index = idx, tie = tie)
}
