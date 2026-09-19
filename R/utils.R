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

#' Seed the RNG for the calling function only (internal)
#'
#' Sets the seed for the rest of the caller's frame and restores the caller's
#' `.Random.seed` on exit (or removes it if none existed), so that no exported
#' function leaves the user's random stream in a different state.
#' @keywords internal
local_seed <- function(seed, envir = parent.frame()) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
  restore <- function() {
    if (had) assign(".Random.seed", old, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  }
  do.call(on.exit, list(bquote(.(restore)()), add = TRUE), envir = envir)
  set.seed(seed)
  invisible(NULL)
}

#' Row-wise argmax with a seeded random tie-break
#'
#' @return A list with `index` (column index per row, NA for all-NA rows) and
#'   `tie` (logical, whether the maximum was tied).
#' @keywords internal
argmax_rows <- function(m, seed = NULL) {
  if (!is.null(seed)) local_seed(seed)
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

#' Adjusted Rand index between two labellings
#'
#' Hubert and Arabie's adjusted Rand index, computed from the contingency
#' table; 1 for identical partitions (up to renaming), about 0 for independent ones.
#'
#' @param a,b Vectors of labels of the same length (NA pairs are dropped).
#' @export
adjusted_rand_index <- function(a, b) {
  if (length(a) != length(b)) stop("`a` and `b` must have the same length.")
  keep <- !is.na(a) & !is.na(b); a <- a[keep]; b <- b[keep]
  if (length(a) < 2) return(NA_real_)
  tab <- table(a, b); n <- sum(tab)
  comb2 <- function(x) x * (x - 1) / 2
  sum_ij <- sum(comb2(tab)); sum_a <- sum(comb2(rowSums(tab))); sum_b <- sum(comb2(colSums(tab)))
  expected <- sum_a * sum_b / comb2(n); maxi <- (sum_a + sum_b) / 2
  if (maxi == expected) return(1)
  (sum_ij - expected) / (maxi - expected)
}

#' Count matrix of one iteration's samples, tolerant to repeated draws
#'
#' Draws with replacement repeat a sample within an iteration. The repeated rows
#' were part of the 2022 fits and must stay, but row names have to be unique for
#' lookups, so they are made unique here and the per-sample `diversity` is
#' aligned to them. Used by [erp()], [erpm()] and the legacy readers.
#'
#' @param abundance Samples x features matrix.
#' @param ids Sample IDs of the iteration, possibly repeated, in fit order.
#' @param diversity Optional named numeric vector.
#' @return A list with `counts`, `diversity` (or `NULL`), `keep` (logical, first
#'   occurrence of each ID) and `ids_unique`.
#' @export
union_counts <- function(abundance, ids, diversity = NULL) {
  X <- abundance[ids, , drop = FALSE]
  rownames(X) <- make.unique(ids)
  div <- if (is.null(diversity)) NULL else stats::setNames(diversity[ids], rownames(X))
  keep <- !duplicated(ids)
  list(counts = X, diversity = div, keep = keep, ids_unique = ids[keep])
}

#' Check that a diversity vector covers a set of sample IDs (internal)
#' @keywords internal
check_diversity <- function(diversity, ids) {
  if (is.null(diversity)) return(invisible(NULL))
  if (!is.numeric(diversity) || is.null(names(diversity))) stop("`diversity` must be a numeric vector named by sample ID.")
  miss <- setdiff(ids, names(diversity))
  if (length(miss)) {
    stop(length(miss), " sample(s) have no `diversity` value, e.g. ", miss[1],
         ". A missing value would silently flag every iteration containing the sample.")
  }
  if (anyNA(diversity[ids])) stop("`diversity` is NA for ", sum(is.na(diversity[ids])), " sample(s).")
  invisible(NULL)
}
