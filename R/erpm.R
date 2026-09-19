#' Enterotype Robustness Process, size-perturbation variant (ERPM)
#'
#' At each of `n_iter` iterations draws a stratified fraction of the cohort
#' (see [stratified_subsample()]), fits a DMM, names the components and records
#' every sample's label and posterior, for every fraction in `fractions`. With
#' the defaults it reproduces the 2022 runs (seeds `1:n_iter` per draw, DMM seed
#' 314, truncation of non-integer cells, rank annotation, the 55 fractions of
#' the paper).
#'
#' The Supplementary Figure S3 sub-cohorts are obtained by drawing the
#' sub-cohort first (`stratified_subsample(classes, M, seed = 505)`) and calling
#' `erpm()` on `abundance[sub, ]` with `classes[sub]`.
#'
#' @inheritParams erp
#' @param classes Character or factor vector of class labels **named by sample
#'   ID** (clinical groups in the paper). Samples absent from `abundance` are
#'   dropped with a message; samples of `abundance` without a class are ignored.
#' @param fractions Fractions of the cohort to draw; default
#'   [erpm_default_fractions()].
#' @param class_order Passed to [stratified_subsample()].
#' @param keep_fits Keep the `DMN` objects in memory (`FALSE` by default: 500
#'   fits at 95 % of 2 022 samples are about 160 MB per fraction).
#' @param fits_dir If not `NULL`, the fits of each fraction are saved there as
#'   `fits_Q<fraction>.rds` instead of being kept in memory.
#' @return An object of class `erpm_result`: a list with `fractions`, `results`
#'   (one `erp_result`-like object per fraction, whose `index` holds the drawn
#'   sample IDs and whose `samples$cohort` is the whole cohort) and `params`.
#' @examples
#' X <- simulate_enterotype_counts(n = 60)
#' classes <- stats::setNames(c("A", "B")[attr(X, "group") %% 2 + 1], rownames(X))
#' res <- erpm(X, classes, fractions = c(0.5, 0.8), n_iter = 2,
#'             marker_features = c(bacteroides = "Bacteroides", prevotella = "Prevotella"),
#'             verbose = FALSE)
#' res
#' head(robustness(res$results[["Q0.80"]]))
#' @seealso [erp()], [robustness()], [merge_pairs()]
#' @export
erpm <- function(abundance, classes, fractions = erpm_default_fractions(), n_iter = 500,
                 marker_features, diversity = NULL, K = 4,
                 seed = 314, iteration_seeds = seq_len(n_iter), dmm_seed = 314,
                 class_order = NULL,
                 annotate = c("rank", "match"), reference = NULL, extra_checks = NULL,
                 integerise = c("trunc", "round", "none"),
                 n_cores = 1, keep_fits = FALSE, fits_dir = NULL, verbose = TRUE) {
  abundance <- check_abundance(abundance)
  if (is.null(names(classes)) || any(names(classes) == "")) stop("`classes` must be named by sample ID.")
  if (anyDuplicated(names(classes))) stop("`classes` has duplicated sample IDs, e.g. ", names(classes)[duplicated(names(classes))][1])
  absent <- setdiff(names(classes), rownames(abundance))
  if (length(absent)) {
    message(length(absent), " classified sample(s) are not in `abundance` and are dropped.")
    classes <- classes[setdiff(names(classes), absent)]
  }
  classes <- stats::setNames(as.character(classes), names(classes))
  if (is.null(class_order)) class_order <- sort(unique(classes), method = "radix")
  n_iter <- check_count(n_iter, "n_iter"); K <- check_count(K, "K", 2); n_cores <- check_count(n_cores, "n_cores")
  if (!is.numeric(fractions) || !length(fractions) || any(!is.finite(fractions)) || any(fractions <= 0 | fractions > 1)) {
    stop("`fractions` must be numbers in (0, 1].")
  }
  keys <- sprintf("Q%.2f", fractions)
  if (anyDuplicated(keys)) {
    stop("Fractions are named to two decimals (", paste(keys[duplicated(keys)], collapse = ", "),
         " would collide); use a grid with at most two decimals.")
  }
  if (length(iteration_seeds) != n_iter) stop("`iteration_seeds` must have length `n_iter`.")
  iteration_seeds <- check_seeds(iteration_seeds, "iteration_seeds", n_iter); dmm_seed <- check_seeds(dmm_seed, "dmm_seed", 1)
  check_diversity(diversity, names(classes))
  annotate <- match.arg(annotate)
  integerise <- match.arg(integerise)
  if (!is.null(fits_dir)) dir.create(fits_dir, showWarnings = FALSE, recursive = TRUE)

  results <- vector("list", length(fractions))
  names(results) <- keys
  for (f in seq_along(fractions)) {
    Q <- fractions[f]
    if (verbose) message("ERPM fraction ", Q, " (", f, "/", length(fractions), ")")
    id_list <- lapply(seq_len(n_iter), function(i) {
      stratified_subsample(classes, Q, seed = iteration_seeds[i], class_order = class_order)
    })
    runs <- run_dmm_iterations(abundance, id_list, K, dmm_seed, integerise,
                               marker_features, diversity, annotate, reference,
                               extra_checks, n_cores, keep_fits || !is.null(fits_dir), verbose)
    res <- assemble_result(runs, all_samples = names(classes), keep_fits = TRUE)
    if (!is.null(fits_dir)) {
      saveRDS(res$fits, file.path(fits_dir, sprintf("fits_Q%.2f.rds", Q)))
      res$fits <- NULL
    } else if (!keep_fits) {
      res$fits <- NULL
    }
    res$index <- id_list
    res$fraction <- Q
    res$samples <- list(cohort = names(classes))
    res$params <- list(procedure = "ERPM", fraction = Q, n_iter = n_iter, K = K, seed = seed,
                       iteration_seeds = iteration_seeds, dmm_seed = dmm_seed, replace = FALSE,
                       annotate = annotate, integerise = integerise, class_order = class_order)
    class(res) <- c("erpm_fraction", "dmm_iterations")
    results[[f]] <- res
  }
  structure(list(
    fractions = fractions,
    results = results,
    params = list(procedure = "ERPM", n_iter = n_iter, K = K, seed = seed,
                  iteration_seeds = iteration_seeds, dmm_seed = dmm_seed,
                  class_order = class_order, annotate = annotate,
                  integerise = integerise, marker_features = marker_features,
                  has_diversity = !is.null(diversity), extra_checks = names(extra_checks),
                  cohort_size = length(classes),
                  class_sizes = table(classes)[class_order])
  ), class = "erpm_result")
}

#' @export
print.erpm_fraction <- function(x, ...) {
  cat("ERPM fraction ", x$params$fraction, ": ", x$params$n_iter, " iterations of ",
      if (length(x$index)) length(x$index[[1]]) else NA, " of ", length(x$samples$cohort),
      " samples, K = ", x$params$K, ", annotation = ", x$params$annotate, "\n", sep = "")
  fl <- x$flags[, setdiff(names(x$flags), "iteration"), drop = FALSE]
  if (ncol(fl)) cat("Flagged iterations:", paste(sprintf("%s=%d", names(fl), colSums(fl, na.rm = TRUE)), collapse = "  "), "\n")
  invisible(x)
}

#' @export
print.erpm_result <- function(x, ...) {
  cat("ERPM result: cohort of ", x$params$cohort_size, " samples, ",
      length(x$fractions), " fractions (", min(x$fractions), " to ", max(x$fractions), "), ",
      x$params$n_iter, " iterations each, K = ", x$params$K,
      ", annotation = ", x$params$annotate, "\n", sep = "")
  invisible(x)
}
