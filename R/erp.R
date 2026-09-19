#' Enterotype Robustness Process (ERP): perturb cohort composition
#'
#' Holds an **anchor set** of samples fixed and, at each iteration, clusters it
#' together with a random draw from a **complementary pool**, fits a DMM,
#' names the components as enterotypes, and records every sample's label and
#' posterior. Aggregated over iterations this gives each anchor sample's
#' assignment frequency under that clinical context (see
#' [assignment_frequency()]).
#'
#' With the defaults (`iteration_seeds = 1:n_iter`, `dmm_seed = 314`,
#' `integerise = "trunc"`, `annotate = "rank"`) and the anchor and pool IDs in
#' the same order as the 2022 scripts, the function reproduces the published
#' experiments. Each ERP condition of the paper was one call: the healthy
#' anchor (142 IDs) with a pool of healthy, obese, T2D, statin- or
#' metformin-stratified samples and `n_complementary = 86`.
#'
#' @param abundance Numeric matrix, samples in rows, features in columns, with
#'   dimnames. Any cohort; the table is never modified.
#' @param anchors Character vector of anchor sample IDs (rows of `abundance`).
#' @param complementary Character vector of pool sample IDs, disjoint from
#'   `anchors`. **Order matters for reproduction**: draws are positions into this
#'   vector.
#' @param n_complementary Number of pool samples drawn per iteration. The 2022
#'   runs used `round(0.3 * cohort_size)` where `cohort_size` is the size of the
#'   cohort the anchors were sampled from (86 for 285 healthy, 135 for 451
#'   T2D-metformin, 125 for 418 obese). Give either this or `p_complementary`
#'   with `cohort_size`.
#' @param p_complementary,cohort_size Alternative to `n_complementary`:
#'   `n_complementary = round(p_complementary * cohort_size)`.
#' @param n_iter Number of iterations (500 in the paper).
#' @param marker_features,diversity,extra_checks Passed to [annotate_clusters()].
#' @param K Number of DMM components (4).
#' @param seed Tie-breaking seed recorded in `params` and used by default by
#'   [robustness()], [labels_by_condition()] and [merge_pairs()] on this result.
#' @param iteration_seeds Integer vector of length `n_iter`: the seed of each
#'   iteration's draw. The default `1:n_iter` reproduces the 2022 runs.
#' @param dmm_seed Seed passed to `dmn()` in every iteration (314 in the paper).
#' @param replace `"auto"` (with replacement only if the pool is smaller than
#'   `n_complementary`, as the 2022 runs did), `"never"` or `"always"`.
#' @param annotate `"rank"` (published heuristic) or `"match"` (centroid
#'   matching against `reference`); see [annotate_clusters()].
#' @param reference Centroids for `annotate = "match"`; see [reference_centroids()].
#' @param integerise See [integerise_counts()]; `"trunc"` reproduces the paper.
#' @param n_cores Number of forked workers ([parallel::mclapply()]). Results do
#'   not depend on it.
#' @param keep_fits Keep the 500 `DMN` objects in the result (about 15 MB per
#'   ERP condition). Needed for centroid matching after the fact.
#' @param verbose Print progress.
#' @return An object of class `erp_result`, a list with:
#' * `assignments`: character matrix, samples x iterations, the enterotype
#'   assigned in that iteration or `NA` when the sample was not drawn;
#' * `posteriors`: list of samples x 4 posterior matrices, one per iteration,
#'   columns always `Bact1, Bact2, Prevotella, Ruminococcus`;
#' * `label_map`: integer matrix, iterations x 4, the DMM component of each name;
#' * `index`: list of the drawn positions into `complementary`;
#' * `flags`: data frame, one row per iteration, the annotation checks;
#' * `gap`: data frame, one row per iteration, the Bact1/Bact2 gap and margins;
#' * `fits`: the `DMN` objects, or `NULL`;
#' * `samples`: list with `anchors` and `complementary`;
#' * `params`: every argument that determines the result.
#' @examples
#' X <- simulate_enterotype_counts(n = 60)
#' markers <- c(bacteroides = "Bacteroides", prevotella = "Prevotella")
#' anchors <- rownames(X)[1:30]; pool <- rownames(X)[31:60]
#' res <- erp(X, anchors, pool, n_complementary = 10, n_iter = 3,
#'            marker_features = markers, diversity = attr(X, "diversity"),
#'            keep_fits = FALSE, verbose = FALSE)
#' res
#' head(assignment_frequency(res))
#' head(robustness(res))
#' @seealso [erpm()], [assignment_frequency()], [robustness()], [core_variable()]
#' @export
erp <- function(abundance, anchors, complementary,
                n_complementary = NULL, p_complementary = NULL, cohort_size = NULL,
                n_iter = 500, marker_features, diversity = NULL, K = 4,
                seed = 314, iteration_seeds = seq_len(n_iter), dmm_seed = 314,
                replace = c("auto", "never", "always"),
                annotate = c("rank", "match"), reference = NULL, extra_checks = NULL,
                integerise = c("trunc", "round", "none"),
                n_cores = 1, keep_fits = TRUE, verbose = TRUE) {
  abundance <- check_abundance(abundance)
  anchors <- check_samples(anchors, abundance, "anchors")
  complementary <- check_samples(complementary, abundance, "complementary")
  overlap <- intersect(anchors, complementary)
  if (length(overlap)) {
    stop(length(overlap), " sample(s) are both anchors and complementary, e.g. ", overlap[1])
  }
  if (is.null(n_complementary)) {
    if (is.null(p_complementary) || is.null(cohort_size)) {
      stop("Give `n_complementary`, or both `p_complementary` and `cohort_size`.")
    }
    n_complementary <- round(p_complementary * cohort_size)
  }
  if (length(iteration_seeds) != n_iter) stop("`iteration_seeds` must have length `n_iter`.")
  check_diversity(diversity, c(anchors, complementary))
  replace <- match.arg(replace)
  annotate <- match.arg(annotate)
  integerise <- match.arg(integerise)
  use_replace <- switch(replace,
    auto = length(complementary) < n_complementary,
    never = FALSE, always = TRUE)
  if (use_replace && verbose) {
    message("Pool (", length(complementary), ") smaller than n_complementary (",
            n_complementary, "): sampling with replacement.")
  }

  index <- lapply(seq_len(n_iter), function(i) {
    draw_complementary(complementary, n_complementary, iteration_seeds[i], use_replace)
  })
  id_list <- lapply(index, function(ix) c(anchors, complementary[ix]))

  runs <- run_dmm_iterations(abundance, id_list, K, dmm_seed, integerise,
                             marker_features, diversity, annotate, reference,
                             extra_checks, n_cores, keep_fits, verbose)
  out <- assemble_result(runs, all_samples = c(anchors, complementary), keep_fits)
  out$index <- index
  out$samples <- list(anchors = anchors, complementary = complementary)
  out$params <- list(procedure = "ERP", n_complementary = n_complementary,
                     replace = use_replace, n_iter = n_iter, K = K, seed = seed,
                     iteration_seeds = iteration_seeds, dmm_seed = dmm_seed,
                     annotate = annotate, integerise = integerise,
                     marker_features = marker_features,
                     has_diversity = !is.null(diversity), extra_checks = names(extra_checks))
  class(out) <- c("erp_result", "dmm_iterations")
  out
}

#' Run DMM fit + annotation over a list of sample-ID sets (internal)
#'
#' Shared by [erp()] and [erpm()]. Iteration `i` clusters `abundance[id_list[[i]], ]`
#' in that row order. Deterministic for any `n_cores` because the draws are
#' pre-computed and the DMM seed is fixed.
#' @keywords internal
run_dmm_iterations <- function(abundance, id_list, K, dmm_seed, integerise,
                               marker_features, diversity, annotate, reference,
                               extra_checks, n_cores, keep_fits, verbose) {
  n_iter <- length(id_list)
  one <- function(i) {
    ids <- id_list[[i]]
    u <- union_counts(abundance, ids, diversity)
    X <- u$counts
    fit <- fit_dmm(X, K = K, seed = dmm_seed, integerise = integerise)
    ann <- annotate_clusters(fit, X, marker_features, u$diversity,
                             method = annotate, reference = reference,
                             extra_checks = extra_checks)
    post <- dmm_posterior(fit)
    colnames(post) <- ann$component_names
    post <- post[, enterotype_names(), drop = FALSE]
    labels <- ann$component_names[dmm_assign(fit)]
    keep <- !duplicated(ids)                     # one entry per sample and iteration
    post <- post[keep, , drop = FALSE]
    rownames(post) <- ids[keep]
    list(ids = ids[keep], labels = labels[keep], posterior = post,
         map = ann$map, flags = ann$flags, gap = ann$gap,
         fit = if (keep_fits) fit else NULL)
  }
  if (verbose) message("Fitting ", n_iter, " DMM models (K = ", K, ") on ", n_cores, " core(s)...")
  t0 <- Sys.time()
  runs <- if (n_cores > 1) {
    # One fork per iteration (no prescheduling). Note that dmn() itself is not
    # run-to-run deterministic: DirichletMultinomial defines calloc as R_alloc
    # (un-zeroed) and kmeans() reads that memory before writing it, so a few
    # percent of fits on multimodal draws reach a different EM optimum from one
    # run to the next, sequential or parallel. A fresh fork per fit reduces the
    # effect; it cannot remove it. Aggregates (assignment frequencies, consensus
    # labels) are unaffected; see inst/scripts/determinism_full.R.
    parallel::mclapply(seq_len(n_iter), one, mc.cores = n_cores,
                       mc.preschedule = FALSE, mc.set.seed = FALSE)
  } else {
    lapply(seq_len(n_iter), one)
  }
  failed <- vapply(runs, inherits, logical(1), "try-error")
  if (any(failed)) {
    first <- which(failed)[1]
    stop("Iteration(s) failed: ", paste(which(failed), collapse = ", "),
         ". First error: ", trimws(as.character(runs[[first]])))
  }
  if (verbose) message("Done in ", format(round(difftime(Sys.time(), t0, units = "mins"), 1)), ".")
  runs
}

#' Assemble per-iteration runs into the result object (internal)
#' @keywords internal
assemble_result <- function(runs, all_samples, keep_fits) {
  n_iter <- length(runs)
  assignments <- matrix(NA_character_, length(all_samples), n_iter,
                        dimnames = list(all_samples, paste0("iter", seq_len(n_iter))))
  for (i in seq_len(n_iter)) assignments[runs[[i]]$ids, i] <- runs[[i]]$labels
  label_map <- do.call(rbind, lapply(runs, function(r) r$map))
  rownames(label_map) <- colnames(assignments)
  bind_rows <- function(field) {          # tolerate empty per-iteration data frames
    d <- do.call(rbind, lapply(runs, function(r) r[[field]]))
    if (is.null(d) || ncol(d) == 0) d <- data.frame(row.names = seq_len(n_iter))
    cbind(iteration = seq_len(n_iter), d)
  }
  list(
    assignments = assignments,
    posteriors = lapply(runs, function(r) r$posterior),
    label_map = label_map,
    flags = bind_rows("flags"),
    gap = bind_rows("gap"),
    fits = if (keep_fits) lapply(runs, function(r) r$fit) else NULL
  )
}

#' @export
print.erp_result <- function(x, ...) {
  cat("ERP result: ", length(x$samples$anchors), " anchors + ", x$params$n_complementary,
      " of ", length(x$samples$complementary), " complementary samples per iteration, ",
      x$params$n_iter, " iterations, K = ", x$params$K,
      ", annotation = ", x$params$annotate, "\n", sep = "")
  fl <- x$flags[, setdiff(names(x$flags), "iteration"), drop = FALSE]
  cat("Flagged iterations:", paste(sprintf("%s=%d", names(fl), colSums(fl)), collapse = "  "), "\n")
  if (!is.null(x$fits)) cat("DMM fits kept:", length(x$fits), "\n")
  invisible(x)
}
