#' Assignment frequency of every sample
#'
#' The fraction of the iterations **in which a sample was drawn** that assigned
#' it to each enterotype. This is the quantity the paper calls "posterior
#' probability" and plots in every figure; it is an empirical frequency, not a
#' DMM posterior (for that see [mean_posterior()]).
#'
#' @param x An `erp_result`, an `erpm_fraction`, or a samples x iterations
#'   character matrix with `NA` where a sample was not drawn.
#' @return A samples x 4 numeric matrix (rows sum to 1; `NA` rows for samples
#'   never drawn), with attribute `n_drawn`.
#' @examples
#' a <- rbind(s1 = c("Bact1", "Bact1", NA), s2 = c("Bact2", "Bact1", "Bact2"))
#' assignment_frequency(a)
#' assignment_entropy(assignment_frequency(a))
#' @export
assignment_frequency <- function(x) {
  if (inherits(x, "erpm_result")) {
    stop("`x` is a whole ERPM result; call assignment_frequency() on one fraction, e.g. x$results[[1]].")
  }
  a <- if (is.matrix(x)) x else x$assignments
  if (!is.matrix(a) || !is.character(a)) stop("`x` must be an erp_result, an erpm_fraction or a character matrix.")
  n_drawn <- rowSums(!is.na(a))
  freq <- vapply(enterotype_names(), function(e) rowSums(a == e, na.rm = TRUE) / n_drawn,
                 numeric(nrow(a)))
  freq <- matrix(freq, nrow(a), 4, dimnames = list(rownames(a), enterotype_names()))
  freq[n_drawn == 0, ] <- NA_real_
  attr(freq, "n_drawn") <- n_drawn
  freq
}

#' Normalised Shannon entropy of assignment frequencies
#'
#' `-sum(p log p) / log(K)`, so that 0 means one enterotype in every iteration
#' and 1 means all four equally often. Base 2 as in the paper (the base cancels
#' after normalisation).
#'
#' @param freq A samples x K frequency matrix, e.g. from [assignment_frequency()].
#'   If it has the four enterotype columns, only those are used, whatever else
#'   is bound to it.
#' @param base Logarithm base.
#' @return A named numeric vector.
#' @export
assignment_entropy <- function(freq, base = 2) {
  if (!is.numeric(base) || length(base) != 1 || !is.finite(base) || base <= 0 || base == 1) stop("`base` must be a positive number other than 1.")
  if (all(enterotype_names() %in% colnames(freq))) freq <- freq[, enterotype_names(), drop = FALSE]
  if (any(freq < 0, na.rm = TRUE) || any(freq > 1, na.rm = TRUE)) stop("`freq` must hold frequencies in [0, 1].")
  rs <- rowSums(freq); bad <- !is.na(rs) & abs(rs - 1) > 1e-6
  if (any(bad)) stop(sum(bad), " row(s) of `freq` do not sum to 1 (first: ", rownames(freq)[which(bad)[1]], ", sum ", signif(rs[which(bad)[1]], 4), ").")
  K <- ncol(freq)
  plogp <- function(p) ifelse(p > 0, p * log(p, base), 0)
  h <- -rowSums(plogp(freq)) / log(K, base)
  stats::setNames(as.numeric(h), rownames(freq))
}

#' Consensus enterotype and robustness of every sample
#'
#' The consensus enterotype is the most frequent one across iterations;
#' robustness is that maximum frequency (the paper's "enterotype robustness").
#' Ties are broken at random with `seed`, and counted (the `W` column of the
#' paper's Supplementary Tables S2 and S3 is `sum(tie)`).
#'
#' @inheritParams assignment_frequency
#' @param seed Seed for tie-breaking. Default: the `seed` recorded in `x$params`
#'   when `x` is a result object, 314 otherwise.
#' @return A data frame with `sample`, `enterotype`, `robustness`, `n_drawn`, `tie`.
#' @export
robustness <- function(x, seed = result_seed(x)) {
  freq <- if (is.matrix(x) && is.numeric(x)) x else assignment_frequency(x)
  am <- argmax_rows(freq, seed = seed)
  n_drawn <- attr(freq, "n_drawn")
  data.frame(
    sample = rownames(freq),
    enterotype = ifelse(is.na(am$index), NA_character_, colnames(freq)[am$index]),
    robustness = ifelse(is.na(am$index), NA_real_, freq[cbind(seq_len(nrow(freq)), am$index)]),
    n_drawn = if (is.null(n_drawn)) NA_integer_ else n_drawn,
    tie = am$tie,
    stringsAsFactors = FALSE, row.names = NULL)
}

#' Mean DMM posterior of every sample across iterations
#'
#' The average of the per-iteration posterior component probabilities (columns
#' aligned by enterotype name). This is the 2022 `ConfidenceSample` object; no
#' paper figure uses it, but it is the quantity the Methods' wording "maximum
#' posterior probability" literally describes.
#'
#' @param x An `erp_result` or `erpm_fraction`.
#' @return A samples x 4 matrix, plus attribute `sd` with the per-cell sample
#'   standard deviation (denominator n - 1, as [stats::sd()]) across the
#'   iterations in which the sample was drawn; `NA` when drawn fewer than twice.
#' @export
mean_posterior <- function(x) {
  samples <- rownames(x$assignments)
  acc <- matrix(0, length(samples), 4, dimnames = list(samples, enterotype_names()))
  acc2 <- acc; n <- numeric(length(samples)); names(n) <- samples
  for (p in x$posteriors) {
    ids <- rownames(p)
    acc[ids, ] <- acc[ids, ] + p
    acc2[ids, ] <- acc2[ids, ] + p^2
    n[ids] <- n[ids] + 1
  }
  m <- acc / n
  v <- (acc2 - n * m^2) / (n - 1)
  v[n < 2, ] <- NA_real_
  m[n == 0, ] <- NA_real_
  attr(m, "sd") <- sqrt(pmax(v, 0))
  m
}

#' Core and variable samples across ERP conditions
#'
#' The published rule (`rule = "unanimous"`, `code_map.md` D-1): for each
#' enterotype *e*, a sample is **core** for *e* if it was assigned *e* in every
#' condition, and **variable** for *e* if it was assigned *e* in at least one
#' but not all conditions. A sample can therefore be variable for two
#' enterotypes. Two sensitivity variants are provided for the reviewers:
#' `"k_minus_1"` relaxes "every" to "all but one", and `"sd"` is the
#' criterion described in the submitted Methods (core if the standard deviation
#' of each of the four assignment frequencies across conditions is below
#' `threshold`), which the published figures did not actually use.
#'
#' @param labels Character matrix, samples x conditions, the consensus
#'   enterotype of each anchor sample in each ERP condition (see
#'   [labels_by_condition()]). Samples with a missing label in any condition
#'   are dropped with a message: the rule is defined on complete rows only.
#' @param rule `"unanimous"`, `"k_minus_1"` or `"sd"`.
#' @param frequencies For `rule = "sd"`: a list, one samples x 4 frequency
#'   matrix per condition (same row order as `labels`).
#' @param threshold For `rule = "sd"`.
#' @return A data frame with `sample`, `enterotype`, `status` (`core` or
#'   `variable`) and, for `"sd"`, `max_sd`.
#' @examples
#' lab <- rbind(a = c("Bact1", "Bact1", "Bact1"), b = c("Bact1", "Bact2", "Bact1"))
#' core_variable(lab)                  # a: core Bact1; b: variable for Bact1 and Bact2
#' core_variable(lab, "k_minus_1")     # b becomes core Bact1
#' @export
core_variable <- function(labels, rule = c("unanimous", "k_minus_1", "sd"),
                          frequencies = NULL, threshold = 0.2) {
  rule <- match.arg(rule)
  if (!is.matrix(labels)) stop("`labels` must be a samples x conditions character matrix.")
  incomplete <- !stats::complete.cases(labels)
  if (any(incomplete)) {
    message(sum(incomplete), " sample(s) without a label in every condition dropped from core_variable().")
    labels <- labels[!incomplete, , drop = FALSE]
  }
  k <- ncol(labels)
  if (rule %in% c("unanimous", "k_minus_1")) {
    need <- if (rule == "unanimous") k else k - 1
    out <- lapply(enterotype_names(), function(e) {
      n_e <- rowSums(labels == e, na.rm = TRUE)
      core <- rownames(labels)[n_e >= need & n_e > 0]
      variable <- rownames(labels)[n_e > 0 & n_e < need]
      rbind(
        if (length(core)) data.frame(sample = core, enterotype = e, status = "core"),
        if (length(variable)) data.frame(sample = variable, enterotype = e, status = "variable"))
    })
    out <- do.call(rbind, out)
    if (is.null(out)) out <- data.frame(sample = character(0), enterotype = character(0), status = character(0))
    rownames(out) <- NULL
    return(out)
  }
  if (is.null(frequencies)) stop("`rule = \"sd\"` needs `frequencies`.")
  arr <- simplify2array(lapply(frequencies, function(f) f[rownames(labels), enterotype_names(), drop = FALSE]))
  if (length(dim(arr)) != 3) stop("`frequencies` must be a list of samples x 4 matrices, one per condition.")
  sds <- apply(arr, c(1, 2), stats::sd)                    # samples x 4
  means <- apply(arr, c(1, 2), mean)
  if (anyNA(means)) stop("`frequencies` contain NA for ", sum(rowSums(is.na(means)) > 0), " sample(s).")
  max_sd <- apply(sds, 1, max)
  data.frame(sample = rownames(labels),
             enterotype = enterotype_names()[apply(means, 1, which.max)],
             status = ifelse(max_sd < threshold, "core", "variable"),
             max_sd = max_sd, stringsAsFactors = FALSE, row.names = NULL)
}

#' Consensus labels of the anchors across several ERP results
#'
#' @param results Named list of `erp_result` objects sharing the same anchors
#'   (one per condition).
#' @param samples Which samples to keep (default: the anchors of the first result).
#' @param seed Tie-breaking seed passed to [robustness()]; `NULL` uses each
#'   result's own recorded seed.
#' @return A samples x conditions character matrix.
#' @export
labels_by_condition <- function(results, samples = results[[1]]$samples$anchors, seed = NULL) {
  if (is.null(names(results))) stop("`results` must be a named list (one name per condition).")
  m <- sapply(results, function(r) {
    rb <- if (is.null(seed)) robustness(r) else robustness(r, seed = seed)
    stats::setNames(rb$enterotype, rb$sample)[samples]
  })
  m <- matrix(m, length(samples), length(results), dimnames = list(samples, names(results)))
  m
}

#' Collapse the four enterotypes into pairs
#'
#' Reproduces Figure 4 of the paper. **Method A** (the one retained in the
#' paper) keeps each sample's four-way consensus enterotype and assigns the
#' sample to the pair that contains it; the pair value is the sum of the two
#' frequencies. **Method B** re-assigns each sample to the pair with the larger
#' summed frequency.
#'
#' @param freq A samples x 4 frequency matrix from [assignment_frequency()].
#' @param pairs A list of two character vectors partitioning the four names.
#'   Default: the Bact1/Bact2 vs Prevotella/Ruminococcus trade-off axes.
#' @param method `"A"` or `"B"`.
#' @param seed Tie-breaking seed.
#' @return A data frame with `sample`, `pair`, `value`, `tie`.
#' @export
merge_pairs <- function(freq, pairs = list(c("Bact1", "Bact2"), c("Prevotella", "Ruminococcus")),
                        method = c("A", "B"), seed = 314) {
  method <- match.arg(method)
  if (length(pairs) != 2 || anyDuplicated(unlist(pairs)) || !setequal(unlist(pairs), enterotype_names())) {
    stop("`pairs` must be two disjoint vectors partitioning ", paste(enterotype_names(), collapse = ", "))
  }
  pair_names <- vapply(pairs, paste, character(1), collapse = "")
  sums <- sapply(pairs, function(p) rowSums(freq[, p, drop = FALSE]))
  sums <- matrix(sums, nrow(freq), 2, dimnames = list(rownames(freq), pair_names))
  if (method == "A") {
    am <- argmax_rows(freq, seed = seed)
    which_pair <- vapply(am$index, function(i) {
      if (is.na(i)) NA_integer_ else which(vapply(pairs, function(p) colnames(freq)[i] %in% p, logical(1)))
    }, integer(1))
    tie <- am$tie
  } else {
    am <- argmax_rows(sums, seed = seed)
    which_pair <- am$index; tie <- am$tie
  }
  data.frame(sample = rownames(freq),
             pair = ifelse(is.na(which_pair), NA_character_, pair_names[which_pair]),
             value = ifelse(is.na(which_pair), NA_real_, sums[cbind(seq_len(nrow(sums)), which_pair)]),
             tie = tie, stringsAsFactors = FALSE, row.names = NULL)
}

#' Tie-breaking seed recorded in a result object (internal)
#' @keywords internal
result_seed <- function(x) {
  s <- if (is.list(x) && !is.null(x$params$seed)) x$params$seed else 314
  as.integer(s)
}
