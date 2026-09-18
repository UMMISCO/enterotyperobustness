#' Name the components of a DMM fit as enterotypes
#'
#' Two methods are available.
#'
#' **`"rank"`** is the heuristic of the published analyses
#' (`GivingNameToDMMGroups.R`): the component with the highest mean
#' *Bacteroides* value over its MAP-assigned samples is `Bact2`, the second
#' highest is `Bact1`, the component with the highest mean *Prevotella* is
#' `Prevotella`, and the remaining one is `Ruminococcus`. The means are taken
#' on `counts` exactly as supplied (the 2022 code used the non-integerised
#' table). Unlike the original, a collision (Prevotella falling on a
#' Bacteroides component) does not abort: it is flagged and Prevotella is then
#' taken as the highest *Prevotella* mean among the two remaining components.
#' The size of the Bact1/Bact2 gap and of the Prevotella margin, which the
#' original never recorded, are returned.
#'
#' **`"match"`** ignores marker ranks and matches the fit's component centroids
#' (see [dmm_centroids()]) to a `reference` set of named centroids by
#' Jensen-Shannon distance, solving the assignment with the Hungarian algorithm
#' ([clue::solve_LSAP()]). This is the decoupling of cluster identity from
#' label that the reviewers asked for.
#'
#' In both cases the annotation **checks** of the original are computed and
#' returned as flags, never as errors: whether `Bact2` has the lowest mean
#' `diversity` and `Ruminococcus` the highest, plus any `extra_checks`.
#'
#' @param fit A `DMN` object from [fit_dmm()].
#' @param counts The matrix the fit was made on, samples in rows, **before**
#'   integerisation, with the same row order as the fit.
#' @param marker_features Named character vector giving the column names of
#'   `counts` used for annotation: at least `bacteroides` and `prevotella`
#'   (e.g. `c(bacteroides = "Bacteroides", prevotella = "Prevotella")`).
#' @param diversity Optional numeric vector named by sample ID (gene richness
#'   in the paper), used for the two fundamental checks.
#' @param method `"rank"` or `"match"`.
#' @param reference For `method = "match"`: a features x enterotypes matrix of
#'   centroids with column names `Bact1`, `Bact2`, `Prevotella`, `Ruminococcus`,
#'   typically [reference_centroids()] of a fit on the full cohort.
#' @param extra_checks Optional named list of character vectors of feature
#'   names. For each, the check is that `Ruminococcus` has the highest mean of
#'   the row sums of those features. The 2022 code used
#'   `list(FS = c("Faecalibacterium", "Subdoligranulum"), FSR = c(..., "Ruminococcus"),
#'   FSRE = c(..., "Eubacterium"))`.
#' @return A list of class `enterotype_annotation`:
#' * `map`: named integer vector, enterotype name -> component index;
#' * `component_names`: the enterotype name of each component, in component order;
#' * `flags`: one-row data frame of logical checks (`TRUE` = anomaly);
#' * `gap`: one-row data frame with the Bact1/Bact2 gap (absolute and relative
#'   to Bact2) and the Prevotella margin, or the matching cost and margin for
#'   `"match"`;
#' * `group_means`: components x markers matrix of the means used.
#' @export
annotate_clusters <- function(fit, counts, marker_features, diversity = NULL,
                              method = c("rank", "match"), reference = NULL,
                              extra_checks = NULL) {
  method <- match.arg(method)
  counts <- check_abundance(counts)
  comp <- dmm_assign(fit)
  K <- ncol(dmm_posterior(fit))
  if (K != 4) stop("The four-enterotype annotation scheme needs K = 4 (got ", K, ").")
  if (length(comp) != nrow(counts)) {
    stop("`counts` does not have one row per sample of `fit`.")
  }
  need <- c("bacteroides", "prevotella")
  if (!all(need %in% names(marker_features))) {
    stop("`marker_features` must name at least `bacteroides` and `prevotella`.")
  }
  absent <- setdiff(marker_features, colnames(counts))
  if (length(absent)) stop("Marker feature(s) not in `counts`: ", paste(absent, collapse = ", "))

  members <- split(seq_len(nrow(counts)), factor(comp, levels = seq_len(K)))
  comp_mean <- function(feature) {
    vapply(members, function(rows) if (length(rows)) mean(counts[rows, feature]) else NA_real_,
           numeric(1))
  }
  bact <- comp_mean(marker_features[["bacteroides"]])
  prev <- comp_mean(marker_features[["prevotella"]])
  empty <- vapply(members, length, integer(1)) == 0

  if (method == "rank") {
    rank_of <- function(x) { x[is.na(x)] <- -Inf; x }
    b <- rank_of(bact); p <- rank_of(prev)
    Bact2 <- which.max(b)
    Bact1 <- setdiff(order(b, decreasing = TRUE), Bact2)[1]
    Prevotella <- which.max(p)
    collision <- Prevotella %in% c(Bact1, Bact2)
    if (collision) {
      rest <- setdiff(seq_len(K), c(Bact1, Bact2))
      Prevotella <- rest[which.max(p[rest])]
    }
    Ruminococcus <- setdiff(seq_len(K), c(Bact1, Bact2, Prevotella))
    map <- c(Bact1 = unname(Bact1), Bact2 = unname(Bact2),
             Prevotella = unname(Prevotella), Ruminococcus = unname(Ruminococcus))
    gap <- data.frame(
      bact_gap = bact[Bact2] - bact[Bact1],
      bact_gap_relative = (bact[Bact2] - bact[Bact1]) / bact[Bact2],
      prevotella_margin = prev[Prevotella] - max(prev[-Prevotella], na.rm = TRUE),
      prevotella_collision = collision)
  } else {
    if (is.null(reference)) stop("`method = \"match\"` needs `reference` centroids.")
    if (!all(enterotype_names() %in% colnames(reference))) {
      stop("`reference` must have columns ", paste(enterotype_names(), collapse = ", "))
    }
    cent <- dmm_centroids(fit)
    common <- intersect(rownames(cent), rownames(reference))
    if (length(common) < 2) stop("`reference` and the fit share fewer than two features.")
    cost <- jsd_matrix(cent[common, , drop = FALSE], reference[common, enterotype_names(), drop = FALSE])
    sol <- clue::solve_LSAP(cost)              # component i -> reference column sol[i]
    map <- stats::setNames(integer(K), enterotype_names())
    map[colnames(cost)[as.integer(sol)]] <- seq_len(K)
    assigned <- cost[cbind(seq_len(K), as.integer(sol))]
    second <- vapply(seq_len(K), function(i) min(cost[i, -as.integer(sol)[i]]), numeric(1))
    gap <- data.frame(
      matching_cost = sum(assigned),
      matching_margin_min = min(second - assigned),
      bact_gap = bact[map[["Bact2"]]] - bact[map[["Bact1"]]],
      bact_gap_relative = (bact[map[["Bact2"]]] - bact[map[["Bact1"]]]) / bact[map[["Bact2"]]])
  }

  component_names <- character(K)
  component_names[map] <- names(map)

  flags <- data.frame(empty_component = any(empty))
  if (!is.null(diversity)) {
    if (is.null(names(diversity))) stop("`diversity` must be named by sample ID.")
    div <- vapply(members, function(rows) {
      if (!length(rows)) return(NA_real_)
      mean(diversity[rownames(counts)[rows]])
    }, numeric(1))
    flags$bact2_not_lowest_diversity <- !isTRUE(div[map[["Bact2"]]] == min(div, na.rm = TRUE))
    flags$rum_not_highest_diversity  <- !isTRUE(div[map[["Ruminococcus"]]] == max(div, na.rm = TRUE))
  }
  if (!is.null(extra_checks)) {
    for (nm in names(extra_checks)) {
      feats <- extra_checks[[nm]]
      miss <- setdiff(feats, colnames(counts))
      if (length(miss)) stop("extra_checks$", nm, ": feature(s) not in `counts`: ", paste(miss, collapse = ", "))
      s <- vapply(members, function(rows) if (length(rows)) mean(rowSums(counts[rows, feats, drop = FALSE])) else NA_real_, numeric(1))
      flags[[paste0("rum_not_highest_", nm)]] <- !isTRUE(s[map[["Ruminococcus"]]] == max(s, na.rm = TRUE))
    }
  }

  structure(list(
    method = method,
    map = map,
    component_names = component_names,
    flags = flags,
    gap = gap,
    group_means = cbind(bacteroides = bact, prevotella = prev)
  ), class = "enterotype_annotation")
}

#' Reference centroids for `annotate_clusters(method = "match")`
#'
#' Fits a DMM on a reference table (typically the full cohort), names its
#' components with the rank heuristic, and returns the named centroids.
#'
#' @inheritParams annotate_clusters
#' @inheritParams fit_dmm
#' @param fit Optional: an existing `DMN` fit on `counts`; fitted if `NULL`.
#' @return A features x 4 matrix of centroids with enterotype column names, with
#'   the annotation attached as attribute `"annotation"`.
#' @export
reference_centroids <- function(counts, marker_features, diversity = NULL, fit = NULL,
                                K = 4, seed = 314, integerise = "trunc") {
  if (is.null(fit)) fit <- fit_dmm(counts, K = K, seed = seed, integerise = integerise)
  ann <- annotate_clusters(fit, counts, marker_features, diversity, method = "rank")
  cent <- dmm_centroids(fit)
  colnames(cent) <- ann$component_names
  cent <- cent[, enterotype_names(), drop = FALSE]
  attr(cent, "annotation") <- ann
  cent
}

#' @export
print.enterotype_annotation <- function(x, ...) {
  cat("Enterotype annotation (", x$method, ")\n", sep = "")
  cat("  components:", paste(sprintf("%d=%s", seq_along(x$component_names), x$component_names), collapse = "  "), "\n")
  fl <- unlist(x$flags)
  cat("  flags     :", if (any(fl)) paste(names(fl)[fl], collapse = ", ") else "none", "\n")
  cat("  gap       :", paste(sprintf("%s=%.3g", names(x$gap), unlist(x$gap)), collapse = "  "), "\n")
  invisible(x)
}
