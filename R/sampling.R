#' Draw the complementary samples of one ERP iteration
#'
#' Reproduces the 2022 draws exactly: `set.seed(seed)` then
#' `sample(seq_along(pool), n, replace)`. The positions are returned rather
#' than the IDs so that the draw can be stored compactly and audited.
#'
#' @param pool Character vector of candidate sample IDs, in a fixed order.
#' @param n Number of samples to draw.
#' @param seed Integer seed for this draw (the 2022 runs used the iteration number).
#' @param replace Sample with replacement? The 2022 runs did so only when the
#'   pool was smaller than `n`.
#' @return Integer vector of positions into `pool`.
#' @export
draw_complementary <- function(pool, n, seed, replace = FALSE) {
  if (!replace && n > length(pool)) {
    stop(sprintf("Cannot draw %d of %d samples without replacement.", n, length(pool)))
  }
  set.seed(seed)
  sample.int(length(pool), n, replace = replace)
}

#' Stratified subsample of a cohort, as in ERPM
#'
#' Draws a fraction of a cohort while preserving the class proportions, exactly
#' as the 2022 ERPM scripts did: the total is `round(fraction * N)`, the quota
#' of class *g* is `round(p_g * total)` with `p_g` the class proportion, and
#' classes are drawn one after another in `class_order` (the random stream is
#' shared, so the order matters for reproduction).
#'
#' @param classes Character or factor vector of class labels **named by sample
#'   ID**. Samples without a name are an error.
#' @param fraction Fraction of the cohort to draw, in (0, 1].
#' @param seed Integer seed for this draw. The 2022 runs used the iteration
#'   number for ERPM iterations and 505 for the sub-cohorts of Supplementary
#'   Figure S3.
#' @param class_order Order in which classes are drawn. The default,
#'   `sort(unique(classes))`, reproduces the 2022 order for the MetaCardis
#'   groups (`"1" < "1 RH" < "2a" < "2b" < "3" < ... < "8"`).
#' @param replace Sample with replacement (never in the published analyses).
#' @return Character vector of drawn sample IDs, in draw order.
#' @export
stratified_subsample <- function(classes, fraction, seed,
                                 class_order = sort(unique(as.character(classes))),
                                 replace = FALSE) {
  if (is.null(names(classes)) || any(names(classes) == "")) {
    stop("`classes` must be named by sample ID.")
  }
  if (!is.numeric(fraction) || length(fraction) != 1 || fraction <= 0 || fraction > 1) {
    stop("`fraction` must be a single number in (0, 1].")
  }
  ids <- names(classes)
  classes <- as.character(classes)          # as.character() drops names: keep them aside
  ids_by_class <- split(ids, factor(classes, levels = class_order))
  N <- length(classes)
  total <- round(fraction * N)
  quota <- vapply(ids_by_class, function(ids) round(length(ids) / N * total), numeric(1))
  set.seed(seed)
  drawn <- lapply(seq_along(ids_by_class), function(g) {
    ids <- ids_by_class[[g]]
    if (quota[g] == 0) return(character(0))
    ids[sample.int(length(ids), quota[g], replace = replace)]
  })
  unlist(drawn, use.names = FALSE)
}

#' The 55 cohort fractions used for ERPM in the paper
#'
#' 5 % to 50 % in steps of 1 %, then 55 % to 95 % in steps of 5 %.
#' @export
erpm_default_fractions <- function() {
  c(seq(0.05, 0.50, by = 0.01), seq(0.55, 0.95, by = 0.05))
}
