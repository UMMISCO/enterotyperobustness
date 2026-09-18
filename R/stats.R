# Statistical tests that respect the design: the same anchor samples appear in
# every ERP condition, so conditions are repeated measures, not independent
# groups. The `legacy_*` functions reproduce the tests of the submitted
# manuscript so that the two can be compared.

#' Friedman test with post hoc comparisons on repeated measures
#'
#' For each stratum, keeps the samples observed in every condition, runs
#' [stats::friedman.test()] (samples as blocks) and an all-pairs post hoc test
#' from \pkg{PMCMRplus} (Conover by default). This replaces the unpaired
#' Kruskal-Wallis + Dunn tests of the submitted manuscript, which treated the
#' ERP conditions as independent groups although every anchor sample appears in
#' each of them.
#'
#' @param data Long data frame, one row per sample and condition.
#' @param sample Column with the sample ID.
#' @param condition Column with the condition (ERP experiment).
#' @param value Column with the value to compare (entropy).
#' @param stratum Optional column name; a **per-sample** grouping (e.g. the
#'   sample's enterotype in the reference condition). It must be constant within
#'   a sample; a per-condition label would break the pairing.
#' @param post_hoc `"conover"` or `"nemenyi"`.
#' @param p_adjust Method for [stats::p.adjust()] applied to the pairwise
#'   p-values within each stratum, and to the omnibus p-values across strata.
#' @return A list with `omnibus` (one row per stratum: `n`, `statistic`, `df`,
#'   `p`, `p_adj`) and `pairwise` (one row per stratum and pair: `cond1`,
#'   `cond2`, `p`, `p_adj`).
#' @export
paired_entropy_test <- function(data, sample = "sample", condition = "condition",
                                value = "entropy", stratum = NULL,
                                post_hoc = c("conover", "nemenyi"), p_adjust = "BH") {
  post_hoc <- match.arg(post_hoc)
  if (!requireNamespace("PMCMRplus", quietly = TRUE)) stop("Package PMCMRplus is needed.")
  strata <- if (is.null(stratum)) list(all = data) else split(data, data[[stratum]])
  omni <- list(); pairs <- list()
  for (s in names(strata)) {
    d <- strata[[s]]
    wide <- stats::reshape(d[, c(sample, condition, value)], idvar = sample, timevar = condition,
                           direction = "wide")
    rownames(wide) <- wide[[sample]]; wide[[sample]] <- NULL
    colnames(wide) <- sub(paste0("^", value, "\\."), "", colnames(wide))
    wide <- as.matrix(wide[stats::complete.cases(wide), , drop = FALSE])
    if (nrow(wide) < 3 || ncol(wide) < 2) next
    ft <- stats::friedman.test(wide)
    omni[[s]] <- data.frame(stratum = s, n = nrow(wide), statistic = unname(ft$statistic),
                            df = unname(ft$parameter), p = ft$p.value)
    ph <- if (post_hoc == "conover") {
      PMCMRplus::frdAllPairsConoverTest(y = wide, p.adjust.method = "none")
    } else PMCMRplus::frdAllPairsNemenyiTest(y = wide)
    pm <- ph$p.value
    rows <- which(!is.na(pm), arr.ind = TRUE)
    pw <- data.frame(stratum = s, cond1 = rownames(pm)[rows[, 1]], cond2 = colnames(pm)[rows[, 2]],
                     p = pm[rows], stringsAsFactors = FALSE)
    pw$p_adj <- stats::p.adjust(pw$p, method = p_adjust)
    pairs[[s]] <- pw
  }
  omni <- do.call(rbind, omni); rownames(omni) <- NULL
  omni$p_adj <- stats::p.adjust(omni$p, method = p_adjust)
  pairs <- do.call(rbind, pairs); rownames(pairs) <- NULL
  list(omnibus = omni, pairwise = pairs, post_hoc = post_hoc, p_adjust = p_adjust)
}

#' Cochran's Q test: does the frequency of one enterotype differ across conditions?
#'
#' Cochran's Q is the Friedman test on a binary outcome; it tests whether the
#' proportion of anchor samples assigned to `enterotype` changes across
#' conditions, with samples as blocks. Post hoc: pairwise McNemar tests.
#'
#' @param labels Character matrix, samples x conditions, e.g. from
#'   [labels_by_condition()].
#' @param enterotype The label to test.
#' @param p_adjust Method for the pairwise p-values.
#' @return A list with `test` (statistic, df, p, prevalence per condition) and
#'   `pairwise` (McNemar p and adjusted p per pair).
#' @export
cochran_q_test <- function(labels, enterotype, p_adjust = "BH") {
  ind <- (labels == enterotype) * 1
  ft <- stats::friedman.test(ind)
  conds <- colnames(labels)
  cmb <- utils::combn(conds, 2)
  pw <- data.frame(cond1 = cmb[1, ], cond2 = cmb[2, ], stringsAsFactors = FALSE)
  pw$p <- apply(cmb, 2, function(cc) {
    tab <- table(factor(ind[, cc[1]], levels = 0:1), factor(ind[, cc[2]], levels = 0:1))
    b <- tab[1, 2]; c <- tab[2, 1]
    if (b + c == 0) return(1)
    stats::mcnemar.test(tab, correct = TRUE)$p.value
  })
  pw$p_adj <- stats::p.adjust(pw$p, method = p_adjust)
  list(test = data.frame(enterotype = enterotype, n = nrow(labels), statistic = unname(ft$statistic),
                         df = unname(ft$parameter), p = ft$p.value),
       prevalence = colMeans(ind), pairwise = pw)
}

#' Permutation test of condition x enterotype tables with anchors as blocks
#'
#' The submitted manuscript compared enterotype distributions across ERP
#' conditions with a Chi-square test on the condition x label table, counting
#' every anchor once per condition as if the conditions were independent
#' samples. Here the Chi-square statistic is kept but its null distribution is
#' obtained by permuting each sample's labels **across its own conditions**
#' (exchangeability within sample under the null of no condition effect), which
#' respects the shared-anchor structure. Cell-level post hoc p-values use the
#' adjusted standardised residuals under the same permutation scheme.
#'
#' @param labels Character matrix, samples x conditions.
#' @param n_perm Number of permutations.
#' @param seed Seed.
#' @param p_adjust Method for the cell-level p-values.
#' @return A list with `table` (condition x label counts), `statistic`, `p`,
#'   `n_perm`, and `cells` (condition, enterotype, observed, expected,
#'   residual, p, p_adj).
#' @export
anchor_permutation_test <- function(labels, n_perm = 9999, seed = 314, p_adjust = "BH") {
  lev <- sort(unique(as.vector(labels)))
  conds <- colnames(labels); n <- nrow(labels); k <- ncol(labels)
  counts <- function(m) {
    t(vapply(seq_len(k), function(j) tabulate(match(m[, j], lev), nbins = length(lev)), integer(length(lev))))
  }
  chisq_stat <- function(tab) {
    E <- outer(rowSums(tab), colSums(tab)) / sum(tab)
    list(stat = sum((tab - E)^2 / E, na.rm = TRUE), E = E)
  }
  std_resid <- function(tab, E) {
    N <- sum(tab); r <- rowSums(tab) / N; c <- colSums(tab) / N
    (tab - E) / sqrt(E * outer(1 - r, 1 - c))
  }
  obs <- counts(labels); dimnames(obs) <- list(conds, lev)
  o <- chisq_stat(obs); res_obs <- std_resid(obs, o$E)
  set.seed(seed)
  stat_perm <- numeric(n_perm); res_ge <- matrix(0, k, length(lev))
  m <- labels
  for (b in seq_len(n_perm)) {
    for (i in seq_len(n)) m[i, ] <- labels[i, sample.int(k)]
    tab <- counts(m); pb <- chisq_stat(tab)
    stat_perm[b] <- pb$stat
    res_ge <- res_ge + (abs(std_resid(tab, pb$E)) >= abs(res_obs))
  }
  p <- (1 + sum(stat_perm >= o$stat)) / (1 + n_perm)
  cells <- data.frame(condition = rep(conds, times = length(lev)), enterotype = rep(lev, each = k),
                      observed = as.vector(obs), expected = as.vector(o$E), residual = as.vector(res_obs),
                      p = as.vector((1 + res_ge) / (1 + n_perm)), stringsAsFactors = FALSE)
  cells$p_adj <- stats::p.adjust(cells$p, method = p_adjust)
  list(table = obs, statistic = o$stat, p = p, n_perm = n_perm, cells = cells)
}

#' Bootstrap confidence interval of a mean
#'
#' Replaces the normal-approximation interval of the submitted manuscript for
#' robustness scores, which are bounded in \[0, 1\] and pile up near 1.
#'
#' @param x Numeric vector.
#' @param R Bootstrap replicates.
#' @param type `"bca"` or `"percentile"`.
#' @param conf Confidence level.
#' @param seed Seed.
#' @return Named numeric: `mean`, `lower`, `upper`, `n`.
#' @export
bootstrap_ci <- function(x, R = 2000, type = c("bca", "percentile"), conf = 0.95, seed = 314) {
  type <- match.arg(type)
  x <- x[!is.na(x)]; n <- length(x)
  if (n < 2 || length(unique(x)) == 1) return(c(mean = mean(x), lower = mean(x), upper = mean(x), n = n))
  if (!requireNamespace("boot", quietly = TRUE)) stop("Package boot is needed.")
  set.seed(seed)
  b <- boot::boot(x, function(d, i) mean(d[i]), R = R)
  ci <- tryCatch(boot::boot.ci(b, conf = conf, type = type), error = function(e) NULL)
  if (is.null(ci)) { type <- "percentile"; ci <- boot::boot.ci(b, conf = conf, type = "perc") }
  int <- if (type == "bca") ci$bca[4:5] else ci$percent[4:5]
  c(mean = mean(x), lower = int[1], upper = int[2], n = n)
}

#' Normal-approximation interval, as the Methods described it
#'
#' `mean +/- 1.96 * sd / sqrt(n)`. The submitted figure scripts divided by `n`
#' instead of `sqrt(n)`; this function exists so the two can be shown side by
#' side in the response to the reviewers.
#' @param x Numeric vector.
#' @param divisor `"sqrt_n"` (correct) or `"n"` (what the 2026 figure scripts did).
#' @export
normal_ci <- function(x, divisor = c("sqrt_n", "n")) {
  divisor <- match.arg(divisor)
  x <- x[!is.na(x)]; n <- length(x); m <- mean(x); s <- stats::sd(x)
  h <- 1.96 * s / if (divisor == "sqrt_n") sqrt(n) else n
  c(mean = m, lower = m - h, upper = m + h, n = n)
}

#' The Chi-square test of the submitted manuscript, with residual post hoc
#'
#' Counts every anchor once per condition and runs [stats::chisq.test()] on the
#' condition x label table, then converts each adjusted standardised residual
#' to a p-value (`chisq(1)` of its square) with BH adjustment across cells,
#' which is what `chisq.posthoc.test(method = "BH")` does. Kept to regenerate
#' the published numbers; see [anchor_permutation_test()] for the design-aware
#' version.
#'
#' @inheritParams anchor_permutation_test
#' @export
legacy_chisq <- function(labels, p_adjust = "BH") {
  tab <- table(condition = rep(colnames(labels), each = nrow(labels)), enterotype = as.vector(labels))
  ct <- suppressWarnings(stats::chisq.test(tab))
  z <- ct$stdres
  cells <- data.frame(condition = rep(rownames(z), times = ncol(z)), enterotype = rep(colnames(z), each = nrow(z)),
                      observed = as.vector(tab), residual = as.vector(z),
                      p = as.vector(stats::pchisq(z^2, 1, lower.tail = FALSE)), stringsAsFactors = FALSE)
  cells$p_adj <- stats::p.adjust(cells$p, method = p_adjust)
  list(table = tab, statistic = unname(ct$statistic), df = unname(ct$parameter), p = ct$p.value, cells = cells)
}

#' The Kruskal-Wallis + Dunn test of the submitted manuscript
#'
#' Unpaired, within each per-condition enterotype stratum, `dunn.test` with
#' `method = "hochberg"` (labelled "FDR" in the manuscript). Requires
#' \pkg{dunn.test}. Kept for comparison only.
#'
#' @inheritParams paired_entropy_test
#' @param method Passed to `dunn.test::dunn.test`.
#' @export
legacy_dunn <- function(data, condition = "condition", value = "entropy", stratum = NULL,
                        method = "hochberg") {
  if (!requireNamespace("dunn.test", quietly = TRUE)) stop("Package dunn.test is needed.")
  strata <- if (is.null(stratum)) list(all = data) else split(data, data[[stratum]])
  out <- lapply(names(strata), function(s) {
    d <- strata[[s]]
    kw <- stats::kruskal.test(d[[value]] ~ factor(d[[condition]]))
    invisible(utils::capture.output(suppressMessages(
      res <- dunn.test::dunn.test(d[[value]], d[[condition]], method = method, kw = FALSE, table = FALSE, list = FALSE))))
    data.frame(stratum = s, comparison = res$comparisons, p = res$P, p_adj = res$P.adjusted,
               kw_p = kw$p.value, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out); rownames(out) <- NULL
  out
}
