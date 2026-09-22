#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 11  Design-aware tests for the ERPM figures (Figures 3A, 3B, 4A, 4B), replacing
#     the per-fraction Kruskal-Wallis / Dunn tests of the submitted manuscript:
#     * 3A  robustness between enterotypes at each fraction: Kruskal-Wallis and
#           all-pairs Dunn (independent samples at a given fraction), BH across
#           the whole family (fractions x pairs); curve-level test on per-sample
#           means over the fractions (sample grouped by its consensus at Q = 0.95).
#     * 3B  core vs variable within enterotype: Wilcoxon per fraction, BH within
#           enterotype, direction checked; curve-level test.
#     * 4A  the three pairwise groupings: every sample has one pair value under
#           each grouping, so groupings are compared with a Friedman test (samples
#           as blocks) and Conover post hoc per fraction, BH across fractions;
#           the two pairs within the trade-off grouping (independent samples):
#           Wilcoxon per fraction.
#     * 4B  core vs variable within each trade-off pair: Wilcoxon per fraction, BH.
#
# Inputs : out/cache/erpm_legacy.rds (script 00), mirror figure objects (published core/variable)
# Outputs: out/11_erpm_tests/  fig3A_tests.csv fig3A_curve.csv fig3B_tests.csv fig3B_curve.csv
#          fig4A_tests.csv fig4A_within.csv fig4A_curve.csv fig4B_tests.csv fig4B_curve.csv summary.txt
# ------------------------------------------------------------------------------
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source("inst/scripts/metacardis_data.R"); source("inst/scripts/cache_utils.R")
mirror <- Sys.getenv("ENTEROTYPES_MIRROR")
if (!nzchar(mirror) || !dir.exists(mirror)) stop("Set ENTEROTYPES_MIRROR to the directory holding the 2022 analyses (private, not distributed).")
if (!requireNamespace("PMCMRplus", quietly = TRUE)) stop("PMCMRplus is needed.")
fig <- file.path(mirror, "analyses/9.analyses_mars2026")
out <- "out/11_erpm_tests"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
cache <- load_erpm_cache(mirror)
ent <- enterotype_names(); seed <- 314
freq_int <- function(m) { nd <- rowSums(!is.na(m)); f <- vapply(1:4, function(k) rowSums(m == k, na.rm = TRUE) / nd, numeric(nrow(m)))
  dimnames(f) <- list(rownames(m), ent); f[nd == 0, ] <- NA; attr(f, "n_drawn") <- nd; f }
consensus <- function(f) { am <- argmax_rows(f, seed = seed)
  data.frame(sample = rownames(f), enterotype = ent[am$index], robustness = f[cbind(seq_len(nrow(f)), am$index)], stringsAsFactors = FALSE)[!is.na(am$index), ] }
status_of <- load_published_core_variable(fig)$table[, c("sample", "enterotype", "status")]
status_of <- status_of[!duplicated(status_of[, c("sample", "enterotype")]), ]
pairs_list <- list(comb1 = list(c("Bact1", "Bact2"), c("Prevotella", "Ruminococcus")),
                   comb2 = list(c("Bact2", "Ruminococcus"), c("Bact1", "Prevotella")),
                   comb3 = list(c("Prevotella", "Bact2"), c("Bact1", "Ruminococcus")))
dunn_pairs <- function(x, g) {
  d <- tryCatch(PMCMRplus::kwAllPairsDunnTest(x, factor(g), p.adjust.method = "none")$p.value, error = function(e) NULL)
  if (is.null(d)) return(NULL)
  rows <- which(!is.na(d), arr.ind = TRUE)
  data.frame(g1 = rownames(d)[rows[, 1]], g2 = colnames(d)[rows[, 2]], p = d[rows], stringsAsFactors = FALSE)
}
wilcox_p <- function(a, b) if (length(a) >= 3 && length(b) >= 3) stats::wilcox.test(a, b, exact = FALSE)$p.value else NA_real_

# ---- per fraction: consensus, pair values ---------------------------------------------------
rob <- list(); pairs <- list()
for (key in cache$fractions$file) {
  Q <- cache$fractions$fraction[cache$fractions$file == key]
  f <- freq_int(cache$assignments[[key]]); cs <- consensus(f)
  rob[[key]] <- cbind(fraction = Q, cs)
  pairs[[key]] <- do.call(rbind, lapply(names(pairs_list), function(cb) {
    mp <- merge_pairs(f[!is.na(f[, 1]), ], pairs = pairs_list[[cb]], method = "A", seed = seed)
    cbind(fraction = Q, combination = cb, mp[, c("sample", "pair", "value")]) }))
}
rob <- do.call(rbind, rob); pairs <- do.call(rbind, pairs); rownames(rob) <- rownames(pairs) <- NULL

# ---- 3A: enterotypes at each fraction -----------------------------------------------------------
t3A <- do.call(rbind, lapply(split(rob, rob$fraction), function(x) {
  kw <- stats::kruskal.test(x$robustness, factor(x$enterotype)); d <- dunn_pairs(x$robustness, x$enterotype)
  if (is.null(d)) return(NULL); cbind(fraction = x$fraction[1], kw_p = kw$p.value, d) }))
t3A$p_adj <- stats::p.adjust(t3A$p, "BH")
t3A$rum_lower <- with(t3A, ifelse(g1 == "Ruminococcus" | g2 == "Ruminococcus", TRUE, FALSE))
write.csv(t3A, file.path(out, "fig3A_tests.csv"), row.names = FALSE)
ref <- rob[rob$fraction == max(rob$fraction), c("sample", "enterotype")]
per3A <- aggregate(robustness ~ sample, rob, mean); per3A$enterotype <- ref$enterotype[match(per3A$sample, ref$sample)]
per3A <- per3A[!is.na(per3A$enterotype), ]
c3A <- cbind(kw_p = stats::kruskal.test(per3A$robustness, factor(per3A$enterotype))$p.value, dunn_pairs(per3A$robustness, per3A$enterotype))
c3A$p_adj <- stats::p.adjust(c3A$p, "BH")
c3A$mean_g1 <- tapply(per3A$robustness, per3A$enterotype, mean)[c3A$g1]; c3A$mean_g2 <- tapply(per3A$robustness, per3A$enterotype, mean)[c3A$g2]
write.csv(c3A, file.path(out, "fig3A_curve.csv"), row.names = FALSE)

# ---- 3B: core vs variable within enterotype -----------------------------------------------------
t3B <- do.call(rbind, lapply(split(rob, list(rob$fraction, rob$enterotype), drop = TRUE), function(x) {
  e <- x$enterotype[1]
  rc <- x$robustness[x$sample %in% status_of$sample[status_of$enterotype == e & status_of$status == "core"]]
  rv <- x$robustness[x$sample %in% status_of$sample[status_of$enterotype == e & status_of$status == "variable"]]
  data.frame(fraction = x$fraction[1], enterotype = e, n_core = length(rc), n_variable = length(rv), mean_core = mean(rc), mean_variable = mean(rv),
             difference = mean(rc) - mean(rv), p = wilcox_p(rc, rv)) }))
t3B$p_adj <- stats::ave(t3B$p, t3B$enterotype, FUN = function(v) { a <- v; ok <- !is.na(v); a[ok] <- stats::p.adjust(v[ok], "BH"); a })
t3B$core_higher_adj <- !is.na(t3B$p_adj) & t3B$p_adj < 0.05 & t3B$difference > 0
t3B <- t3B[order(t3B$enterotype, t3B$fraction), ]
write.csv(t3B, file.path(out, "fig3B_tests.csv"), row.names = FALSE)
c3B <- do.call(rbind, lapply(ent, function(e) {
  x <- rob[rob$enterotype == e, ]; per <- tapply(x$robustness, x$sample, mean)
  rc <- per[names(per) %in% status_of$sample[status_of$enterotype == e & status_of$status == "core"]]
  rv <- per[names(per) %in% status_of$sample[status_of$enterotype == e & status_of$status == "variable"]]
  data.frame(enterotype = e, n_core = length(rc), n_variable = length(rv), mean_core = mean(rc), mean_variable = mean(rv), difference = mean(rc) - mean(rv), p = wilcox_p(rc, rv)) }))
c3B$p_adj <- stats::p.adjust(c3B$p, "BH"); write.csv(c3B, file.path(out, "fig3B_curve.csv"), row.names = FALSE)

# ---- 4A: the three groupings (paired by sample) and the two pairs within comb1 -----------------
t4A <- do.call(rbind, lapply(split(pairs, pairs$fraction), function(x) {
  w <- stats::reshape(x[, c("sample", "combination", "value")], idvar = "sample", timevar = "combination", direction = "wide")
  w <- as.matrix(w[stats::complete.cases(w), -1]); colnames(w) <- sub("^value\\.", "", colnames(w))
  ft <- stats::friedman.test(w)
  ph <- tryCatch(PMCMRplus::frdAllPairsConoverTest(y = w, p.adjust.method = "none")$p.value, error = function(e) NULL)
  if (is.null(ph)) return(NULL)
  rows <- which(!is.na(ph), arr.ind = TRUE)
  data.frame(fraction = x$fraction[1], n = nrow(w), friedman_p = ft$p.value, g1 = rownames(ph)[rows[, 1]], g2 = colnames(ph)[rows[, 2]], p = ph[rows],
             mean_g1 = colMeans(w)[rownames(ph)[rows[, 1]]], mean_g2 = colMeans(w)[colnames(ph)[rows[, 2]]], stringsAsFactors = FALSE) }))
t4A$p_adj <- stats::p.adjust(t4A$p, "BH")
write.csv(t4A, file.path(out, "fig4A_tests.csv"), row.names = FALSE)
w4A <- do.call(rbind, lapply(split(pairs[pairs$combination == "comb1", ], pairs$fraction[pairs$combination == "comb1"]), function(x) {
  a <- x$value[x$pair == "Bact1Bact2"]; b <- x$value[x$pair == "PrevotellaRuminococcus"]
  data.frame(fraction = x$fraction[1], mean_Bact1Bact2 = mean(a), mean_PrevotellaRuminococcus = mean(b), p = wilcox_p(a, b)) }))
w4A$p_adj <- stats::p.adjust(w4A$p, "BH"); write.csv(w4A, file.path(out, "fig4A_within.csv"), row.names = FALSE)
per4A <- aggregate(value ~ sample + combination, pairs, mean)
w <- stats::reshape(per4A, idvar = "sample", timevar = "combination", direction = "wide"); w <- as.matrix(w[stats::complete.cases(w), -1]); colnames(w) <- sub("^value\\.", "", colnames(w))
c4A <- data.frame(n = nrow(w), friedman_p = stats::friedman.test(w)$p.value, mean_comb1 = mean(w[, "comb1"]), mean_comb2 = mean(w[, "comb2"]), mean_comb3 = mean(w[, "comb3"]))
ph <- PMCMRplus::frdAllPairsConoverTest(y = w, p.adjust.method = "none")$p.value; rows <- which(!is.na(ph), arr.ind = TRUE)
c4A_pairs <- data.frame(g1 = rownames(ph)[rows[, 1]], g2 = colnames(ph)[rows[, 2]], p = ph[rows]); c4A_pairs$p_adj <- stats::p.adjust(c4A_pairs$p, "BH")
write.csv(cbind(c4A, c4A_pairs), file.path(out, "fig4A_curve.csv"), row.names = FALSE)

# ---- 4B: core vs variable within each trade-off pair -------------------------------------------
p1 <- pairs[pairs$combination == "comb1", ]
status_pair <- function(pr, st) { ents <- pairs_list$comb1[[which(vapply(pairs_list$comb1, function(z) paste(z, collapse = "") == pr, logical(1)))]]
  status_of$sample[status_of$status == st & status_of$enterotype %in% ents] }
t4B <- do.call(rbind, lapply(split(p1, list(p1$fraction, p1$pair), drop = TRUE), function(x) {
  pr <- x$pair[1]; rc <- x$value[x$sample %in% status_pair(pr, "core")]; rv <- x$value[x$sample %in% status_pair(pr, "variable")]
  data.frame(fraction = x$fraction[1], pair = pr, n_core = length(rc), n_variable = length(rv), mean_core = mean(rc), mean_variable = mean(rv),
             difference = mean(rc) - mean(rv), p = wilcox_p(rc, rv)) }))
t4B$p_adj <- stats::ave(t4B$p, t4B$pair, FUN = function(v) { a <- v; ok <- !is.na(v); a[ok] <- stats::p.adjust(v[ok], "BH"); a })
t4B$core_higher_adj <- !is.na(t4B$p_adj) & t4B$p_adj < 0.05 & t4B$difference > 0
t4B <- t4B[order(t4B$pair, t4B$fraction), ]; write.csv(t4B, file.path(out, "fig4B_tests.csv"), row.names = FALSE)
c4B <- do.call(rbind, lapply(unique(p1$pair), function(pr) { x <- p1[p1$pair == pr, ]; per <- tapply(x$value, x$sample, mean)
  rc <- per[names(per) %in% status_pair(pr, "core")]; rv <- per[names(per) %in% status_pair(pr, "variable")]
  data.frame(pair = pr, n_core = length(rc), n_variable = length(rv), mean_core = mean(rc), mean_variable = mean(rv), difference = mean(rc) - mean(rv), p = wilcox_p(rc, rv)) }))
c4B$p_adj <- stats::p.adjust(c4B$p, "BH"); write.csv(c4B, file.path(out, "fig4B_curve.csv"), row.names = FALSE)

# ---- 4B (as the text claims it): among variable samples, and among core samples, does the trade-off grouping
#      beat the two other groupings? Friedman across the three groupings per fraction (samples as blocks), BH; curve level
core_ids <- unique(status_of$sample[status_of$status == "core"]); var_ids <- setdiff(unique(status_of$sample[status_of$status == "variable"]), core_ids)
group_test <- function(x) {
  w <- stats::reshape(x[, c("sample", "combination", "value")], idvar = "sample", timevar = "combination", direction = "wide")
  w <- as.matrix(w[stats::complete.cases(w), -1]); colnames(w) <- sub("^value\\.", "", colnames(w))
  if (nrow(w) < 5) return(NULL)
  ft <- stats::friedman.test(w)
  ph <- tryCatch(PMCMRplus::frdAllPairsConoverTest(y = w, p.adjust.method = "none")$p.value, error = function(e) NULL)
  if (is.null(ph)) return(NULL)
  rows <- which(!is.na(ph), arr.ind = TRUE)
  data.frame(n = nrow(w), friedman_p = ft$p.value, g1 = rownames(ph)[rows[, 1]], g2 = colnames(ph)[rows[, 2]], p = ph[rows],
             mean_g1 = colMeans(w)[rownames(ph)[rows[, 1]]], mean_g2 = colMeans(w)[colnames(ph)[rows[, 2]]], stringsAsFactors = FALSE)
}
g4B <- do.call(rbind, lapply(c("variable", "core"), function(st) {
  ids <- if (st == "variable") var_ids else core_ids
  do.call(rbind, lapply(split(pairs[pairs$sample %in% ids, ], pairs$fraction[pairs$sample %in% ids]), function(x) {
    r <- group_test(x); if (is.null(r)) NULL else cbind(status = st, fraction = x$fraction[1], r) })) }))
g4B$p_adj <- stats::ave(g4B$p, g4B$status, FUN = function(v) stats::p.adjust(v, "BH"))
g4B$comb1_higher_adj <- with(g4B, p_adj < 0.05 & ((g1 == "comb1" & mean_g1 > mean_g2) | (g2 == "comb1" & mean_g2 > mean_g1)))
write.csv(g4B, file.path(out, "fig4B_groupings.csv"), row.names = FALSE)
c4Bg <- do.call(rbind, lapply(c("variable", "core"), function(st) {
  ids <- if (st == "variable") var_ids else core_ids
  per <- aggregate(value ~ sample + combination, pairs[pairs$sample %in% ids, ], mean); r <- group_test(per); cbind(status = st, r) }))
c4Bg$p_adj <- stats::ave(c4Bg$p, c4Bg$status, FUN = function(v) stats::p.adjust(v, "BH"))
write.csv(c4Bg, file.path(out, "fig4B_groupings_curve.csv"), row.names = FALSE)

# ---- summary --------------------------------------------------------------------------------
sink(file.path(out, "summary.txt"))
cat("11  ERPM tests, run", format(Sys.time()), "\n")
cat(sprintf("%d fractions; families BH-adjusted as described in the header\n\n", length(unique(rob$fraction))))
cat("3A  Curve level (per-sample mean over fractions, grouped by consensus at Q = 0.95): Kruskal-Wallis and Dunn pairs\n"); print(c3A, row.names = FALSE, digits = 3)
cat("\n3A  Per fraction: number of fractions (of ", length(unique(t3A$fraction)), ") at which each pair differs (BH over the family):\n", sep = "")
print(aggregate(cbind(sig = p_adj < 0.05) ~ g1 + g2, t3A, sum), row.names = FALSE)
cat("    Ruminococcus lowest: fractions at which every pair involving Ruminococcus is significant and its mean is the lowest:\n")
low <- do.call(rbind, lapply(split(rob, rob$fraction), function(x) data.frame(fraction = x$fraction[1], lowest = names(which.min(tapply(x$robustness, x$enterotype, mean))))))
cat("    Ruminococcus has the lowest mean in", sum(low$lowest == "Ruminococcus"), "of", nrow(low), "fractions\n")
cat("\n3B  Curve level, core vs variable per enterotype:\n"); print(c3B, row.names = FALSE, digits = 3)
cat("    Per fraction: fractions with core > variable (BH within enterotype):\n")
print(aggregate(cbind(core_higher_adj, tested = !is.na(p)) ~ enterotype, t3B, sum), row.names = FALSE)
cat("\n4A  Curve level, the three groupings (Friedman, samples as blocks) and Conover pairs:\n"); print(cbind(c4A, c4A_pairs), row.names = FALSE, digits = 3)
cat("    Per fraction: fractions (of ", length(unique(t4A$fraction)), ") at which comb1 beats each other grouping (BH over the family, higher mean):\n", sep = "")
print(aggregate(cbind(sig_and_higher = p_adj < 0.05 & ((g1 == "comb1" & mean_g1 > mean_g2) | (g2 == "comb1" & mean_g2 > mean_g1))) ~ g1 + g2, t4A, sum), row.names = FALSE)
cat("    Within comb1, Bact1Bact2 vs PrevotellaRuminococcus differ (BH) at", sum(w4A$p_adj < 0.05, na.rm = TRUE), "of", nrow(w4A), "fractions\n")
cat("\n4B  Curve level, core vs variable within each trade-off pair:\n"); print(c4B, row.names = FALSE, digits = 3)
print(aggregate(cbind(core_higher_adj, tested = !is.na(p)) ~ pair, t4B, sum), row.names = FALSE)
cat("\n4B  Among variable and among core samples: the three groupings compared (Friedman, samples as blocks), curve level:\n")
print(c4Bg, row.names = FALSE, digits = 3)
cat("    Per fraction: fractions at which the trade-off grouping (comb1) beats the other grouping (BH within status):\n")
print(aggregate(cbind(comb1_higher_adj, tested = 1) ~ status + g1 + g2, g4B, sum), row.names = FALSE)
cat("\nSession:\n"); print(sessionInfo()); sink()
cat(readLines(file.path(out, "summary.txt"), n = 60), sep = "\n")
