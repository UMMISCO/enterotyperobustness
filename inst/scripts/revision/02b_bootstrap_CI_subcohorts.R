#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 02b BCa bootstrap intervals and design-aware tests for Supplementary Figure S5
#     (ERPM on sub-cohorts of the MetaCardis dataset), the counterpart of script
#     02 for Figure 3A. For every sub-cohort, fraction and enterotype: mean
#     robustness with a BCa bootstrap interval over samples. Tests: at each
#     fraction and enterotype, Kruskal-Wallis across sub-cohorts and all-pairs
#     Dunn tests, BH-adjusted across the whole family (fractions x pairs) within
#     enterotype, plus a curve-level test on per-sample means over the fractions.
#
# Inputs : out/cache/erpm_subcohorts.rds (script 00b)
# Outputs: out/02b_subcohorts/  figS3_ci.csv, figS3_tests.csv, figS3_curve_level.csv,
#          FigureS3_draft.pdf, summary.txt  (file names predate the renumbering:
#          the sub-cohorts are Supplementary Figure S5, drawn by reproduce_figures.R)
# ------------------------------------------------------------------------------
suppressMessages({ pkgload::load_all(".", quiet = TRUE); library(ggplot2) })
source("inst/scripts/cache_utils.R")
mirror <- Sys.getenv("ENTEROTYPES_MIRROR")
if (!nzchar(mirror) || !dir.exists(mirror)) stop("Set ENTEROTYPES_MIRROR to the directory holding the 2022 analyses (private, not distributed).")
if (!requireNamespace("PMCMRplus", quietly = TRUE)) stop("PMCMRplus is needed.")
out <- "out/02b_subcohorts"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
cache <- load_subcohort_cache(mirror)
ent <- enterotype_names(); R_samples <- 2000; seed <- 314
freq_int <- function(m) { nd <- rowSums(!is.na(m)); f <- vapply(1:4, function(k) rowSums(m == k, na.rm = TRUE) / nd, numeric(nrow(m)))
  dimnames(f) <- list(rownames(m), ent); f[nd == 0, ] <- NA; attr(f, "n_drawn") <- nd; f }
consensus <- function(f) { am <- argmax_rows(f, seed = seed)
  data.frame(sample = rownames(f), enterotype = ent[am$index], robustness = f[cbind(seq_len(nrow(f)), am$index)], stringsAsFactors = FALSE)[!is.na(am$index), ] }
runs <- cache$runs
runs$label <- sprintf("%d%% of the cohort (n = %d)", round(100 * runs$n_samples / 2022), runs$n_samples)   # ConfidenceN=2000 is the 90 % sub-cohort (1 820 samples)
rob <- do.call(rbind, lapply(seq_len(nrow(runs)), function(i) {
  key <- file.path(runs$subcohort[i], runs$file[i]); cs <- consensus(freq_int(cache$assignments[[key]]))
  cbind(subcohort = runs$label[i], fraction = runs$fraction[i], cs) }))
ci <- do.call(rbind, lapply(split(rob, list(rob$subcohort, rob$fraction, rob$enterotype), drop = TRUE), function(x) {
  b <- bootstrap_ci(x$robustness, R = R_samples, type = "bca", seed = seed)
  data.frame(subcohort = x$subcohort[1], fraction = x$fraction[1], group = x$enterotype[1], n = length(x$robustness),
             mean = b[["mean"]], lower = b[["lower"]], upper = b[["upper"]]) }))

##subset rob and ci to common sample fractions (exclude sub-cohort fractions <50% in in 30% set)
tfractions.table <- table(rob$subcohort, rob$fraction) ##table of subcohort x fractions
tfractions.core <- names(colSums(tfractions.table>0)[colSums(tfractions.table>0)==nrow(tfractions.table)]) ##get fractions common to all subcohorts
rob <- rob[rob$fraction %in% as.numeric(tfractions.core),] ##subset rob to tfractions.core
ci <- ci[ci$fraction %in% as.numeric(tfractions.core),] ##subset ci to tfractions.core
##continue with the analyses
ci <- ci[order(ci$group, ci$subcohort, ci$fraction), ]; rownames(ci) <- NULL
write.csv(ci, file.path(out, "figS3_ci.csv"), row.names = FALSE)
# tests across sub-cohorts at each fraction, per enterotype
tests <- do.call(rbind, lapply(split(rob, list(rob$fraction, rob$enterotype), drop = TRUE), function(x) {
  if (length(unique(x$subcohort)) < 2) return(NULL)
  g <- factor(x$subcohort); kw <- stats::kruskal.test(x$robustness, g)
  d <- tryCatch(PMCMRplus::kwAllPairsDunnTest(x$robustness, g, p.adjust.method = "none")$p.value, error = function(e) NULL)
  if (is.null(d)) return(NULL)
  rows <- which(!is.na(d), arr.ind = TRUE)
  data.frame(enterotype = x$enterotype[1], fraction = x$fraction[1], kw_p = kw$p.value,
             sub1 = rownames(d)[rows[, 1]], sub2 = colnames(d)[rows[, 2]], p = d[rows], stringsAsFactors = FALSE) }))
tests$p_adj <- stats::ave(tests$p, tests$enterotype, FUN = function(v) stats::p.adjust(v, "BH"))
write.csv(tests, file.path(out, "figS3_tests.csv"), row.names = FALSE)
curve <- do.call(rbind, lapply(split(rob, rob$enterotype), function(x) {
  per <- aggregate(robustness ~ sample + subcohort, x, mean)
  kw <- stats::kruskal.test(per$robustness, factor(per$subcohort))
  data.frame(enterotype = x$enterotype[1], n_subcohorts = length(unique(per$subcohort)), kw_p = kw$p.value) }))
write.csv(curve, file.path(out, "figS3_curve_level.csv"), row.names = FALSE)
p <- ggplot(ci, aes(fraction, mean, colour = subcohort, fill = subcohort)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, colour = NA) + geom_line() + geom_point(size = 0.7) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey40") + facet_wrap(~group, nrow = 1) +
  coord_cartesian(ylim = c(min(ci$lower), 1)) + labs(x = "fraction of the sub-cohort", y = "mean robustness", colour = "sub-cohort", fill = "sub-cohort") + theme_bw()
ggsave(file.path(out, "FigureS3_draft.pdf"), p, width = 13, height = 4)
sink(file.path(out, "summary.txt"))
cat("02b  Sub-cohort bootstrap intervals and tests, run", format(Sys.time()), "\n\n")
print(unique(runs[, c("subcohort", "n_samples")]), row.names = FALSE)
cat("\nFirst fraction at which the BCa lower bound exceeds 0.95, per sub-cohort and enterotype:\n")
q <- aggregate(fraction ~ subcohort + group, ci[ci$lower > 0.95, ], min)
print(stats::reshape(q, idvar = "subcohort", timevar = "group", direction = "wide"), row.names = FALSE)
cat("\nCurve-level Kruskal-Wallis across sub-cohorts on per-sample mean robustness:\n"); print(curve, row.names = FALSE, digits = 3)
cat("\nPairs of sub-cohorts differing (BH within enterotype, p_adj < 0.05), count of fractions:\n")
print(aggregate(cbind(sig = p_adj < 0.05, n = 1) ~ enterotype + sub1 + sub2, tests, sum), row.names = FALSE)
cat("\nSession:\n"); print(sessionInfo()); sink()
cat(readLines(file.path(out, "summary.txt"), n = 40), sep = "\n")
