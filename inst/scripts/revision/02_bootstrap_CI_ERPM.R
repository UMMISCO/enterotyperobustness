#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 02  Bootstrap confidence intervals for ERPM robustness (Reviewer 1, major point 6)
#
# The submitted figures drew mean +/- 1.96 * sd / n (sic: n, not sqrt(n)) on a
# score bounded in [0, 1]. Here, per cohort fraction and enterotype (Fig. 3A),
# per core/variable status (3B), per enterotype pair (4A) and pair x status (4B):
#   * BCa bootstrap over SAMPLES (primary; n = samples of the group at that fraction),
#   * bootstrap over ITERATIONS (supplementary; resamples the 500 DMM fits and
#     captures the algorithm's own run-to-run variability),
#   * the submitted /n ribbon and the correct /sqrt(n) ribbon, for the disclosure.
#
# Inputs : out/cache/erpm_legacy.rds (script 00), the samplesList_alluvials_*.rda
#          core/variable lists (mirror).
# Outputs: out/02_bootstrap_CI/  fig3A_ci.csv fig3B_ci.csv fig4A_ci.csv fig4B_ci.csv,
#          legacy_check.csv, Figure3A.pdf Figure3B.pdf Figure4A.pdf Figure4B.pdf,
#          CI_comparison.pdf, summary.txt
# ------------------------------------------------------------------------------
suppressMessages({ pkgload::load_all(".", quiet = TRUE); library(ggplot2) })
mirror <- Sys.getenv("ENTEROTYPES_MIRROR")
if (!nzchar(mirror) || !dir.exists(mirror)) stop("Set ENTEROTYPES_MIRROR to the directory holding the 2022 analyses (private, not distributed).")
fig <- file.path(mirror, "analyses/9.analyses_mars2026")
out <- "out/02_bootstrap_CI"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
cache <- readRDS("out/cache/erpm_legacy.rds")
ent <- enterotype_names()                      # Bact1 Bact2 Prevotella Ruminococcus = codes 1..4
R_samples <- 2000; R_iter <- 200; seed <- 314

# ---- 0. Fast summaries on the integer-coded assignment matrices --------------------
freq_int <- function(m) {                      # samples x 4 frequencies over iterations drawn
  nd <- rowSums(!is.na(m))
  f <- vapply(1:4, function(k) rowSums(m == k, na.rm = TRUE) / nd, numeric(nrow(m)))
  dimnames(f) <- list(rownames(m), ent); f[nd == 0, ] <- NA; attr(f, "n_drawn") <- nd; f
}
consensus <- function(f, seed = 314) {         # consensus enterotype + robustness, ties seeded
  am <- argmax_rows(f, seed = seed)
  data.frame(sample = rownames(f), enterotype = ent[am$index], robustness = f[cbind(seq_len(nrow(f)), am$index)],
             tie = am$tie, stringsAsFactors = FALSE)[!is.na(am$index), ]
}

# ---- 1. Core / variable lists of the paper (union over the three anchors) ---------
source("inst/scripts/metacardis_data.R")
status_of <- load_published_core_variable(fig)$table[, c("sample", "enterotype", "status")]
status_of <- status_of[!duplicated(status_of[, c("sample", "enterotype")]), ]   # union over the three anchor sets
cat("core/variable entries:", table(status_of$status), "(", length(unique(status_of$sample)), "unique samples )\n")

# ---- 2. Per fraction: robustness by enterotype, status, and pair ---------------------
pairs_list <- list(comb1 = list(c("Bact1", "Bact2"), c("Prevotella", "Ruminococcus")),
                   comb2 = list(c("Bact2", "Ruminococcus"), c("Bact1", "Prevotella")),
                   comb3 = list(c("Prevotella", "Bact2"), c("Bact1", "Ruminococcus")))
ci_row <- function(x, label) {
  b <- bootstrap_ci(x, R = R_samples, type = "bca", seed = seed); p <- bootstrap_ci(x, R = R_samples, type = "percentile", seed = seed)
  ln <- normal_ci(x, "n"); ls <- normal_ci(x, "sqrt_n")
  cbind(label, data.frame(n = length(x), mean = mean(x), sd = stats::sd(x), bca_lower = b[["lower"]], bca_upper = b[["upper"]],
                          perc_lower = p[["lower"]], perc_upper = p[["upper"]],
                          submitted_n_lower = ln[["lower"]], submitted_n_upper = ln[["upper"]],
                          sqrt_n_lower = ls[["lower"]], sqrt_n_upper = ls[["upper"]]))
}
res3A <- res3B <- res4A <- res4B <- list(); legacy_check <- list(); iter_boot <- list()
t0 <- Sys.time()
for (key in cache$fractions$file) {
  Q <- cache$fractions$fraction[cache$fractions$file == key]
  m <- cache$assignments[[key]]
  f <- freq_int(m); cs <- consensus(f)
  # legacy check: stored Prob_Enterotype = E_ProbMax split by consensus enterotype
  leg <- cache$prob_enterotype[[key]]
  if (!is.null(leg) && length(leg) == 4) {
    names(leg) <- c("Bact1", "Prevotella", "Bact2", "Ruminococcus")   # legacy order (Figure3 script line 28)
    legacy_check[[key]] <- do.call(rbind, lapply(ent, function(e) data.frame(fraction = Q, enterotype = e,
      n_new = sum(cs$enterotype == e), n_legacy = length(leg[[e]]),
      mean_new = mean(cs$robustness[cs$enterotype == e]), mean_legacy = mean(leg[[e]]))))
  }
  # 3A: by enterotype
  res3A[[key]] <- do.call(rbind, lapply(ent, function(e) ci_row(cs$robustness[cs$enterotype == e], data.frame(fraction = Q, group = e))))
  # 3B: by enterotype x status (all / core / variable)
  res3B[[key]] <- do.call(rbind, lapply(ent, function(e) {
    x_all <- cs$robustness[cs$enterotype == e]
    rows <- list(ci_row(x_all, data.frame(fraction = Q, group = "all", facet = e)))
    for (st in c("core", "variable")) {
      ids <- status_of$sample[status_of$enterotype == e & status_of$status == st]
      x <- cs$robustness[cs$enterotype == e & cs$sample %in% ids]
      if (length(x) >= 3) rows[[st]] <- ci_row(x, data.frame(fraction = Q, group = st, facet = e))
    }
    do.call(rbind, rows)
  }))
  # 4A / 4B: pairs, algorithm A
  for (cb in names(pairs_list)) {
    mp <- merge_pairs(f[!is.na(f[, 1]), ], pairs = pairs_list[[cb]], method = "A", seed = seed)
    for (pr in unique(mp$pair)) {
      x <- mp$value[mp$pair == pr]
      res4A[[paste(key, cb, pr)]] <- ci_row(x, data.frame(fraction = Q, combination = cb, group = pr))
      if (cb == "comb1") for (st in c("core", "variable")) {
        ids <- status_of$sample[status_of$status == st & status_of$enterotype %in% pairs_list[[cb]][[which(vapply(pairs_list[[cb]], function(z) paste(z, collapse = "") == pr, logical(1)))]]]
        xs <- mp$value[mp$pair == pr & mp$sample %in% ids]
        if (length(xs) >= 3) res4B[[paste(key, pr, st)]] <- ci_row(xs, data.frame(fraction = Q, group = st, facet = pr))
      }
    }
  }
  # bootstrap over iterations: resample the 500 columns, recompute consensus + mean robustness per enterotype
  set.seed(seed)
  bi <- replicate(R_iter, {
    mb <- m[, sample.int(ncol(m), replace = TRUE)]
    fb <- freq_int(mb); cb <- consensus(fb, seed = NULL)
    vapply(ent, function(e) mean(cb$robustness[cb$enterotype == e]), numeric(1))
  })
  iter_boot[[key]] <- data.frame(fraction = Q, group = ent, iter_lower = apply(bi, 1, stats::quantile, 0.025, na.rm = TRUE),
                                 iter_upper = apply(bi, 1, stats::quantile, 0.975, na.rm = TRUE), row.names = NULL)
  cat(sprintf("Q=%.2f done (%.1f min)\n", Q, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
res3A <- do.call(rbind, res3A); res3B <- do.call(rbind, res3B); res4A <- do.call(rbind, res4A); res4B <- do.call(rbind, res4B)
iter_boot <- do.call(rbind, iter_boot); legacy_check <- do.call(rbind, legacy_check)
res3A <- merge(res3A, iter_boot, by = c("fraction", "group"), all.x = TRUE)
for (nm in c("res3A", "res3B", "res4A", "res4B")) { d <- get(nm); d <- d[order(d$fraction, d$group), ]; rownames(d) <- NULL; assign(nm, d) }
write.csv(res3A, file.path(out, "fig3A_ci.csv"), row.names = FALSE); write.csv(res3B, file.path(out, "fig3B_ci.csv"), row.names = FALSE)
write.csv(res4A, file.path(out, "fig4A_ci.csv"), row.names = FALSE); write.csv(res4B, file.path(out, "fig4B_ci.csv"), row.names = FALSE)
write.csv(legacy_check, file.path(out, "legacy_check.csv"), row.names = FALSE)

# ---- 3. Figures -----------------------------------------------------------------------
band <- function(d, lo, hi) { d$lower <- d[[lo]]; d$upper <- d[[hi]]; d }
ggsave(file.path(out, "Figure3A.pdf"), plot_robustness(band(res3A, "bca_lower", "bca_upper")) +
         labs(title = "Figure 3A: mean robustness by enterotype, BCa bootstrap over samples (95%)"), width = 8, height = 5)
ggsave(file.path(out, "Figure3B.pdf"), plot_robustness(band(res3B, "bca_lower", "bca_upper"),
         colours = c(all = "grey40", core = "#1b9e77", variable = "#d95f02")) +
         labs(title = "Figure 3B: core vs variable, BCa bootstrap over samples (95%)"), width = 14, height = 4.5)
p4 <- band(res4A, "bca_lower", "bca_upper"); p4$facet <- p4$combination
ggsave(file.path(out, "Figure4A.pdf"), plot_robustness(p4) + labs(title = "Figure 4A: pairwise groupings (algorithm A), BCa bootstrap (95%)"), width = 14, height = 4.5)
ggsave(file.path(out, "Figure4B.pdf"), plot_robustness(band(res4B, "bca_lower", "bca_upper"),
         colours = c(core = "#1b9e77", variable = "#d95f02")) + labs(title = "Figure 4B: trade-off pairs, core vs variable, BCa bootstrap (95%)"), width = 10, height = 4.5)
# CI comparison for the disclosure: submitted /n, correct /sqrt(n), BCa, iteration bootstrap
cmp <- rbind(
  cbind(kind = "submitted: 1.96 sd / n", band(res3A, "submitted_n_lower", "submitted_n_upper")),
  cbind(kind = "normal: 1.96 sd / sqrt(n)", band(res3A, "sqrt_n_lower", "sqrt_n_upper")),
  cbind(kind = "BCa bootstrap over samples", band(res3A, "bca_lower", "bca_upper")),
  cbind(kind = "bootstrap over iterations", band(res3A, "iter_lower", "iter_upper")))
cmp$facet <- factor(cmp$kind, levels = unique(cmp$kind))
ggsave(file.path(out, "CI_comparison.pdf"), plot_robustness(cmp) + labs(title = "Figure 3A under four interval definitions"), width = 16, height = 4.5)

# ---- 4. Summary -------------------------------------------------------------------------
sink(file.path(out, "summary.txt"))
cat("02  Bootstrap CIs, run", format(Sys.time()), "\n")
cat(sprintf("%d fractions; %d sample-bootstrap replicates (BCa); %d iteration-bootstrap replicates\n\n", nrow(cache$fractions), R_samples, R_iter))
cat("Legacy check (recomputed robustness vs stored Prob_Enterotype; ties broken differently):\n")
lc <- legacy_check; lc$d_mean <- lc$mean_new - lc$mean_legacy; lc$d_n <- lc$n_new - lc$n_legacy
cat(sprintf("  max |d mean| = %.4f, max |d n| = %d, over %d fraction x enterotype cells\n\n", max(abs(lc$d_mean)), max(abs(lc$d_n)), nrow(lc)))
w <- res3A; w$submitted <- w$submitted_n_upper - w$submitted_n_lower; w$sqrt_n <- w$sqrt_n_upper - w$sqrt_n_lower
w$bca <- w$bca_upper - w$bca_lower; w$iter <- w$iter_upper - w$iter_lower
cat("Interval widths (median over fractions) by enterotype:\n")
print(aggregate(cbind(submitted, sqrt_n, bca, iter) ~ group, w, stats::median), row.names = FALSE)
cat(sprintf("\nmedian ratio sqrt(n)/submitted width = %.1f ; BCa upper bound > 1 in %d cells (submitted: %d cells above 1)\n",
            stats::median(w$sqrt_n / w$submitted), sum(w$bca_upper > 1), sum(w$submitted_n_upper > 1)))
cat("\nFirst fraction at which every enterotype's BCa lower bound exceeds 0.95:\n")
q95 <- aggregate(bca_lower ~ fraction, res3A, min); cat("  Q =", q95$fraction[which(q95$bca_lower > 0.95)[1]], "\n")
cat("\nSession:\n"); print(sessionInfo())
sink()
cat(readLines(file.path(out, "summary.txt"), n = 30), sep = "\n")
