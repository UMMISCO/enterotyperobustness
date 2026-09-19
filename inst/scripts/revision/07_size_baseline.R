#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 07  Size baseline for the ERP conditions (Reviewer 1, major 4)
#     Each ERP condition clusters n = anchors + complementary samples (228 healthy,
#     361 T2D-metformin, 334 obese), far below the 2 022 of the full cohort. How
#     much of the anchors' instability is expected from that size alone? The ERPM
#     runs give the answer: the same anchor samples, re-clustered in random
#     stratified draws of the cohort at the fraction Q that matches n
#     (Q = 0.11 -> 222, 0.18 -> 364, 0.17 -> 344 samples). Compared here:
#       (a) each anchor's robustness and entropy under its own-group ERP condition
#           versus under ERPM at the matched Q (paired, Wilcoxon signed rank);
#       (b) the spread between ERP conditions (composition, at fixed n) versus the
#           spread along Q in ERPM (size, random composition);
#       (c) a figure: ERPM robustness of the anchors along n, with the ERP
#           conditions placed at their n.
#
# Inputs : figure objects (mirror), out/cache/erpm_legacy.rds (script 00)
# Outputs: out/07_size_baseline/  anchors_erp_vs_erpm.csv, paired_tests.csv,
#          condition_vs_size_spread.csv, erpm_curve_anchors.csv, Figure_size_baseline.pdf,
#          summary.txt
# ------------------------------------------------------------------------------
suppressMessages({ pkgload::load_all(".", quiet = TRUE); library(ggplot2) })
mirror <- Sys.getenv("ENTEROTYPES_MIRROR")
if (!nzchar(mirror) || !dir.exists(mirror)) stop("Set ENTEROTYPES_MIRROR to the directory holding the 2022 analyses (private, not distributed).")
fig <- file.path(mirror, "analyses/9.analyses_mars2026")
out <- "out/07_size_baseline"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
ent <- enterotype_names()
N_cohort <- 2022

# ---- 1. ERP designs: per-condition frequencies of every anchor ------------------------------
load_design <- function(file, object, cond_levels) {
  e <- new.env(); load(file.path(fig, file), envir = e); d <- get(object, e)$data
  d <- data.frame(sample = d$Row.names, condition = as.character(d$Experiment), label = as.character(d$enterotype),
                  Bact1 = d$Bact1, Bact2 = d$Bact2, Prevotella = d$Prevotella, Ruminococcus = d$Ruminococcus, stringsAsFactors = FALSE)
  d[d$condition %in% cond_levels, ]
}
designs <- list(
  healthy = list(own = "Healthy-Healthy", n_anchor = 142, n_comp = 86,
                 d = load_design("Figure1/Figure1_objects.Rda", "plot.alluvial_5experiments",
                                 c("Healthy-DiabetesT", "Healthy-DiabetesNT", "Healthy-Healthy", "Healthy-ObeseNT", "Healthy-ObeseT"))),
  t2d = list(own = "DiabeticsT-DiabeticsT", n_anchor = 226, n_comp = 135,
             d = load_design("Figure2/Figure2ABC_T2Dcore/Figure2t2dCore_plotObjects.rda", "plot.alluvial",
                             c("DiabeticsT-DiabeticsNT", "DiabeticsT-DiabeticsT", "DiabeticsT-Healthy", "DiabeticsT-Obese"))),
  obese = list(own = "Obese-Obese", n_anchor = 209, n_comp = 125,
               d = load_design("Figure2/Figure2DEF_ObeseCore/Figure2_obeseCoreObjects.rda", "plot.alluvial_4experiments",
                               c("Obese-DiabeticsNT", "Obese-DiabeticsT", "Obese-Healthy", "Obese-Obese"))))
erp_rows <- do.call(rbind, lapply(names(designs), function(a) {
  x <- designs[[a]]; d <- x$d; f <- as.matrix(d[, ent])
  am <- argmax_rows(f, seed = 314)
  data.frame(anchor_set = a, condition = d$condition, sample = d$sample, n_clustered = x$n_anchor + x$n_comp,
             enterotype = ent[am$index], robustness = f[cbind(seq_len(nrow(f)), am$index)],
             entropy = assignment_entropy(f), own = d$condition == x$own, stringsAsFactors = FALSE)
}))

# ---- 2. ERPM: the same anchors along Q --------------------------------------------------------
cache <- readRDS("out/cache/erpm_legacy.rds")
freq_int <- function(m) { nd <- rowSums(!is.na(m)); f <- vapply(1:4, function(k) rowSums(m == k, na.rm = TRUE) / nd, numeric(nrow(m)))
  dimnames(f) <- list(rownames(m), ent); f[nd == 0, ] <- NA; attr(f, "n_drawn") <- nd; f }
anchor_ids <- lapply(designs, function(x) unique(x$d$sample))
erpm_rows <- do.call(rbind, lapply(seq_len(nrow(cache$fractions)), function(i) {
  key <- cache$fractions$file[i]; Q <- cache$fractions$fraction[i]
  f <- freq_int(cache$assignments[[key]]); nd <- attr(f, "n_drawn")
  do.call(rbind, lapply(names(anchor_ids), function(a) {
    ids <- intersect(anchor_ids[[a]], rownames(f)); fa <- f[ids, , drop = FALSE]
    am <- argmax_rows(fa, seed = 314); ok <- !is.na(am$index)
    data.frame(anchor_set = a, fraction = Q, n_clustered = round(Q * N_cohort), sample = ids,
               enterotype = ifelse(ok, ent[am$index], NA), robustness = ifelse(ok, fa[cbind(seq_len(nrow(fa)), am$index)], NA),
               entropy = assignment_entropy(fa), n_drawn = nd[ids], stringsAsFactors = FALSE)
  }))
}))
write.csv(erpm_rows, file.path(out, "erpm_curve_anchors.csv"), row.names = FALSE)

# ---- 3. (a) Own-group ERP condition vs ERPM at the matched fraction, paired by sample -----------
matched_Q <- vapply(designs, function(x) { n <- x$n_anchor + x$n_comp
  cache$fractions$fraction[which.min(abs(cache$fractions$fraction * N_cohort - n))] }, numeric(1))
paired <- do.call(rbind, lapply(names(designs), function(a) {
  Q <- matched_Q[[a]]
  e <- erp_rows[erp_rows$anchor_set == a & erp_rows$own, ]
  m <- erpm_rows[erpm_rows$anchor_set == a & erpm_rows$fraction == Q, ]
  m <- m[match(e$sample, m$sample), ]
  data.frame(anchor_set = a, sample = e$sample, enterotype_erp = e$enterotype, enterotype_erpm = m$enterotype,
             robustness_erp = e$robustness, robustness_erpm = m$robustness, entropy_erp = e$entropy, entropy_erpm = m$entropy,
             n_drawn_erpm = m$n_drawn, matched_fraction = Q, n_erp = e$n_clustered[1], n_erpm = round(Q * N_cohort))
}))
write.csv(paired, file.path(out, "anchors_erp_vs_erpm.csv"), row.names = FALSE)
wtest <- function(x, y) { ok <- !is.na(x) & !is.na(y); if (sum(ok) < 5) return(NA_real_)
  suppressWarnings(stats::wilcox.test(x[ok], y[ok], paired = TRUE, exact = FALSE)$p.value) }
tests <- do.call(rbind, lapply(names(designs), function(a) {
  p <- paired[paired$anchor_set == a, ]
  by_e <- do.call(rbind, lapply(c("all", ent), function(e) {
    q <- if (e == "all") p else p[p$enterotype_erp == e, ]
    data.frame(anchor_set = a, enterotype = e, n = nrow(q), n_drawn_erpm_median = stats::median(q$n_drawn_erpm, na.rm = TRUE),
               robustness_erp = mean(q$robustness_erp, na.rm = TRUE), robustness_erpm = mean(q$robustness_erpm, na.rm = TRUE),
               entropy_erp = mean(q$entropy_erp, na.rm = TRUE), entropy_erpm = mean(q$entropy_erpm, na.rm = TRUE),
               same_consensus = mean(q$enterotype_erp == q$enterotype_erpm, na.rm = TRUE),
               p_robustness = wtest(q$robustness_erp, q$robustness_erpm), p_entropy = wtest(q$entropy_erp, q$entropy_erpm))
  }))
  by_e
}))
write.csv(tests, file.path(out, "paired_tests.csv"), row.names = FALSE)

# ---- 4. (b) Composition spread (between ERP conditions, fixed n) vs size spread (ERPM along Q) ---
cond_means <- aggregate(cbind(robustness, entropy) ~ anchor_set + condition + enterotype, erp_rows, mean)
size_means <- aggregate(cbind(robustness, entropy) ~ anchor_set + fraction + enterotype, erpm_rows, mean)
spread <- do.call(rbind, lapply(names(designs), function(a) do.call(rbind, lapply(ent, function(e) {
  cm <- cond_means[cond_means$anchor_set == a & cond_means$enterotype == e, ]
  sm <- size_means[size_means$anchor_set == a & size_means$enterotype == e, ]
  Q <- matched_Q[[a]]
  at <- function(q) sm$robustness[which.min(abs(sm$fraction - q))]
  data.frame(anchor_set = a, enterotype = e, n_erp = designs[[a]]$n_anchor + designs[[a]]$n_comp, matched_fraction = Q,
             erp_min = min(cm$robustness), erp_max = max(cm$robustness), erp_spread = max(cm$robustness) - min(cm$robustness),
             erp_own = cm$robustness[cm$condition == designs[[a]]$own],
             erpm_at_n = at(Q), erpm_at_0.05 = at(0.05), erpm_at_0.50 = at(0.50), erpm_at_0.95 = at(0.95),
             size_gain_to_full = at(0.95) - at(Q))
}))))
write.csv(spread, file.path(out, "condition_vs_size_spread.csv"), row.names = FALSE)

# ---- 5. Figure: ERPM curve of the anchors along n, ERP conditions at their n ------------------
size_means$anchor_set <- factor(size_means$anchor_set, levels = names(designs)); cond_means$anchor_set <- factor(cond_means$anchor_set, levels = names(designs))
size_means$n <- round(size_means$fraction * N_cohort)
cond_means$n <- vapply(as.character(cond_means$anchor_set), function(a) designs[[a]]$n_anchor + designs[[a]]$n_comp, numeric(1))
cond_means$own <- mapply(function(a, cc) designs[[a]]$own == cc, as.character(cond_means$anchor_set), cond_means$condition)
p <- ggplot() +
  geom_line(data = size_means, aes(n, robustness, colour = enterotype)) +
  geom_point(data = cond_means, aes(n, robustness, colour = enterotype, shape = own), size = 2.4, position = position_jitter(width = 6, height = 0, seed = 1)) +
  geom_vline(data = data.frame(anchor_set = factor(names(designs), levels = names(designs)), n = vapply(designs, function(x) x$n_anchor + x$n_comp, numeric(1))),
             aes(xintercept = n), linetype = "dotted", colour = "grey40") +
  facet_wrap(~anchor_set, nrow = 1, labeller = labeller(anchor_set = c(healthy = "healthy anchors (n = 228)", t2d = "T2D-metformin anchors (n = 361)", obese = "obese anchors (n = 334)"))) +
  scale_colour_manual(values = enterotype_colours()) + scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "own group", `FALSE` = "other condition"), name = "ERP condition") +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey50") +
  labs(x = "number of samples clustered", y = "mean robustness of the anchors",
       title = "Size baseline: ERPM robustness of the anchor samples along cohort size (lines), ERP conditions at their size (points)") +
  theme_bw() + theme(legend.position = "bottom")
ggsave(file.path(out, "Figure_size_baseline.pdf"), p, width = 13, height = 5)

# ---- 6. Summary ----------------------------------------------------------------------------------
sink(file.path(out, "summary.txt"))
cat("07  Size baseline, run", format(Sys.time()), "\n\n")
cat("Matched ERPM fractions:", paste(sprintf("%s: n_ERP = %d -> Q = %.2f (n = %d)", names(designs), vapply(designs, function(x) x$n_anchor + x$n_comp, numeric(1)),
                                            matched_Q, round(matched_Q * N_cohort)), collapse = "; "), "\n\n")
cat("(a) Anchors under their own-group ERP condition (500 iterations) vs under ERPM at the matched fraction\n")
cat("    (each anchor drawn in ~n_drawn of 500 iterations); paired Wilcoxon signed-rank tests:\n")
print(tests, row.names = FALSE, digits = 3)
cat("\n(b) Spread of mean robustness between ERP conditions (composition, fixed n) vs along ERPM (size):\n")
print(spread, row.names = FALSE, digits = 3)
cat("\nReading: 'erp_spread' is what composition does at fixed n; 'size_gain_to_full' is what size does from n_ERP to 95 % of the cohort.\n")
cat("Caveat: ERPM draws vary in composition too (random stratified subsets), so the ERPM baseline is 'size + random composition',\n",
    "and each ERPM robustness rests on ~n_drawn draws, a coarser frequency than the 500 ERP iterations.\n")
cat("\nSession:\n"); print(sessionInfo())
sink()
cat(readLines(file.path(out, "summary.txt"), n = 60), sep = "\n")
