#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 01  Tests that respect the shared anchor set (Reviewer 1, major point 5)
#
# The same anchors appear in every ERP condition. The submitted manuscript used
# a Chi-square on the condition x label table and unpaired Kruskal-Wallis + Dunn
# tests on entropy (dunn.test, method = "hochberg", labelled FDR). Here:
#   * entropy: Friedman test with samples as blocks + Conover post hoc, BH,
#     within strata defined per sample (the sample's enterotype when clustered
#     with its own group, the "reference" condition);
#   * label distributions: Chi-square statistic with a within-sample
#     permutation null, cell post hoc on standardised residuals; Cochran's Q per
#     enterotype as the parametric alternative;
#   * the legacy tests recomputed for a side-by-side comparison.
#
# Inputs : the figure objects saved by the 2026 figure scripts (mirror).
# Outputs: out/01_paired_tests/  tableS4_paired.{csv,tex}, label_tests.csv,
#          legacy_vs_paired.csv, Fig1C.pdf Fig1F.pdf Fig2C.pdf Fig2F.pdf, summary.txt
# ------------------------------------------------------------------------------
suppressMessages({ pkgload::load_all(".", quiet = TRUE); library(ggplot2) })
mirror <- Sys.getenv("ENTEROTYPES_MIRROR", "/Users/jdz/PROGRAMMATION/R/ENTEROTYPES")
fig <- file.path(mirror, "analyses/9.analyses_mars2026")
out <- "out/01_paired_tests"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
ent_levels <- c("Bact2", "Prevotella", "Ruminococcus", "Bact1")

# ---- 1. The four anchor designs -----------------------------------------------
load_frame <- function(file, object, reference, cond_levels) {
  e <- new.env(); load(file.path(fig, file), envir = e)
  d <- get(object, e)$data
  d <- data.frame(sample = d$Row.names, condition = as.character(d$Experiment), label = as.character(d$enterotype),
                  Bact1 = d$Bact1, Bact2 = d$Bact2, Prevotella = d$Prevotella, Ruminococcus = d$Ruminococcus,
                  entropy_stored = d$entropy, stringsAsFactors = FALSE)
  freq <- as.matrix(d[, enterotype_names()])
  d$entropy <- assignment_entropy(freq)                       # recomputed, base 2 / log2(4)
  ref <- d[d$condition == reference, ]
  d$reference_label <- ref$label[match(d$sample, ref$sample)]
  d$condition <- factor(d$condition, levels = cond_levels)
  attr(d, "reference") <- reference
  d
}
designs <- list(
  healthy3 = load_frame("Figure1/Figure1_objects.Rda", "plot.alluvial_3experiments", "Healthy-Healthy",
                        c("Healthy-Diabetes", "Healthy-Healthy", "Healthy-Obese")),
  healthy5 = load_frame("Figure1/Figure1_objects.Rda", "plot.alluvial_5experiments", "Healthy-Healthy",
                        c("Healthy-DiabetesT", "Healthy-DiabetesNT", "Healthy-Healthy", "Healthy-ObeseNT", "Healthy-ObeseT")),
  t2d4     = load_frame("Figure2/Figure2ABC_T2Dcore/Figure2t2dCore_plotObjects.rda", "plot.alluvial", "DiabeticsT-DiabeticsT",
                        c("DiabeticsT-DiabeticsNT", "DiabeticsT-DiabeticsT", "DiabeticsT-Healthy", "DiabeticsT-Obese")),
  obese4   = load_frame("Figure2/Figure2DEF_ObeseCore/Figure2_obeseCoreObjects.rda", "plot.alluvial_4experiments", "Obese-Obese",
                        c("Obese-DiabeticsNT", "Obese-DiabeticsT", "Obese-Healthy", "Obese-Obese")))
panel <- c(healthy3 = "Figure 1C", healthy5 = "Figure 1F", t2d4 = "Figure 2C", obese4 = "Figure 2F")

# stored entropy vs recomputed (the stored one exceeds 1 in places)
ent_check <- do.call(rbind, lapply(names(designs), function(k) {
  d <- designs[[k]]; data.frame(design = k, n = nrow(d), max_stored = max(d$entropy_stored), max_recomputed = max(d$entropy),
                                max_abs_diff = max(abs(d$entropy_stored - d$entropy)), cor = stats::cor(d$entropy_stored, d$entropy))
}))

# ---- 2. Entropy: Friedman + Conover, per reference stratum ----------------------
paired <- lapply(designs, function(d) paired_entropy_test(d, stratum = "reference_label"))
tab_s4 <- do.call(rbind, lapply(names(paired), function(k) {
  pw <- paired[[k]]$pairwise; om <- paired[[k]]$omnibus
  pw$design <- k; pw$panel <- panel[[k]]
  pw$friedman_p_adj <- om$p_adj[match(pw$stratum, om$stratum)]
  pw$n <- om$n[match(pw$stratum, om$stratum)]
  pw[, c("design", "panel", "stratum", "n", "friedman_p_adj", "cond1", "cond2", "p", "p_adj")]
}))
tab_s4$stratum <- factor(tab_s4$stratum, levels = ent_levels); tab_s4 <- tab_s4[order(tab_s4$design, tab_s4$stratum), ]
write.csv(tab_s4, file.path(out, "tableS4_paired.csv"), row.names = FALSE)

# LaTeX longtable in the layout of the current Supplementary Table S4
fmt_p <- function(p) ifelse(p < 1e-3, sprintf("$%.2f \\times 10^{%d}$", p / 10^floor(log10(p)), floor(log10(p))), sprintf("%.4f", p))
tex <- c("\\begin{longtable}{p{7.5cm}lrrr}",
  "\\caption{\\textbf{Paired comparisons of anchor entropy across ERP conditions.}",
  "For each anchor set and each enterotype (the anchor's enterotype when clustered with its own clinical group),",
  "a Friedman test with samples as blocks was followed by Conover all-pairs tests; $p$-values were adjusted with the",
  "Benjamini--Hochberg method within each stratum. $n$ is the number of anchors in the stratum.",
  "T/NT: metformin-treated / not treated (T2D) or statin-treated / not treated (obese).}",
  "\\label{tab:paired_entropy} \\\\", "\\toprule",
  "\\textbf{Comparison} & \\textbf{Enterotype} & \\textbf{$n$} & \\textbf{$p$} & \\textbf{FDR} \\\\", "\\midrule", "\\endfirsthead",
  "\\multicolumn{5}{l}{\\small\\textit{(continued from previous page)}} \\\\", "\\toprule",
  "\\textbf{Comparison} & \\textbf{Enterotype} & \\textbf{$n$} & \\textbf{$p$} & \\textbf{FDR} \\\\", "\\midrule", "\\endhead",
  "\\midrule", "\\multicolumn{5}{r}{\\small\\textit{Continued on next page}} \\\\", "\\endfoot", "\\bottomrule", "\\endlastfoot")
for (k in names(designs)) {
  tex <- c(tex, sprintf("\\multicolumn{5}{l}{\\textbf{ERP %s anchors; %s}} \\\\[2pt]",
                        c(healthy3 = "healthy", healthy5 = "healthy", t2d4 = "T2D", obese4 = "obese")[[k]], panel[[k]]))
  tk <- tab_s4[tab_s4$design == k, ]
  tex <- c(tex, sprintf("  %s -- %s & %s & %d & %s & %s \\\\", tk$cond1, tk$cond2, tk$stratum, tk$n, fmt_p(tk$p), fmt_p(tk$p_adj)), "\\\\")
}
tex <- c(tex, "\\end{longtable}")
writeLines(tex, file.path(out, "tableS4_paired.tex"))

# ---- 3. Label distributions: permutation Chi-square + Cochran's Q ---------------
label_tests <- do.call(rbind, lapply(names(designs), function(k) {
  d <- designs[[k]]
  lab <- stats::reshape(d[, c("sample", "condition", "label")], idvar = "sample", timevar = "condition", direction = "wide")
  rownames(lab) <- lab$sample; lab$sample <- NULL; colnames(lab) <- sub("^label\\.", "", colnames(lab))
  lab <- as.matrix(lab)[, levels(d$condition)]
  perm <- anchor_permutation_test(lab, n_perm = 9999)
  leg <- legacy_chisq(lab)
  cq <- lapply(ent_levels, function(e) cochran_q_test(lab, e))
  sig_cells <- perm$cells[perm$cells$p_adj < 0.05, ]
  rbind(
    data.frame(design = k, panel = sub("[CF]$", "B/E", panel[[k]]), test = "Chi-square, independent (submitted)", statistic = leg$statistic, p = leg$p,
               detail = paste(sprintf("%s/%s FDR=%.3g", leg$cells$condition, leg$cells$enterotype, leg$cells$p_adj)[leg$cells$p_adj < 0.05], collapse = "; ")),
    data.frame(design = k, panel = sub("[CF]$", "B/E", panel[[k]]), test = "Chi-square, within-sample permutation", statistic = perm$statistic, p = perm$p,
               detail = paste(sprintf("%s/%s FDR=%.3g", sig_cells$condition, sig_cells$enterotype, sig_cells$p_adj), collapse = "; ")),
    do.call(rbind, lapply(cq, function(q) data.frame(design = k, panel = sub("[CF]$", "B/E", panel[[k]]), test = paste("Cochran Q,", q$test$enterotype),
               statistic = q$test$statistic, p = q$test$p,
               detail = paste(sprintf("%s=%.1f%%", names(q$prevalence), 100 * q$prevalence), collapse = "; ")))))
}))
write.csv(label_tests, file.path(out, "label_tests.csv"), row.names = FALSE)

# ---- 4. Legacy Dunn (unpaired, per-condition stratum, Hochberg) side by side -----
legacy <- do.call(rbind, lapply(names(designs), function(k) {
  d <- designs[[k]]; d$condition <- as.character(d$condition)
  l <- legacy_dunn(d, stratum = "label"); l$design <- k
  l$cond1 <- trimws(sub(" - .*", "", l$comparison)); l$cond2 <- trimws(sub(".* - ", "", l$comparison)); l
}))
key <- function(x) paste(x$design, x$stratum, pmin(x$cond1, x$cond2), pmax(x$cond1, x$cond2))
cmp <- merge(data.frame(key = key(legacy), dunn_p = legacy$p, dunn_hochberg = legacy$p_adj),
             data.frame(key = key(tab_s4), design = tab_s4$design, enterotype = tab_s4$stratum, cond1 = tab_s4$cond1, cond2 = tab_s4$cond2,
                        friedman_conover_p = tab_s4$p, friedman_conover_bh = tab_s4$p_adj), by = "key", all = TRUE)
cmp$key <- NULL
cmp$verdict <- with(cmp, ifelse(dunn_hochberg < 0.05 & friedman_conover_bh < 0.05, "significant in both",
                     ifelse(dunn_hochberg < 0.05, "only submitted", ifelse(friedman_conover_bh < 0.05, "only paired", "neither"))))
write.csv(cmp, file.path(out, "legacy_vs_paired.csv"), row.names = FALSE)

# ---- 5. Panels C and F with brackets from the paired post hoc -------------------
for (k in names(designs)) {
  d <- designs[[k]]; d$reference_label <- factor(d$reference_label, levels = ent_levels)
  pw <- paired[[k]]$pairwise
  p <- plot_entropy(d, value = "entropy", stratum = "reference_label", pairwise = pw, also_raw = FALSE,
                    condition_levels = levels(d$condition)) +
    ggplot2::labs(title = sprintf("%s: anchors stratified by their enterotype in %s", panel[[k]], attr(d, "reference")),
                  subtitle = "Friedman + Conover, ** = BH-adjusted p < 0.05")
  ggplot2::ggsave(file.path(out, paste0(gsub(" ", "", panel[[k]]), ".pdf")), p, width = 11, height = 4)
}

# ---- 6. Summary ------------------------------------------------------------------
sink(file.path(out, "summary.txt"))
cat("01  Paired tests, run", format(Sys.time()), "\n\n")
cat("Entropy stored in the figure objects vs recomputed (base 2, normalised):\n"); print(ent_check, row.names = FALSE)
cat("\nFriedman omnibus per design and reference stratum:\n")
for (k in names(paired)) { cat("--", k, panel[[k]], "\n"); print(paired[[k]]$omnibus, row.names = FALSE) }
cat("\nLabel distribution tests:\n"); print(label_tests[, c("design", "test", "statistic", "p", "detail")], row.names = FALSE)
cat("\nSubmitted Dunn/Hochberg vs paired Friedman/Conover/BH, pair by pair:\n"); print(table(cmp$design, cmp$verdict))
cat("\nSession:\n"); print(sessionInfo())
sink()
cat(readLines(file.path(out, "summary.txt"), n = 60), sep = "\n")
