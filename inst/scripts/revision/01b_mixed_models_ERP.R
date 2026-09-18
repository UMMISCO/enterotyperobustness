#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 01b  Mixed-effects analyses of the ERP designs (Reviewer 1, major point 5)
#      Primary analysis, as decided with E. Belda (18 Sep 2026):
#        entropy : lmer(entropy ~ experiment * label + (1 | individual)), type III,
#                  emmeans pairwise contrasts between experiments within label, BH
#        labels  : mblogit(enterotype ~ experiment, random = ~1 | individual), joint Wald;
#                  fallback MQL -> within-individual permutation of a multinom Wald
#                  statistic -> tests of script 01
#      Friedman/Conover by reference stratum (script 01) is the sensitivity check.
#
# Inputs : the saved figure objects (mirror) — same as script 01.
# Outputs: out/01b_mixed_models/  mixed_entropy_pairwise.csv, mixed_entropy_omnibus.csv,
#          mixed_labels.csv, three_way_comparison.csv, tableS4_mixed.tex,
#          Figure1C/1F/2C/2F.pdf (LMM brackets), Figure1_BCEF.pdf, Figure2_BCEF.pdf, summary.txt
# ------------------------------------------------------------------------------
suppressMessages({ pkgload::load_all(".", quiet = TRUE); library(ggplot2) })
mirror <- Sys.getenv("ENTEROTYPES_MIRROR", "/Users/jdz/PROGRAMMATION/R/ENTEROTYPES")
fig <- file.path(mirror, "analyses/9.analyses_mars2026")
out <- "out/01b_mixed_models"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
ent_levels <- c("Bact2", "Prevotella", "Ruminococcus", "Bact1")

load_frame <- function(file, object, reference, cond_levels) {
  e <- new.env(); load(file.path(fig, file), envir = e)
  d <- get(object, e)$data
  d <- data.frame(sample = d$Row.names, condition = as.character(d$Experiment), label = as.character(d$enterotype),
                  Bact1 = d$Bact1, Bact2 = d$Bact2, Prevotella = d$Prevotella, Ruminococcus = d$Ruminococcus,
                  stringsAsFactors = FALSE)
  d$entropy <- assignment_entropy(as.matrix(d[, enterotype_names()]))
  d$condition <- factor(d$condition, levels = cond_levels)
  d$label <- factor(d$label, levels = ent_levels)
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
alluvial_panel <- c(healthy3 = "Figure 1B", healthy5 = "Figure 1E", t2d4 = "Figure 2B", obese4 = "Figure 2E")
anchor_name <- c(healthy3 = "healthy", healthy5 = "healthy", t2d4 = "T2D", obese4 = "obese")

# ---- 1. Entropy: LMM with per-condition label ------------------------------------
lmm <- lapply(designs, function(d) mixed_entropy_model(d, label = "label"))
omni <- do.call(rbind, lapply(names(lmm), function(k) cbind(design = k, panel = panel[[k]], lmm[[k]]$omnibus,
                                                             p_rank = lmm[[k]]$omnibus_rank$p)))
pw <- do.call(rbind, lapply(names(lmm), function(k) cbind(design = k, panel = panel[[k]], lmm[[k]]$pairwise)))
pw$stratum <- factor(pw$stratum, levels = ent_levels); pw <- pw[order(pw$design, pw$stratum), ]
write.csv(omni, file.path(out, "mixed_entropy_omnibus.csv"), row.names = FALSE)
write.csv(pw, file.path(out, "mixed_entropy_pairwise.csv"), row.names = FALSE)

# ---- 2. Labels: mixed multinomial logit, with fallbacks ---------------------------
label_matrix <- function(d) {
  lab <- stats::reshape(d[, c("sample", "condition", "label")], idvar = "sample", timevar = "condition", direction = "wide")
  rownames(lab) <- lab$sample; lab$sample <- NULL; colnames(lab) <- sub("^label\\.", "", colnames(lab))
  m <- as.matrix(lab)[, levels(d$condition)]; mode(m) <- "character"; m
}
label_tests <- do.call(rbind, lapply(names(designs), function(k) {
  d <- designs[[k]]; lab <- label_matrix(d); ref <- attr(d, "reference")
  r <- multinomial_mixed_test(lab, reference_condition = ref, method = "PQL")
  used <- "mblogit PQL"
  if (!r$converged) { r <- multinomial_mixed_test(lab, reference_condition = ref, method = "MQL"); used <- "mblogit MQL" }
  if (r$converged) {
    per <- r$per_enterotype
    data.frame(design = k, panel = alluvial_panel[[k]], method = used, statistic = r$wald$statistic, df = r$wald$df, p = r$wald$p,
               detail = paste(sprintf("%s: W=%.1f p_adj=%.2g", per$enterotype, per$statistic, per$p_adj), collapse = "; "),
               messages = paste(r$messages, collapse = " | "))
  } else {
    pt <- multinomial_permutation_test(lab, reference_condition = ref, n_perm = 999)
    data.frame(design = k, panel = alluvial_panel[[k]], method = "multinom Wald, within-individual permutation (fallback)",
               statistic = pt$statistic, df = pt$df, p = pt$p, detail = sprintf("%d permutations", pt$n_perm),
               messages = paste(r$messages, collapse = " | "))
  }
}))
write.csv(label_tests, file.path(out, "mixed_labels.csv"), row.names = FALSE)
# coefficient tables (log-odds vs own-group condition) where the mixed logit converged
coef_tables <- lapply(names(designs), function(k) {
  d <- designs[[k]]; r <- multinomial_mixed_test(label_matrix(d), reference_condition = attr(d, "reference"))
  if (r$converged) cbind(design = k, r$coefficients) else NULL
})
coef_tables <- do.call(rbind, coef_tables)
if (!is.null(coef_tables)) write.csv(coef_tables, file.path(out, "mixed_labels_coefficients.csv"), row.names = FALSE)

# ---- 3. Three-way comparison with scripts 01 outputs -----------------------------
s01 <- "out/01_paired_tests"
three <- NULL
if (file.exists(file.path(s01, "legacy_vs_paired.csv"))) {
  prev <- read.csv(file.path(s01, "legacy_vs_paired.csv"), stringsAsFactors = FALSE)
  key <- function(x, a, b, e) paste(x$design, x[[e]], pmin(x[[a]], x[[b]]), pmax(x[[a]], x[[b]]))
  prev$key <- key(prev, "cond1", "cond2", "enterotype")
  pwk <- pw; pwk$key <- key(pwk, "cond1", "cond2", "stratum")
  three <- merge(prev, data.frame(key = pwk$key, lmm_p = pwk$p, lmm_bh = pwk$p_adj), by = "key", all = TRUE)
  three$key <- NULL
  three$lmm_verdict <- ifelse(three$lmm_bh < 0.05, "significant", "not")
  write.csv(three, file.path(out, "three_way_comparison.csv"), row.names = FALSE)
}

# ---- 4. Table S4 replacement (LMM contrasts) -------------------------------------
fmt_p <- function(p) ifelse(is.na(p), "--", ifelse(p < 1e-3, sprintf("$%.2f \\times 10^{%d}$", p / 10^floor(log10(p)), floor(log10(p))), sprintf("%.4f", p)))
tex <- c("\\begin{longtable}{p{7.2cm}lrrrr}",
  "\\caption{\\textbf{Pairwise differences in anchor entropy between ERP conditions (linear mixed model).}",
  "For each anchor set, entropy was modelled as experiment $\\times$ enterotype with a random intercept per individual;",
  "contrasts between experiments within each enterotype (the enterotype assigned in that experiment) are estimated marginal-mean",
  "differences with Satterthwaite degrees of freedom, adjusted with the Benjamini--Hochberg method within enterotype.",
  "T/NT: metformin-treated / not treated (T2D) or statin-treated / not treated (obese).}",
  "\\label{tab:mixed_entropy} \\\\", "\\toprule",
  "\\textbf{Comparison} & \\textbf{Enterotype} & \\textbf{Difference} & \\textbf{SE} & \\textbf{$p$} & \\textbf{FDR} \\\\", "\\midrule", "\\endfirsthead",
  "\\multicolumn{6}{l}{\\small\\textit{(continued from previous page)}} \\\\", "\\toprule",
  "\\textbf{Comparison} & \\textbf{Enterotype} & \\textbf{Difference} & \\textbf{SE} & \\textbf{$p$} & \\textbf{FDR} \\\\", "\\midrule", "\\endhead",
  "\\midrule", "\\multicolumn{6}{r}{\\small\\textit{Continued on next page}} \\\\", "\\endfoot", "\\bottomrule", "\\endlastfoot")
for (k in names(designs)) {
  tex <- c(tex, sprintf("\\multicolumn{6}{l}{\\textbf{ERP %s anchors; %s}} \\\\[2pt]", anchor_name[[k]], panel[[k]]))
  tk <- pw[pw$design == k, ]
  tex <- c(tex, sprintf("  %s -- %s & %s & %.3f & %.3f & %s & %s \\\\", tk$cond1, tk$cond2, tk$stratum, tk$estimate, tk$se, fmt_p(tk$p), fmt_p(tk$p_adj)), "\\\\")
}
writeLines(c(tex, "\\end{longtable}"), file.path(out, "tableS4_mixed.tex"))

# ---- 5. Panels with LMM brackets, and the assembled B/C/E/F pages -----------------
ent_plots <- list(); all_plots <- list()
for (k in names(designs)) {
  d <- designs[[k]]
  pe <- plot_entropy(d, stratum = "label", pairwise = lmm[[k]]$pairwise, condition_levels = levels(d$condition)) +
    ggplot2::labs(title = NULL)
  ent_plots[[k]] <- pe
  ggplot2::ggsave(file.path(out, paste0(gsub(" ", "", panel[[k]]), ".pdf")),
                  pe + ggplot2::labs(title = panel[[k]], subtitle = "LMM: experiment x enterotype + (1 | individual); ** BH-adjusted p < 0.05"),
                  width = 11, height = 4)
  all_plots[[k]] <- plot_alluvial(d, label = "label", condition_levels = levels(d$condition),
                                  ylab = paste(length(unique(d$sample)), "anchor samples"))
}
fig1 <- plot_figure_panels(all_plots$healthy3, ent_plots$healthy3, all_plots$healthy5, ent_plots$healthy5)
fig2 <- plot_figure_panels(all_plots$t2d4, ent_plots$t2d4, all_plots$obese4, ent_plots$obese4)
ggplot2::ggsave(file.path(out, "Figure1_BCEF.pdf"), fig1, width = 13.6, height = 11)
ggplot2::ggsave(file.path(out, "Figure2_BCEF.pdf"), fig2, width = 13.6, height = 11)

# ---- 6. Summary ----------------------------------------------------------------------
sink(file.path(out, "summary.txt"))
cat("01b  Mixed-effects analyses, run", format(Sys.time()), "\n")
cat("Stars on the panels come from the LMM contrasts (experiment x label + (1 | individual)); Friedman/Conover by\n",
    "reference stratum (script 01) is the non-parametric sensitivity check.\n\n")
for (k in names(lmm)) { cat("--", k, panel[[k]], "\n"); print(lmm[[k]]); cat("   type III p on rank-transformed entropy:", format(lmm[[k]]$omnibus_rank$p, digits = 3), "\n") }
cat("\nLabel distributions (mixed multinomial logit or fallback):\n"); print(label_tests[, c("design", "method", "statistic", "df", "p", "detail")], row.names = FALSE)
if (any(nzchar(label_tests$messages))) { cat("\nmessages:\n"); print(label_tests[nzchar(label_tests$messages), c("design", "messages")], row.names = FALSE) }
if (!is.null(three)) {
  cat("\nPair-by-pair verdicts (submitted Dunn/Hochberg vs Friedman/Conover vs LMM/emmeans):\n")
  print(table(design = three$design, submitted = three$dunn_hochberg < 0.05, friedman = three$friedman_conover_bh < 0.05, lmm = three$lmm_verdict))
}
cat("\nSession:\n"); print(sessionInfo())
sink()
cat(readLines(file.path(out, "summary.txt"), n = 70), sep = "\n")
