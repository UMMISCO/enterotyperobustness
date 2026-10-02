#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# reproduce_figures.R  Journal-ready figures of the revision from the outputs of
#     the revision scripts (no data access needed once out/ exists). Panels are
#     drawn with ggplot2, assembled with patchwork, no in-plot titles (captions
#     live in the manuscript), panel letters where the figure has several parts.
#
#   Supplementary Figure S6  ERP-core vs ERP-variable under DMM and under PAM/JSD   (script 09)
#   Supplementary Figure S7  ERPM robustness with K = 3, 4, 5                        (script 10)
#   Supplementary Figure S8  core/variable rule sensitivity (Figure 3B under rules)  (script 06)
#   Figures 3 and 4          BCa ribbons (script 02) and test bars (script 11)       [when 11 exists]
#   Supplementary Figure S4  sub-cohorts, BCa ribbons and tests (script 02b)          [when 02b exists]
#   (Supplementary Figure S2, the Bacteroides gap, is out/03_cluster_matching/Bacteroides_gap.pdf)
#   Figures 1 and 2, panels B to F: produced by script 01b (Figure1_BCEF.pdf, Figure2_BCEF.pdf);
#   the flower panels A and D are added in Inkscape by E.B.
#
# Outputs: out/figures/<name>.pdf ; copied into the manuscript with the compiled
#          numbers by the caller (see inst/scripts/README.md).
# ------------------------------------------------------------------------------
suppressMessages({ library(ggplot2); library(patchwork) })
out <- "out/figures"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
ent <- c("Bact1", "Bact2", "Prevotella", "Ruminococcus")
ent_cols <- c(Bact2 = "#F2C57C", Bact1 = "#E08E45", Prevotella = "#FF6B6B", Ruminococcus = "#6E2C00")
status_cols <- c(core = "#1b9e77", variable = "#d95f02", `all samples` = "grey35")
th <- theme_bw(base_size = 9) + theme(strip.background = element_rect(fill = "grey93"), legend.position = "bottom",
                                      panel.grid.minor = element_blank(), plot.tag = element_text(face = "bold"))
line95 <- geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey45", linewidth = 0.3)

# ---- S5: PAM / JSD ----------------------------------------------------------------------------
if (file.exists("out/09_pam_jsd/core_vs_variable.csv")) {
  cv <- read.csv("out/09_pam_jsd/core_vs_variable.csv")
  long <- rbind(data.frame(cv[, c("clustering", "fraction", "enterotype")], status = "core", mean = cv$mean_core),
                data.frame(cv[, c("clustering", "fraction", "enterotype")], status = "variable", mean = cv$mean_variable))
  long$clustering <- factor(long$clustering, levels = c("dmm", "pam"), labels = c("DMM (Figure 3B)", "PAM on Jensen-Shannon distance"))
  long$enterotype <- factor(long$enterotype, levels = ent)
  p <- ggplot(long, aes(fraction, mean, colour = status)) + geom_line(linewidth = 0.5) + geom_point(size = 0.9) + line95 +
    facet_grid(enterotype ~ clustering) + scale_colour_manual(values = status_cols, name = NULL) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "fraction of the cohort drawn at each ERPM iteration", y = "mean robustness") + th
  ggsave(file.path(out, "FigureS6_PAM_core_variable.pdf"), p, width = 6.5, height = 7.5)
}

# ---- S6: alternative K -------------------------------------------------------------------------
if (file.exists("out/10_alt_K/core_vs_variable_by_K.csv")) {
  cv <- read.csv("out/10_alt_K/core_vs_variable_by_K.csv"); rb <- read.csv("out/10_alt_K/robustness_by_K.csv")
  curve <- aggregate(robustness ~ K + fraction, rb, mean)
  long <- rbind(data.frame(cv[, c("K", "fraction")], status = "core", mean = cv$mean_core),
                data.frame(cv[, c("K", "fraction")], status = "variable", mean = cv$mean_variable),
                data.frame(curve[, c("K", "fraction")], status = "all samples", mean = curve$robustness))
  long$K <- factor(paste0("K = ", long$K)); long$status <- factor(long$status, levels = c("all samples", "core", "variable"))
  p <- ggplot(long, aes(fraction, mean, colour = status)) + geom_line(linewidth = 0.5) + geom_point(size = 0.9) + line95 +
    facet_wrap(~K, nrow = 1) + scale_colour_manual(values = status_cols, name = NULL) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "fraction of the cohort drawn at each ERPM iteration", y = "mean robustness") + th
  ggsave(file.path(out, "FigureS7_alt_K.pdf"), p, width = 7.5, height = 3.2)
}

# ---- S7: criterion sensitivity ---------------------------------------------------------------
if (file.exists("out/06_core_variable/fig3B_by_rule.csv")) {
  f <- read.csv("out/06_core_variable/fig3B_by_rule.csv")
  keep <- c("unanimous", "k_minus_1", "sd<0.05", "sd<0.10", "sd<0.20", "sd<0.30")
  f <- f[f$rule %in% keep, ]
  long <- rbind(data.frame(f[, c("rule", "fraction", "enterotype")], status = "core", mean = f$mean_core),
                data.frame(f[, c("rule", "fraction", "enterotype")], status = "variable", mean = f$mean_variable))
  lab <- c(unanimous = "unanimity (main analysis)", k_minus_1 = "all but one condition", `sd<0.05` = "SD < 0.05", `sd<0.10` = "SD < 0.10", `sd<0.20` = "SD < 0.20", `sd<0.30` = "SD < 0.30")
  long$rule <- factor(lab[long$rule], levels = lab[keep]); long$enterotype <- factor(long$enterotype, levels = ent)
  p <- ggplot(long, aes(fraction, mean, colour = status)) + geom_line(linewidth = 0.45) + geom_point(size = 0.6) + line95 +
    facet_grid(enterotype ~ rule) + scale_colour_manual(values = status_cols, name = NULL) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1), breaks = c(0.25, 0.5, 0.75)) +
    labs(x = "fraction of the cohort drawn at each ERPM iteration", y = "mean robustness") + th +
    theme(strip.text.x = element_text(size = 7))
  ggsave(file.path(out, "FigureS8_criterion_sensitivity.pdf"), p, width = 9, height = 7)
}

# ---- Figures 3 and 4: BCa ribbons + test bars (needs scripts 02 and 11) --------------------------
bars <- function(tests, groups, xcol = "fraction") {   # one row per significant pair, drawn as short vertical ticks under the panel
  tests <- tests[!is.na(tests$p_adj) & tests$p_adj < 0.05, ]
  tests$pair <- factor(paste(tests$g1, "vs", tests$g2)); tests
}
if (file.exists("out/02_bootstrap_CI/fig3A_ci.csv") && file.exists("out/11_erpm_tests/fig3A_tests.csv")) {
  ci3A <- read.csv("out/02_bootstrap_CI/fig3A_ci.csv"); t3A <- read.csv("out/11_erpm_tests/fig3A_tests.csv")
  ci3A$group <- factor(ci3A$group, levels = ent)
  pA <- ggplot(ci3A, aes(fraction, mean, colour = group, fill = group)) +
    geom_ribbon(aes(ymin = bca_lower, ymax = bca_upper), alpha = 0.25, colour = NA) + geom_line(linewidth = 0.5) + geom_point(size = 0.7) + line95 +
    scale_colour_manual(values = ent_cols, name = NULL) + scale_fill_manual(values = ent_cols, name = NULL) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) + coord_cartesian(ylim = c(min(ci3A$bca_lower), 1)) +
    labs(x = NULL, y = "mean robustness") + th + theme(legend.position = "right")
  sig <- bars(t3A); sig$pair <- factor(sig$pair, levels = rev(sort(unique(sig$pair))))
  pAt <- ggplot(sig, aes(fraction, pair)) + geom_point(shape = 15, size = 1.6, colour = "grey25") +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = range(ci3A$fraction)) +
    labs(x = "fraction of the cohort drawn at each ERPM iteration", y = NULL) + th + theme(legend.position = "none", panel.grid.major.y = element_blank())
  ci3B <- read.csv("out/02_bootstrap_CI/fig3B_ci.csv"); t3B <- read.csv("out/11_erpm_tests/fig3B_tests.csv")
  ci3B <- ci3B[ci3B$group != "all", ]; ci3B$facet <- factor(ci3B$facet, levels = ent)
  pB <- ggplot(ci3B, aes(fraction, mean, colour = group, fill = group)) +
    geom_ribbon(aes(ymin = bca_lower, ymax = bca_upper), alpha = 0.25, colour = NA) + geom_line(linewidth = 0.5) + geom_point(size = 0.7) + line95 +
    facet_wrap(~facet, nrow = 1) + scale_colour_manual(values = status_cols, name = NULL) + scale_fill_manual(values = status_cols, name = NULL) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) + coord_cartesian(ylim = c(min(ci3B$bca_lower), 1)) +
    labs(x = NULL, y = "mean robustness") + th
  sigB <- t3B[t3B$core_higher_adj, ]; sigB$enterotype <- factor(sigB$enterotype, levels = ent)
  pBt <- ggplot(sigB, aes(fraction, "core > variable")) + geom_point(shape = 15, size = 1.6, colour = "grey25") + facet_wrap(~enterotype, nrow = 1) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = range(ci3B$fraction)) +
    labs(x = "fraction of the cohort drawn at each ERPM iteration", y = NULL) + th + theme(strip.text = element_blank(), panel.grid.major.y = element_blank())
  fig3 <- wrap_plots(pA, pB, pAt, pBt, ncol = 2, widths = c(1, 2.2), heights = c(4, 1.4)) +
    plot_annotation(tag_levels = list(c("A", "B", "", "")))
  ggsave(file.path(out, "Figure3_ERPM.pdf"), fig3, width = 12, height = 5)
}
if (file.exists("out/02_bootstrap_CI/fig4A_ci.csv") && file.exists("out/11_erpm_tests/fig4A_tests.csv")) {
  ci4A <- read.csv("out/02_bootstrap_CI/fig4A_ci.csv"); t4A <- read.csv("out/11_erpm_tests/fig4A_tests.csv"); w4A <- read.csv("out/11_erpm_tests/fig4A_within.csv")
  comb_lab <- c(comb1 = "Bact1+Bact2 / Prevotella+Ruminococcus", comb2 = "Bact2+Ruminococcus / Bact1+Prevotella", comb3 = "Prevotella+Bact2 / Bact1+Ruminococcus")
  ci4A$combination <- factor(comb_lab[ci4A$combination], levels = comb_lab)
  pA <- ggplot(ci4A, aes(fraction, mean, colour = group, fill = group)) +
    geom_ribbon(aes(ymin = bca_lower, ymax = bca_upper), alpha = 0.25, colour = NA) + geom_line(linewidth = 0.5) + geom_point(size = 0.7) + line95 +
    facet_wrap(~combination, nrow = 1) + scale_x_continuous(labels = scales::percent_format(accuracy = 1)) + coord_cartesian(ylim = c(min(ci4A$bca_lower), 1)) +
    labs(x = NULL, y = "mean robustness of the pair", colour = NULL, fill = NULL) + th + guides(colour = guide_legend(nrow = 2), fill = guide_legend(nrow = 2))
  sig <- bars(t4A); sig$pair <- gsub("comb", "grouping ", sig$pair)
  pAt <- ggplot(sig, aes(fraction, pair)) + geom_point(shape = 15, size = 1.6, colour = "grey25") +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = range(ci4A$fraction)) +
    labs(x = "fraction of the cohort drawn at each ERPM iteration", y = NULL) + th + theme(panel.grid.major.y = element_blank())
  ci4B <- read.csv("out/02_bootstrap_CI/fig4B_ci.csv"); t4B <- read.csv("out/11_erpm_tests/fig4B_tests.csv")
  pB <- ggplot(ci4B, aes(fraction, mean, colour = group, fill = group)) +
    geom_ribbon(aes(ymin = bca_lower, ymax = bca_upper), alpha = 0.25, colour = NA) + geom_line(linewidth = 0.5) + geom_point(size = 0.7) + line95 +
    facet_wrap(~facet, nrow = 1) + scale_colour_manual(values = status_cols, name = NULL) + scale_fill_manual(values = status_cols, name = NULL) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) + coord_cartesian(ylim = c(min(ci4B$bca_lower), 1)) +
    labs(x = NULL, y = "mean robustness of the pair") + th
  g4B <- read.csv("out/11_erpm_tests/fig4B_groupings.csv")
  sigB <- g4B[g4B$comb1_higher_adj, ]; sigB$other <- ifelse(sigB$g1 == "comb1", sigB$g2, sigB$g1)
  sigB$row <- paste0("trade-off grouping > ", gsub("comb", "grouping ", sigB$other), " (", sigB$status, ")")
  sigB$row <- factor(sigB$row, levels = rev(sort(unique(sigB$row))))
  pBt <- ggplot(sigB, aes(fraction, row)) + geom_point(shape = 15, size = 1.6, colour = "grey25") +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = range(ci4B$fraction)) +
    labs(x = "fraction of the cohort drawn at each ERPM iteration", y = NULL) + th + theme(panel.grid.major.y = element_blank())
  fig4 <- wrap_plots(pA, pAt, pB, pBt, ncol = 1, heights = c(4, 1.2, 4, 1.4)) + plot_annotation(tag_levels = list(c("A", "", "B", "")))
  ggsave(file.path(out, "Figure4_pairs.pdf"), fig4, width = 10, height = 9)
}

# ---- S3: sub-cohorts (needs 02b) ------------------------------------------------------------------
if (file.exists("out/02b_subcohorts/figS3_ci.csv")) {
  ci <- read.csv("out/02b_subcohorts/figS3_ci.csv"); t <- read.csv("out/02b_subcohorts/figS3_tests.csv")
  ci$group <- factor(ci$group, levels = ent)
  pA <- ggplot(ci, aes(fraction, mean, colour = subcohort, fill = subcohort)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, colour = NA) + geom_line(linewidth = 0.5) + geom_point(size = 0.6) + line95 +
    facet_wrap(~group, nrow = 1) + scale_x_continuous(labels = scales::percent_format(accuracy = 1)) + coord_cartesian(ylim = c(min(ci$lower), 1)) +
    labs(x = NULL, y = "mean robustness", colour = "sub-cohort", fill = "sub-cohort") + th
  pB <- ggplot(ci, aes(fraction, mean, colour = group, fill = group)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25, colour = NA) + geom_line(linewidth = 0.5) + geom_point(size = 0.6) + line95 +
    facet_wrap(~subcohort, nrow = 1) + scale_colour_manual(values = ent_cols, name = NULL) + scale_fill_manual(values = ent_cols, name = NULL) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) + coord_cartesian(ylim = c(min(ci$lower), 1)) +
    labs(x = "fraction of the sub-cohort drawn at each ERPM iteration", y = "mean robustness") + th
  figS3 <- wrap_plots(pA, pB, ncol = 1) + plot_annotation(tag_levels = "A")
  ggsave(file.path(out, "FigureS4_subcohorts.pdf"), figS3, width = 11, height = 7)
}
cat("figures written to", out, ":\n"); print(list.files(out))
