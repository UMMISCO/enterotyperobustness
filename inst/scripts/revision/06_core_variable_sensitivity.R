#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 06  Sensitivity of the core / variable classification (Reviewer 1, minor 1)
#     Rules compared: "unanimous" (published: same enterotype in every ERP condition
#     of the anchor), "k_minus_1" (all but one), and the continuous criterion the
#     submitted Methods described (core if the SD of each enterotype's assignment
#     frequency across conditions is below a threshold), thresholds 0.05 to 0.40.
#     For each rule: how many anchors are core / variable per enterotype, and the
#     Figure 3B contrast (ERPM robustness of core vs variable) across fractions.
#
# Inputs : figure objects (mirror), out/cache/erpm_legacy.rds (script 00)
# Outputs: out/06_core_variable/  core_variable_counts.csv, fig3B_by_rule.csv,
#          overlap_with_published.csv, Figure3B_sensitivity.pdf, summary.txt
# ------------------------------------------------------------------------------
suppressMessages({ pkgload::load_all(".", quiet = TRUE); library(ggplot2) })
mirror <- Sys.getenv("ENTEROTYPES_MIRROR", "/Users/jdz/PROGRAMMATION/R/ENTEROTYPES")
fig <- file.path(mirror, "analyses/9.analyses_mars2026")
out <- "out/06_core_variable"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
ent <- enterotype_names()

# ---- 1. Anchor designs: per-condition labels and frequencies ---------------------------
load_design <- function(file, object, cond_levels) {
  e <- new.env(); load(file.path(fig, file), envir = e); d <- get(object, e)$data
  d <- data.frame(sample = d$Row.names, condition = as.character(d$Experiment), label = as.character(d$enterotype),
                  Bact1 = d$Bact1, Bact2 = d$Bact2, Prevotella = d$Prevotella, Ruminococcus = d$Ruminococcus, stringsAsFactors = FALSE)
  samples <- unique(d$sample)
  labels <- sapply(cond_levels, function(cc) { x <- d[d$condition == cc, ]; stats::setNames(x$label, x$sample)[samples] })
  labels <- matrix(labels, length(samples), length(cond_levels), dimnames = list(samples, cond_levels))
  freqs <- lapply(cond_levels, function(cc) { x <- d[d$condition == cc, ]; m <- as.matrix(x[, ent]); rownames(m) <- x$sample; m[samples, ] })
  list(labels = labels, freqs = freqs)
}
designs <- list(
  healthy = load_design("Figure1/Figure1_objects.Rda", "plot.alluvial_5experiments",
                        c("Healthy-DiabetesT", "Healthy-DiabetesNT", "Healthy-Healthy", "Healthy-ObeseNT", "Healthy-ObeseT")),
  t2d = load_design("Figure2/Figure2ABC_T2Dcore/Figure2t2dCore_plotObjects.rda", "plot.alluvial",
                    c("DiabeticsT-DiabeticsNT", "DiabeticsT-DiabeticsT", "DiabeticsT-Healthy", "DiabeticsT-Obese")),
  obese = load_design("Figure2/Figure2DEF_ObeseCore/Figure2_obeseCoreObjects.rda", "plot.alluvial_4experiments",
                      c("Obese-DiabeticsNT", "Obese-DiabeticsT", "Obese-Healthy", "Obese-Obese")))

# ---- 2. Rules ---------------------------------------------------------------------------------
rules <- c(list(list(name = "unanimous", rule = "unanimous", threshold = NA), list(name = "k_minus_1", rule = "k_minus_1", threshold = NA)),
           lapply(c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40), function(t) list(name = sprintf("sd<%.2f", t), rule = "sd", threshold = t)))
classify <- function(rl) do.call(rbind, lapply(names(designs), function(a) {
  d <- designs[[a]]
  cv <- if (rl$rule == "sd") core_variable(d$labels, "sd", frequencies = d$freqs, threshold = rl$threshold) else core_variable(d$labels, rl$rule)
  cbind(rule = rl$name, anchor = a, cv[, c("sample", "enterotype", "status")])
}))
cls <- lapply(rules, classify); names(cls) <- vapply(rules, `[[`, character(1), "name")
counts <- do.call(rbind, lapply(names(cls), function(r) {
  x <- cls[[r]]; t <- as.data.frame(table(anchor = x$anchor, enterotype = x$enterotype, status = x$status))
  cbind(rule = r, stats::reshape(t, idvar = c("anchor", "enterotype"), timevar = "status", direction = "wide"))
}))
names(counts) <- sub("^Freq\\.", "", names(counts)); counts$rule <- factor(counts$rule, levels = names(cls))
write.csv(counts, file.path(out, "core_variable_counts.csv"), row.names = FALSE)

# overlap of each rule's (sample, enterotype, status) with the published rule
pub <- cls$unanimous; pub_key <- paste(pub$sample, pub$enterotype)
overlap <- do.call(rbind, lapply(names(cls), function(r) {
  x <- cls[[r]]; k <- paste(x$sample, x$enterotype)
  common <- intersect(k, pub_key)
  agree <- mean(x$status[match(common, k)] == pub$status[match(common, pub_key)])
  data.frame(rule = r, n_entries = nrow(x), n_core = sum(x$status == "core"), n_variable = sum(x$status == "variable"),
             shared_with_published = length(common), status_agreement = agree)
}))
write.csv(overlap, file.path(out, "overlap_with_published.csv"), row.names = FALSE)

# ---- 3. Figure 3B contrast under each rule, from the ERPM cache ---------------------------------
cache <- readRDS("out/cache/erpm_legacy.rds")
freq_int <- function(m) { nd <- rowSums(!is.na(m)); f <- vapply(1:4, function(k) rowSums(m == k, na.rm = TRUE) / nd, numeric(nrow(m)))
  dimnames(f) <- list(rownames(m), ent); f[nd == 0, ] <- NA; f }
rob <- lapply(cache$fractions$file, function(k) { f <- freq_int(cache$assignments[[k]]); am <- argmax_rows(f, seed = 314)
  data.frame(sample = rownames(f), enterotype = ent[am$index], robustness = f[cbind(seq_len(nrow(f)), am$index)])[!is.na(am$index), ] })
names(rob) <- cache$fractions$file
fig3b <- do.call(rbind, lapply(names(cls), function(r) {
  st <- cls[[r]]; st <- st[!duplicated(st[, c("sample", "enterotype")]), ]
  do.call(rbind, lapply(names(rob), function(k) {
    Q <- cache$fractions$fraction[cache$fractions$file == k]; x <- rob[[k]]
    do.call(rbind, lapply(ent, function(e) {
      ids_c <- st$sample[st$enterotype == e & st$status == "core"]; ids_v <- st$sample[st$enterotype == e & st$status == "variable"]
      rc <- x$robustness[x$enterotype == e & x$sample %in% ids_c]; rv <- x$robustness[x$enterotype == e & x$sample %in% ids_v]
      p <- if (length(rc) >= 3 && length(rv) >= 3) stats::wilcox.test(rc, rv, exact = FALSE)$p.value else NA
      data.frame(rule = r, fraction = Q, enterotype = e, n_core = length(rc), n_variable = length(rv),
                 mean_core = mean(rc), mean_variable = mean(rv), difference = mean(rc) - mean(rv), p_wilcoxon = p)
    }))
  }))
}))
fig3b$rule <- factor(fig3b$rule, levels = names(cls))
write.csv(fig3b, file.path(out, "fig3B_by_rule.csv"), row.names = FALSE)

# ---- 4. Figure and summary ----------------------------------------------------------------------
long <- rbind(data.frame(fig3b[, c("rule", "fraction", "enterotype")], status = "core", mean = fig3b$mean_core),
              data.frame(fig3b[, c("rule", "fraction", "enterotype")], status = "variable", mean = fig3b$mean_variable))
p <- ggplot(long, aes(fraction, mean, colour = status)) + geom_line() + geom_point(size = 0.6) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey40") +
  facet_grid(enterotype ~ rule) + scale_colour_manual(values = c(core = "#1b9e77", variable = "#d95f02")) +
  labs(x = "fraction of the cohort", y = "mean robustness", title = "Figure 3B under each core / variable rule") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6))
ggsave(file.path(out, "Figure3B_sensitivity.pdf"), p, width = 18, height = 9)
q95 <- aggregate(fraction ~ rule + enterotype, fig3b[fig3b$mean_core >= 0.95, ], min)
sig <- aggregate(p_wilcoxon ~ rule + enterotype, fig3b, function(v) sprintf("%d/%d", sum(v < 0.05, na.rm = TRUE), sum(!is.na(v))))
sink(file.path(out, "summary.txt"))
cat("06  Core / variable sensitivity, run", format(Sys.time()), "\n\n")
cat("Entries per rule (union of the three anchors) and agreement with the published rule on shared (sample, enterotype) pairs:\n")
print(overlap, row.names = FALSE, digits = 3)
cat("\nCore / variable counts per anchor and enterotype:\n"); print(counts, row.names = FALSE)
cat("\nFirst fraction at which mean core robustness reaches 0.95, per rule and enterotype:\n")
print(stats::reshape(q95, idvar = "rule", timevar = "enterotype", direction = "wide"), row.names = FALSE)
cat("\nFractions (of 55) at which core > variable robustness is significant (Wilcoxon p < 0.05):\n")
print(stats::reshape(sig, idvar = "rule", timevar = "enterotype", direction = "wide"), row.names = FALSE)
cat("\nMean difference core - variable at Q = 0.25 / 0.50 / 0.75:\n")
print(fig3b[fig3b$fraction %in% c(0.25, 0.5, 0.75), c("rule", "fraction", "enterotype", "n_core", "n_variable", "difference")], row.names = FALSE, digits = 3)
cat("\nSession:\n"); print(sessionInfo())
sink()
cat(readLines(file.path(out, "summary.txt"), n = 40), sep = "\n")
