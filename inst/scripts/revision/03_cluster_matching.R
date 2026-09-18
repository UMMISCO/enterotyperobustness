#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 03 + 04  Decouple cluster identity from label (Reviewer 1, major point 1)
#
# On the 500 stored DMM fits of each ERP experiment:
#   * re-annotate every fit by Hungarian matching of its component centroids to a
#     reference fit on the full cohort (Jensen-Shannon distance), instead of the
#     Bacteroides-rank heuristic; compare the two labellings per iteration (ARI,
#     identical name -> component map) and per anchor (consensus enterotype);
#   * record, per iteration, the Bacteroides gap between the Bact1 and Bact2
#     clusters (absolute and relative), the Prevotella margin, the matching margin,
#     and the JSD between the Bact1 and Bact2 centroids vs the closest other pair;
#   * for the anchors, the within-fit posterior margin (assigned minus second
#     component), summarised per anchor and split by the published core/variable
#     status: do variable individuals sit between two distinct structures?
#   * recompute the Bact1/Bact2 and Prevotella/Ruminococcus trade-offs under the
#     matched labels.
#
# Inputs : mirror (stored fits), METACARDIS_DATA (counts of the drawn samples and
#          the reference fit), the samplesList_alluvials_*.rda core/variable lists.
# Outputs: out/03_cluster_matching/  iterations.csv anchors.csv tradeoffs.csv
#          prevalence_rank_vs_matched.csv, Bacteroides_gap.pdf ARI.pdf posterior_margin.pdf,
#          summary.txt ; out/cache/reference_fit.rds
# ------------------------------------------------------------------------------
suppressMessages({ pkgload::load_all(".", quiet = TRUE); library(ggplot2) })
source("inst/scripts/metacardis_data.R")
mirror <- Sys.getenv("ENTEROTYPES_MIRROR", "/Users/jdz/PROGRAMMATION/R/ENTEROTYPES")
exp_dir <- file.path(mirror, "analyses/2.Estimating_DMM/Article/Experiments")
fig <- file.path(mirror, "analyses/9.analyses_mars2026")
out <- "out/03_cluster_matching"; dir.create(out, showWarnings = FALSE, recursive = TRUE); dir.create("out/cache", showWarnings = FALSE)
ent <- enterotype_names()
md <- load_metacardis_data(); A <- md$abundance

# ---- 1. Reference fit on the full cohort (rank-annotated) --------------------------------
ref_file <- "out/cache/reference_fit.rds"
if (file.exists(ref_file)) ref <- readRDS(ref_file) else {
  t0 <- Sys.time()
  ref_fit <- fit_dmm(A, K = 4, seed = 314, integerise = "trunc")
  ref <- reference_centroids(A, md$marker_features, md$richness, fit = ref_fit)
  attr(ref, "fit_time_min") <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  saveRDS(ref, ref_file)
}
ref_ann <- attr(ref, "annotation")
cat("reference fit on", nrow(A), "samples: components", paste(ref_ann$component_names, collapse = "/"),
    "| flags:", paste(names(ref_ann$flags)[unlist(ref_ann$flags)], collapse = ","), "\n")

# ---- 2. Published core / variable status per anchor set -----------------------------------
status_tab <- load_published_core_variable(fig)$table
anchor_set_of <- c(Healthy = "healthy", DiabeticsT = "t2d", Obese = "obese")

# ---- 3. Every experiment, every iteration ------------------------------------------------------
exps <- legacy_erp_experiments()
iter_rows <- list(); anchor_rows <- list(); consensus <- list()
t0 <- Sys.time()
for (i in seq_len(nrow(exps))) {
  lab <- exps$label[i]
  r <- read_legacy_erp(exp_dir, lab, abundance = A, marker_features = md$marker_features, diversity = md$richness, keep_fits = TRUE)
  anchors <- r$samples$anchors; n_iter <- ncol(r$assignments)
  matched <- matrix(NA_character_, nrow(r$assignments), n_iter, dimnames = dimnames(r$assignments))
  margin <- matrix(NA_real_, length(anchors), n_iter, dimnames = list(anchors, NULL))
  it <- vector("list", n_iter)
  for (j in seq_len(n_iter)) {
    fit <- r$fits[[j]]; ids <- rownames(dmm_posterior(fit))
    u <- union_counts(A, ids); ids_u <- u$ids_unique
    am <- annotate_clusters(fit, u$counts, md$marker_features, method = "match", reference = ref)
    comp <- dmm_assign(fit)
    matched[ids_u, j] <- am$component_names[comp][u$keep]
    post <- dmm_posterior(fit)[u$keep, , drop = FALSE]; rownames(post) <- ids_u
    pa <- post[anchors, , drop = FALSE]
    top2 <- t(apply(pa, 1, function(v) sort(v, decreasing = TRUE)[1:2])); margin[, j] <- top2[, 1] - top2[, 2]
    cent <- dmm_centroids(fit); colnames(cent) <- am$component_names
    jsd_b12 <- jsd(cent[, "Bact1"], cent[, "Bact2"])
    others <- utils::combn(ent, 2); others <- others[, !(others[1, ] %in% c("Bact1", "Bact2") & others[2, ] %in% c("Bact1", "Bact2"))]
    jsd_other <- min(apply(others, 2, function(pp) jsd(cent[, pp[1]], cent[, pp[2]])))
    it[[j]] <- data.frame(experiment = lab, iteration = j,
                          map_identical = identical(unname(r$label_map[j, ent]), unname(am$map[ent])),
                          ari = adjusted_rand_index(r$assignments[ids_u, j], matched[ids_u, j]),
                          matching_cost = am$gap$matching_cost, matching_margin = am$gap$matching_margin_min,
                          bact_gap = r$gap$bact_gap[j], bact_gap_relative = r$gap$bact_gap_relative[j],
                          prevotella_margin = r$gap$prevotella_margin[j], prevotella_collision = r$gap$prevotella_collision[j],
                          jsd_bact1_bact2 = jsd_b12, jsd_min_other_pair = jsd_other,
                          flagged_2022 = r$flags[["bact2_not_lowest_diversity"]][j] | r$flags[["rum_not_highest_diversity"]][j])
  }
  iter_rows[[lab]] <- do.call(rbind, it)
  f_rank <- assignment_frequency(r$assignments[anchors, ]); f_match <- assignment_frequency(matched[anchors, ])
  c_rank <- robustness(f_rank); c_match <- robustness(f_match)
  aset <- anchor_set_of[[sub("-.*", "", lab)]]
  st <- status_tab[status_tab$anchor_set == aset, ]
  anchor_rows[[lab]] <- data.frame(experiment = lab, anchor_set = aset, sample = anchors,
    consensus_rank = c_rank$enterotype, consensus_matched = c_match$enterotype,
    robustness_rank = c_rank$robustness, robustness_matched = c_match$robustness,
    max_abs_dfreq = apply(abs(f_rank - f_match), 1, max),
    median_posterior_margin = apply(margin, 1, stats::median), q10_posterior_margin = apply(margin, 1, stats::quantile, 0.10),
    status_published = ifelse(anchors %in% st$sample[st$status == "variable"], "variable",
                              ifelse(anchors %in% st$sample[st$status == "core"], "core", "unclassified")))
  consensus[[lab]] <- list(rank = stats::setNames(c_rank$enterotype, anchors), matched = stats::setNames(c_match$enterotype, anchors))
  cat(sprintf("%-24s %d iterations, maps identical %.1f%%, mean ARI %.3f, consensus changed %d / %d anchors  (%.1f min)\n",
              lab, n_iter, 100 * mean(iter_rows[[lab]]$map_identical), mean(iter_rows[[lab]]$ari, na.rm = TRUE),
              sum(c_rank$enterotype != c_match$enterotype), length(anchors), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  rm(r); invisible(gc())
}
iters <- do.call(rbind, iter_rows); anchors_df <- do.call(rbind, anchor_rows); rownames(iters) <- rownames(anchors_df) <- NULL
write.csv(iters, file.path(out, "iterations.csv"), row.names = FALSE)
write.csv(anchors_df, file.path(out, "anchors.csv"), row.names = FALSE)

# ---- 4. Prevalence and trade-offs under both labellings --------------------------------------
prev <- do.call(rbind, lapply(names(consensus), function(lab) do.call(rbind, lapply(ent, function(e) data.frame(
  experiment = lab, enterotype = e, prevalence_rank = 100 * mean(consensus[[lab]]$rank == e),
  prevalence_matched = 100 * mean(consensus[[lab]]$matched == e))))))
write.csv(prev, file.path(out, "prevalence_rank_vs_matched.csv"), row.names = FALSE)
refs <- c(Healthy = "Healthy-Healthy", DiabeticsT = "DiabeticsT-DiabeticsT", Obese = "Obese-Obese")
trade <- do.call(rbind, lapply(names(consensus), function(lab) {
  ref_lab <- refs[[sub("-.*", "", lab)]]; if (lab == ref_lab) return(NULL)
  do.call(rbind, lapply(c("rank", "matched"), function(m) {
    from <- consensus[[ref_lab]][[m]]; to <- consensus[[lab]][[m]][names(from)]
    mv <- from != to; if (!any(mv)) return(NULL)
    t <- as.data.frame(table(from = from[mv], to = to[mv])); t <- t[t$Freq > 0, ]
    cbind(experiment = lab, labelling = m, reference = ref_lab, t)
  }))
}))
write.csv(trade, file.path(out, "tradeoffs.csv"), row.names = FALSE)

# ---- 5. Figures ----------------------------------------------------------------------------------
iters$experiment <- factor(iters$experiment, levels = exps$label)
ggsave(file.path(out, "Bacteroides_gap.pdf"),
  ggplot(iters, aes(experiment, bact_gap_relative)) + geom_boxplot(outlier.size = 0.5) +
    geom_hline(yintercept = c(0.1, 0.2), linetype = "dashed", colour = "grey50") +
    labs(y = "relative gap in mean Bacteroides, (Bact2 - Bact1) / Bact2", x = NULL,
         title = "Per-iteration Bacteroides gap between the Bact2 and Bact1 clusters (500 fits per experiment)") +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)), width = 11, height = 5)
ggsave(file.path(out, "ARI.pdf"),
  ggplot(iters, aes(experiment, ari)) + geom_boxplot(outlier.size = 0.5) +
    labs(y = "adjusted Rand index, rank vs matched labels", x = NULL, title = "Agreement between the published (rank) and the matched labelling, per iteration") +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)), width = 11, height = 5)
am <- anchors_df[anchors_df$status_published != "unclassified", ]
ggsave(file.path(out, "posterior_margin.pdf"),
  ggplot(am, aes(status_published, median_posterior_margin, fill = status_published)) + geom_boxplot(outlier.size = 0.5) +
    facet_wrap(~experiment, nrow = 2) + scale_fill_manual(values = c(core = "#1b9e77", variable = "#d95f02"), guide = "none") +
    labs(x = NULL, y = "median within-fit posterior margin (assigned - second component)",
         title = "Are variable anchors between distinct structures? Posterior margin per anchor, core vs variable") +
    theme_bw(), width = 14, height = 7)

# ---- 6. Summary -----------------------------------------------------------------------------------
sink(file.path(out, "summary.txt"))
cat("03/04  Cluster matching and Bacteroides gap, run", format(Sys.time()), "\n")
fit_min <- attr(ref, "fit_time_min"); if (is.null(fit_min)) fit_min <- NA_real_
cat(sprintf("reference fit: %d samples, K = 4, seed 314, %.1f min; components %s\n\n", nrow(A), fit_min, paste(ref_ann$component_names, collapse = "/")))
cat("Per experiment: name->component map identical (rank vs matched), ARI, relative Bacteroides gap, JSD Bact1-Bact2 vs closest other pair:\n")
s <- aggregate(cbind(map_identical, ari, bact_gap_relative, jsd_bact1_bact2, jsd_min_other_pair, flagged_2022) ~ experiment, iters, mean)
s$gap_below_0.1 <- aggregate(bact_gap_relative ~ experiment, iters, function(v) mean(v < 0.1))$bact_gap_relative
print(s, row.names = FALSE, digits = 3)
cat("\nAmong iterations whose maps differ: mean relative gap vs identical-map iterations:\n")
print(aggregate(bact_gap_relative ~ map_identical, iters, function(v) c(mean = mean(v), q10 = stats::quantile(v, 0.1))), row.names = FALSE, digits = 3)
cat("\nAnchors whose consensus enterotype differs between rank and matched labelling:\n")
print(aggregate(cbind(changed = consensus_rank != consensus_matched) ~ experiment, anchors_df, function(v) sprintf("%d / %d", sum(v), length(v))), row.names = FALSE)
cat("\nPrevalence (%) under rank vs matched labelling:\n"); print(prev, row.names = FALSE, digits = 4)
cat("\nTrade-offs (anchors whose consensus differs from the own-group condition), rank vs matched:\n"); print(trade, row.names = FALSE)
cat("\nWithin-fit posterior margin (median per anchor) by published status:\n")
print(aggregate(median_posterior_margin ~ status_published + anchor_set, am, function(v) c(median = stats::median(v), q25 = stats::quantile(v, .25), q75 = stats::quantile(v, .75))), row.names = FALSE, digits = 3)
cat("\nSession:\n"); print(sessionInfo())
sink()
cat(readLines(file.path(out, "summary.txt"), n = 60), sep = "\n")
