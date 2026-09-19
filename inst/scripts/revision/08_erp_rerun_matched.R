#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 08  Re-run of the ERP conditions with centroid matching built in (R1, major 1)
#     Confirmatory after script 03 (which re-annotated the stored fits): here every
#     ERP condition is refitted with the package on the master table, with the
#     legacy draws (same anchors, same pool order, seeds 1:500, dmn seed 314,
#     truncation) and `annotate = "match"` against the full-cohort reference fit
#     (Hungarian assignment on Jensen-Shannon distance between centroids), so
#     that no Bacteroides rank is ever used to name a cluster. Compared with the
#     published rank labelling: consensus enterotype per anchor, prevalences,
#     and the core / variable lists under the unanimity rule.
#
# Inputs : METACARDIS_DATA (master table, richness), mirror (anchors, pool order,
#          stored labels), out/cache/reference_fit.rds (script 03)
# Outputs: out/08_erp_rerun_matched/  anchors.csv, prevalence.csv, core_variable.csv,
#          summary.txt ; out/cache/erp_matched_<label>.rds (assignments only)
# Usage  : Rscript inst/scripts/revision/08_erp_rerun_matched.R [n_cores]
# ------------------------------------------------------------------------------
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source("inst/scripts/metacardis_data.R"); source("inst/scripts/cache_utils.R")
mirror <- Sys.getenv("ENTEROTYPES_MIRROR")
if (!nzchar(mirror) || !dir.exists(mirror)) stop("Set ENTEROTYPES_MIRROR to the directory holding the 2022 analyses (private, not distributed).")
args <- commandArgs(trailingOnly = TRUE); n_cores <- if (length(args)) as.integer(args[1]) else 12L
exp_dir <- file.path(mirror, "analyses/2.Estimating_DMM/Article/Experiments")
fig <- file.path(mirror, "analyses/9.analyses_mars2026")
out <- "out/08_erp_rerun_matched"; dir.create(out, showWarnings = FALSE, recursive = TRUE); dir.create("out/cache", showWarnings = FALSE)
ent <- enterotype_names()
md <- load_metacardis_data(); A <- md$abundance
ref <- cache_read("out/cache/reference_fit.rds", reference_fit_manifest(A, md))
if (is.null(ref)) stop("No valid reference fit cache: run inst/scripts/revision/03_cluster_matching.R first.")
status_tab <- load_published_core_variable(fig)$table
anchor_set_of <- c(Healthy = "healthy", DiabeticsT = "t2d", Obese = "obese")

exps <- legacy_erp_experiments()
anchor_rows <- list(); consensus <- list(); timing <- list()
t0 <- Sys.time()
for (i in seq_len(nrow(exps))) {
  lab <- exps$label[i]
  env_file <- list.files(exp_dir, pattern = paste0("^", lab, "_TotalEnvironment\\.RData$"), recursive = TRUE, full.names = TRUE)
  legacy <- read_legacy_erp(exp_dir, lab, keep_fits = FALSE)
  n_comp <- length(legacy$index[[1]]); n_iter <- ncol(legacy$assignments)
  if (length(env_file)) {
    # anchors and pool in the 2022 table order: the draws are then reproduced exactly
    e <- new.env(); load(env_file[1], envir = e)
    anchors <- rownames(e$Fixed_Dataset); pool <- setdiff(rownames(e$Remaining_Dataset), anchors); pool_order <- "2022 table"
    table_ids <- rownames(e$Remaining_Dataset)
  } else {
    # no stored environment: same anchors and same pool, but the pool order (hence the draws) differs
    anchors <- legacy$samples$anchors; pool <- legacy$samples$complementary; pool_order <- "first appearance in the stored draws"
    table_ids <- NULL
  }
  stopifnot(all(c(anchors, pool) %in% rownames(A)))
  cache_file <- file.path("out/cache", paste0("erp_matched_", lab, ".rds"))
  erp_manifest <- cache_manifest(params = list(algorithm = "erp v1, annotate = match", label = lab, n_complementary = n_comp, n_iter = n_iter,
                                               iteration_seeds = "1:n_iter", dmm_seed = 314, integerise = "trunc", pool_order = pool_order),
                                 inputs = list(A = A, anchors = anchors, pool = pool, richness = md$richness[c(anchors, pool)], reference = ref))
  new <- cache_read(cache_file, erp_manifest)
  if (is.null(new)) {
    t1 <- Sys.time()
    r <- erp(A, anchors, pool, n_complementary = n_comp, n_iter = n_iter, marker_features = md$marker_features,
             diversity = md$richness, K = 4, iteration_seeds = seq_len(n_iter), dmm_seed = 314,
             annotate = "match", reference = ref, integerise = "trunc", n_cores = n_cores, keep_fits = FALSE, verbose = FALSE)
    new <- list(assignments = r$assignments, flags = r$flags, gap = r$gap, index = r$index, label_map = r$label_map,
                minutes = as.numeric(difftime(Sys.time(), t1, units = "mins")))
    cache_write(new, cache_file, erp_manifest); rm(r)
  }
  # same draws as 2022?
  same_draws <- if (is.null(table_ids)) NA_real_ else {
    ids_new <- lapply(new$index, function(ix) pool[ix]); ids_old <- lapply(legacy$index, function(ix) table_ids[ix])
    mean(mapply(identical, ids_new, ids_old)) }
  f_new <- assignment_frequency(new$assignments[anchors, ]); f_old <- assignment_frequency(legacy$assignments[anchors, ])
  c_new <- robustness(f_new); c_old <- robustness(f_old)
  aset <- anchor_set_of[[sub("-.*", "", lab)]]; st <- status_tab[status_tab$anchor_set == aset, ]
  anchor_rows[[lab]] <- data.frame(experiment = lab, anchor_set = aset, sample = anchors,
    consensus_published = c_old$enterotype, consensus_matched_rerun = c_new$enterotype,
    robustness_published = c_old$robustness, robustness_matched_rerun = c_new$robustness,
    max_abs_dfreq = apply(abs(f_new - f_old), 1, max),
    status_published = ifelse(anchors %in% st$sample[st$status == "variable"], "variable",
                              ifelse(anchors %in% st$sample[st$status == "core"], "core", "unclassified")))
  consensus[[lab]] <- list(published = stats::setNames(c_old$enterotype, anchors), rerun = stats::setNames(c_new$enterotype, anchors))
  timing[[lab]] <- data.frame(experiment = lab, n_anchors = length(anchors), n_pool = length(pool), n_complementary = n_comp,
                              pool_order = pool_order, same_draws = same_draws, minutes = new$minutes,
                              label_agreement = mean(new$assignments[anchors, ] == legacy$assignments[anchors, ]),
                              consensus_changed = sum(c_old$enterotype != c_new$enterotype),
                              matching_margin_min = min(new$gap$matching_margin_min), flagged = sum(new$flags$bact2_not_lowest_diversity | new$flags$rum_not_highest_diversity))
  cat(sprintf("%-24s draws identical %s, per-iteration label agreement %.1f%%, consensus changed %d / %d  (%.1f min)\n",
              lab, if (is.na(same_draws)) "n/a (pool order unknown)" else sprintf("%.0f%%", 100 * same_draws), 100 * timing[[lab]]$label_agreement, timing[[lab]]$consensus_changed, length(anchors),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
anchors_df <- do.call(rbind, anchor_rows); rownames(anchors_df) <- NULL
write.csv(anchors_df, file.path(out, "anchors.csv"), row.names = FALSE)
timing_df <- do.call(rbind, timing); rownames(timing_df) <- NULL
prev <- do.call(rbind, lapply(names(consensus), function(lab) do.call(rbind, lapply(ent, function(e) data.frame(
  experiment = lab, enterotype = e, prevalence_published = 100 * mean(consensus[[lab]]$published == e),
  prevalence_matched_rerun = 100 * mean(consensus[[lab]]$rerun == e))))))
write.csv(prev, file.path(out, "prevalence.csv"), row.names = FALSE)

# core / variable under the unanimity rule from the matched re-run, per anchor set
cv <- do.call(rbind, lapply(unique(anchors_df$anchor_set), function(a) {
  labs <- unique(anchors_df$experiment[anchors_df$anchor_set == a])
  if (a == "obese") labs <- setdiff(labs, "Obese-Diabetics")     # not part of the published obese design (4 conditions)
  m <- sapply(labs, function(l) consensus[[l]]$rerun); rownames(m) <- names(consensus[[labs[1]]]$rerun)
  cvr <- core_variable(m, "unanimous"); cvr$anchor_set <- a
  pub <- status_tab[status_tab$anchor_set == a, ]
  key_r <- paste(cvr$sample, cvr$enterotype); key_p <- paste(pub$sample, pub$enterotype)
  cvr$status_published <- pub$status[match(key_r, key_p)]
  cvr
}))
write.csv(cv, file.path(out, "core_variable.csv"), row.names = FALSE)

sink(file.path(out, "summary.txt"))
cat("08  ERP re-run with centroid matching, run", format(Sys.time()), "\n\n")
print(timing_df, row.names = FALSE, digits = 3)
cat("\nPrevalence (%) per experiment, published rank labelling vs matched re-run:\n"); print(prev, row.names = FALSE, digits = 3)
cat("\nCore / variable (unanimity) from the matched re-run vs published lists:\n")
print(with(cv, table(anchor_set, rerun = status, published = ifelse(is.na(status_published), "absent", status_published))))
cat("\nAnchors whose consensus enterotype changes, by published status:\n")
print(with(anchors_df, table(status_published, changed = consensus_published != consensus_matched_rerun)))
cat("\nSession:\n"); print(sessionInfo())
sink()
cat(readLines(file.path(out, "summary.txt"), n = 50), sep = "\n")
