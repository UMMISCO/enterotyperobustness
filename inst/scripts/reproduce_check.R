#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Reproduction check (todo item E4): does erp() with the legacy seeds reproduce
# the stored Healthy-Healthy experiment of the paper?
#
# Inputs (from the ENTEROTYPES mirror, never from the master data):
#   Healthy-Healthy_TotalEnvironment.RData   -> Remaining_Dataset (the 285 healthy
#                                               samples x 117 genera), Fixed_Dataset
#                                               (142 anchors), GenRichness
#   Healthy-Healthy_*                        -> read by read_legacy_erp()
#
# Output: a printed report and out/reproduce_check/report.rds
#
# Usage:  ENTEROTYPES_MIRROR=/path/to/ENTEROTYPES Rscript reproduce_check.R [n_cores]
# ------------------------------------------------------------------------------
suppressMessages({
  if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
    pkgload::load_all(".", quiet = TRUE)
  } else library(enterotyperobustness)
})
args <- commandArgs(trailingOnly = TRUE)
n_cores <- if (length(args)) as.integer(args[1]) else 10L
mirror <- Sys.getenv("ENTEROTYPES_MIRROR", "/Users/jdz/PROGRAMMATION/R/ENTEROTYPES")
exp_dir <- file.path(mirror, "analyses/2.Estimating_DMM/Article/Experiments")
out_dir <- file.path("out", "reproduce_check"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Legacy inputs -------------------------------------------------------
message("Reading the stored Healthy-Healthy experiment ...")
legacy <- read_legacy_erp(exp_dir, "Healthy-Healthy", keep_fits = TRUE)
e <- new.env()
load(file.path(exp_dir, "Healthy-Healthy_TotalEnvironment.RData"), envir = e)
abundance <- e$Remaining_Dataset                       # 285 healthy x 117 genera (non-integer)
anchors <- rownames(e$Fixed_Dataset)
complementary <- setdiff(rownames(abundance), anchors)  # the other 143, in table order
richness <- if (!is.null(e$GenRichness)) e$GenRichness else NULL
stopifnot(length(anchors) == 142, length(complementary) == 143, ncol(abundance) == 117)
message(sprintf("abundance %d x %d, %d anchors, %d complementary, richness %s",
                nrow(abundance), ncol(abundance), length(anchors), length(complementary),
                if (is.null(richness)) "absent" else "present"))

# ---- 2. Re-run with the legacy settings --------------------------------------
t0 <- Sys.time()
new <- erp(abundance, anchors, complementary, n_complementary = 86, n_iter = 500,
           marker_features = c(bacteroides = "Bacteroides", prevotella = "Prevotella"),
           diversity = richness, K = 4, iteration_seeds = 1:500, dmm_seed = 314,
           annotate = "rank", integerise = "trunc", n_cores = n_cores, keep_fits = TRUE)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

# ---- 3. Compare ---------------------------------------------------------------
same_index <- all(mapply(identical, new$index, legacy$index))
lab_new <- new$assignments[anchors, ]; lab_old <- legacy$assignments[anchors, ]
label_agreement <- mean(lab_new == lab_old)
per_iter_agreement <- colMeans(lab_new == lab_old)
map_same <- mean(apply(new$label_map == legacy$label_map, 1, all))
post_diff <- vapply(seq_len(500), function(i) {
  p_new <- new$posteriors[[i]]; p_old <- legacy$posteriors[[i]]
  max(abs(p_new[rownames(p_old), ] - p_old))
}, numeric(1))
alpha_diff <- vapply(seq_len(500), function(i) {
  a_new <- dmm_alpha(new$fits[[i]]); a_old <- dmm_alpha(legacy$fits[[i]])
  max(abs(a_new[, new$label_map[i, ]] - a_old[, legacy$label_map[i, ]]))
}, numeric(1))
flag_same <- c(bact2 = mean(new$flags$bact2_not_lowest_diversity == legacy$flags$bact2_not_lowest_diversity),
               rum = mean(new$flags$rum_not_highest_diversity == legacy$flags$rum_not_highest_diversity))

report <- list(
  elapsed_min = elapsed, n_cores = n_cores, sec_per_fit = elapsed * 60 * n_cores / 500,
  same_index = same_index, label_agreement = label_agreement,
  iterations_fully_identical = sum(per_iter_agreement == 1),
  worst_iteration = which.min(per_iter_agreement), worst_iteration_agreement = min(per_iter_agreement),
  label_map_identical = map_same,
  posterior_max_abs_diff = max(post_diff), posterior_median_max_diff = stats::median(post_diff),
  alpha_max_abs_diff = max(alpha_diff),
  flag_agreement = flag_same,
  n_flagged_new = colSums(new$flags[, -1]), n_flagged_legacy = colSums(legacy$flags[, 2:3]),
  DirichletMultinomial = as.character(utils::packageVersion("DirichletMultinomial")),
  R = R.version.string)
saveRDS(report, file.path(out_dir, "report.rds"))
cat("\n================ REPRODUCTION REPORT: Healthy-Healthy ================\n")
cat(sprintf("run time            : %.1f min on %d cores (%.1f s per fit)\n", elapsed, n_cores, report$sec_per_fit))
cat(sprintf("draws identical     : %s\n", same_index))
cat(sprintf("anchor labels equal : %.4f%%  (%d / 500 iterations fully identical; worst iteration %d at %.3f)\n",
            100 * label_agreement, report$iterations_fully_identical, report$worst_iteration, report$worst_iteration_agreement))
cat(sprintf("name->component map : identical in %.1f%% of iterations\n", 100 * map_same))
cat(sprintf("posterior |diff|    : max %.2e, median-of-max %.2e\n", max(post_diff), stats::median(post_diff)))
cat(sprintf("alpha |diff|        : max %.2e\n", max(alpha_diff)))
cat(sprintf("flags (bact2, rum)  : agreement %.3f, %.3f | new %s | legacy %s\n", flag_same[1], flag_same[2],
            paste(report$n_flagged_new[1:2], collapse = "/"), paste(report$n_flagged_legacy, collapse = "/")))
cat(sprintf("versions            : DirichletMultinomial %s, %s\n", report$DirichletMultinomial, report$R))
cat("=======================================================================\n")
