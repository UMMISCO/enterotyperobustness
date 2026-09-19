#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 10  Alternative numbers of DMM components (Reviewer 1, minor 4: "fixed K")
#     ERPM repeated with K = 3 and K = 5 (K = 4 from the stored runs). With K != 4
#     the rank naming does not apply, so components are matched to a reference
#     DMM fit of the same K on the full cohort (Hungarian assignment on the
#     Jensen-Shannon distance between centroids) and named C1..CK. For every
#     sample and fraction: robustness = frequency of its consensus component over
#     the draws in which it appeared; entropy normalised by log K. Reported: the
#     robustness curve along the fraction for each K, the fraction at which mean
#     robustness reaches 0.95, and whether ERP-core individuals (K = 4 lists) are
#     more robust than ERP-variable ones under K = 3 and K = 5.
#
# Inputs : METACARDIS_DATA, out/cache/erpm_legacy.rds (K = 4), mirror figure objects
# Outputs: out/10_alt_K/  robustness_by_K.csv, core_vs_variable_by_K.csv,
#          reference_crosstab.csv, Figure_alt_K.pdf, summary.txt ; out/cache/altK_*.rds
# Usage  : Rscript inst/scripts/revision/10_alt_K.R [n_cores]
# ------------------------------------------------------------------------------
suppressMessages({ pkgload::load_all(".", quiet = TRUE); library(ggplot2) })
source("inst/scripts/metacardis_data.R"); source("inst/scripts/cache_utils.R")
mirror <- Sys.getenv("ENTEROTYPES_MIRROR")
if (!nzchar(mirror) || !dir.exists(mirror)) stop("Set ENTEROTYPES_MIRROR to the directory holding the 2022 analyses (private, not distributed).")
args <- commandArgs(trailingOnly = TRUE); n_cores <- if (length(args)) as.integer(args[1]) else 12L
fig <- file.path(mirror, "analyses/9.analyses_mars2026")
out <- "out/10_alt_K"; dir.create(out, showWarnings = FALSE, recursive = TRUE); dir.create("out/cache", showWarnings = FALSE)
ent <- enterotype_names()
Ks <- c(3, 5); fractions <- c(0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.70, 0.90); n_iter <- 500
md <- load_metacardis_data(); A <- md$abundance; groups <- md$groups[rownames(A)]

# ---- 1. Reference fits on the full cohort, K = 3, 4, 5 ----------------------------------------
ref_file <- "out/cache/altK_reference_fits.rds"
inputs <- list(ids = rownames(A), groups = groups, dims = dim(A), total = sum(A))
ref_manifest <- cache_manifest(params = list(K = c(3, 4, 5), seed = 314, integerise = "trunc"), inputs = inputs)
refs <- cache_read(ref_file, ref_manifest)
if (is.null(refs)) {
  refs <- lapply(c(3, 4, 5), function(K) { t1 <- Sys.time(); f <- fit_dmm(A, K = K, seed = 314, integerise = "trunc")
    cat(sprintf("reference fit K = %d: %.1f min\n", K, as.numeric(difftime(Sys.time(), t1, units = "mins")))); f })
  names(refs) <- paste0("K", c(3, 4, 5)); cache_write(refs, ref_file, ref_manifest)
}
ann4 <- annotate_clusters(refs$K4, A, md$marker_features, md$richness)
assign4 <- ann4$component_names[dmm_assign(refs$K4)]
crosstab <- do.call(rbind, lapply(Ks, function(K) { t <- table(K4 = assign4, alt = paste0("C", dmm_assign(refs[[paste0("K", K)]])))
  d <- as.data.frame(t); d$K <- K; d }))
write.csv(crosstab, file.path(out, "reference_crosstab.csv"), row.names = FALSE)
ref_cent <- lapply(Ks, function(K) { cent <- dmm_centroids(refs[[paste0("K", K)]]); colnames(cent) <- paste0("C", seq_len(K)); cent }); names(ref_cent) <- paste0("K", Ks)

# ---- 2. ERPM with K = 3 and 5, components matched to the reference ----------------------------
one_draw <- function(i, Q, K) {
  ids <- stratified_subsample(groups, Q, seed = i)
  fit <- fit_dmm(A[ids, ], K = K, seed = 314, integerise = "trunc")
  cent <- dmm_centroids(fit); rc <- ref_cent[[paste0("K", K)]]
  cost <- jsd_matrix(cent, rc); sol <- as.integer(clue::solve_LSAP(cost))
  assigned <- cost[cbind(seq_len(K), sol)]; second <- vapply(seq_len(K), function(k) min(cost[k, -sol[k]]), numeric(1))
  list(ids = ids, labels = colnames(rc)[sol][dmm_assign(fit)], margin = min(second - assigned))
}
alt <- list()
for (K in Ks) {
  cache_file <- sprintf("out/cache/altK_K%d.rds", K)
  k_manifest <- cache_manifest(params = list(K = K, fractions = fractions, n_iter = n_iter, dmm_seed = 314, integerise = "trunc", matching = "Hungarian on JSD to the full-cohort reference"), inputs = inputs)
  cached <- cache_read(cache_file, k_manifest)
  if (!is.null(cached)) { alt[[paste0("K", K)]] <- cached; next }
  res <- list()
  for (Q in fractions) {
    t1 <- Sys.time()
    runs <- parallel::mclapply(seq_len(n_iter), one_draw, Q = Q, K = K, mc.cores = n_cores, mc.preschedule = FALSE, mc.set.seed = FALSE)
    failed <- vapply(runs, inherits, logical(1), "try-error"); if (any(failed)) stop("K = ", K, " Q = ", Q, ": ", sum(failed), " fits failed: ", as.character(runs[[which(failed)[1]]]))
    m <- matrix(NA_character_, nrow(A), n_iter, dimnames = list(rownames(A), paste0("iter", seq_len(n_iter))))
    for (j in seq_len(n_iter)) m[runs[[j]]$ids, j] <- runs[[j]]$labels
    res[[sprintf("Q%.2f", Q)]] <- list(fraction = Q, assignments = m, margin = vapply(runs, `[[`, numeric(1), "margin"))
    cat(sprintf("K = %d, Q = %.2f: %d fits of %d samples, matching margin median %.3f (%.1f min)\n", K, Q, n_iter, length(runs[[1]]$ids),
                stats::median(res[[sprintf("Q%.2f", Q)]]$margin), as.numeric(difftime(Sys.time(), t1, units = "mins"))))
  }
  alt[[paste0("K", K)]] <- res; cache_write(res, cache_file, k_manifest)
}

# ---- 3. Robustness per sample, fraction and K (K = 4 from the stored runs) --------------------
freq_generic <- function(m, levels) { nd <- rowSums(!is.na(m)); f <- vapply(levels, function(k) rowSums(m == k, na.rm = TRUE) / nd, numeric(nrow(m)))
  f <- matrix(f, nrow(m), length(levels), dimnames = list(rownames(m), levels)); f[nd == 0, ] <- NA; attr(f, "n_drawn") <- nd; f }
rob_of <- function(m, levels, K) {
  f <- freq_generic(m, levels); am <- argmax_rows(f, seed = 314); ok <- !is.na(am$index)
  data.frame(sample = rownames(f)[ok], label = levels[am$index[ok]], robustness = f[cbind(seq_len(nrow(f)), am$index)][ok],
             entropy = assignment_entropy(f[ok, , drop = FALSE]), n_drawn = attr(f, "n_drawn")[ok], stringsAsFactors = FALSE)
}
cache <- readRDS("out/cache/erpm_legacy.rds")
rob <- do.call(rbind, c(
  lapply(fractions, function(Q) { key <- cache$fractions$file[abs(cache$fractions$fraction - Q) < 1e-9]
    cbind(K = 4, fraction = Q, rob_of(cache$assignments[[key]], 1:4, 4)) }),
  unlist(lapply(Ks, function(K) lapply(alt[[paste0("K", K)]], function(r) cbind(K = K, fraction = r$fraction, rob_of(r$assignments, paste0("C", seq_len(K)), K)))), recursive = FALSE)))
rob$label <- ifelse(rob$K == 4, ent[as.integer(rob$label)], rob$label)
write.csv(rob, file.path(out, "robustness_by_K.csv"), row.names = FALSE)
curve <- aggregate(cbind(robustness, entropy) ~ K + fraction, rob, mean)
q95 <- aggregate(fraction ~ K, curve[curve$robustness >= 0.95, ], min)

# ---- 4. ERP-core vs ERP-variable (published K = 4 lists) under each K ---------------------------
st <- load_published_core_variable(fig)$table[, c("sample", "enterotype", "status")]
core_ids <- unique(st$sample[st$status == "core"]); var_ids <- setdiff(unique(st$sample[st$status == "variable"]), core_ids)
cv <- do.call(rbind, lapply(split(rob, list(rob$K, rob$fraction)), function(x) {
  rc <- x$robustness[x$sample %in% core_ids]; rv <- x$robustness[x$sample %in% var_ids]
  data.frame(K = x$K[1], fraction = x$fraction[1], n_core = length(rc), n_variable = length(rv), mean_core = mean(rc), mean_variable = mean(rv),
             difference = mean(rc) - mean(rv), p_wilcoxon = stats::wilcox.test(rc, rv, exact = FALSE)$p.value)
}))
cv <- cv[order(cv$K, cv$fraction), ]
cv$p_adj <- stats::ave(cv$p_wilcoxon, cv$K, FUN = function(v) stats::p.adjust(v, "BH"))    # one family per K (the fractions)
cv$core_higher_adj <- cv$p_adj < 0.05 & cv$difference > 0
write.csv(cv, file.path(out, "core_vs_variable_by_K.csv"), row.names = FALSE)
# curve-level test, one per K: each sample's robustness averaged over the fractions
curve_test <- do.call(rbind, lapply(split(rob, rob$K), function(x) {
  per_sample <- tapply(x$robustness, x$sample, mean)
  rc <- per_sample[names(per_sample) %in% core_ids]; rv <- per_sample[names(per_sample) %in% var_ids]
  data.frame(K = x$K[1], n_core = length(rc), n_variable = length(rv), mean_core = mean(rc), mean_variable = mean(rv),
             difference = mean(rc) - mean(rv), p_wilcoxon = stats::wilcox.test(rc, rv, exact = FALSE)$p.value)
}))
curve_test$p_adj <- stats::p.adjust(curve_test$p_wilcoxon, "BH")
write.csv(curve_test, file.path(out, "core_vs_variable_curve_level.csv"), row.names = FALSE)

# ---- 5. Figure and summary -------------------------------------------------------------------
long <- rbind(data.frame(cv[, c("K", "fraction")], status = "core", mean = cv$mean_core),
              data.frame(cv[, c("K", "fraction")], status = "variable", mean = cv$mean_variable),
              data.frame(curve[, c("K", "fraction")], status = "all samples", mean = curve$robustness))
long$K <- factor(paste("K =", long$K))
p <- ggplot(long, aes(fraction, mean, colour = status)) + geom_line() + geom_point(size = 0.8) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey40") + facet_wrap(~K, nrow = 1) +
  scale_colour_manual(values = c(core = "#1b9e77", variable = "#d95f02", `all samples` = "grey30")) +
  labs(x = "fraction of the cohort", y = "mean robustness", title = "ERPM robustness with 3, 4 and 5 DMM components (core / variable status from the K = 4 ERP)") + theme_bw()
ggsave(file.path(out, "Figure_alt_K.pdf"), p, width = 12, height = 4.5)
sink(file.path(out, "summary.txt"))
cat("10  Alternative K, run", format(Sys.time()), "\n\n")
cat("Reference fits on the full cohort: K = 4 enterotypes vs K = 3 / 5 components (samples):\n")
print(stats::reshape(crosstab, idvar = c("K", "K4"), timevar = "alt", direction = "wide"), row.names = FALSE)
cat("\nMatching margin (JSD, second best minus assigned) per K and fraction, median:\n")
for (K in Ks) cat(sprintf("  K = %d: %s\n", K, paste(vapply(alt[[paste0("K", K)]], function(r) sprintf("%.2f:%.3f", r$fraction, stats::median(r$margin)), character(1)), collapse = " ")))
cat("\nMean robustness and normalised entropy along the fraction:\n"); print(stats::reshape(curve, idvar = "fraction", timevar = "K", direction = "wide"), row.names = FALSE, digits = 3)
cat("\nFirst fraction at which mean robustness reaches 0.95:\n"); print(q95, row.names = FALSE)
cat("\nERP-core vs ERP-variable under each K, curve level (one test per K, BH across the three):\n"); print(curve_test, row.names = FALSE, digits = 3)
cat("\nPer fraction (BH within K; core_higher_adj = difference > 0 and adjusted p < 0.05):\n"); print(cv, row.names = FALSE, digits = 3)
cat("\nSession:\n"); print(sessionInfo())
sink()
cat(readLines(file.path(out, "summary.txt"), n = 70), sep = "\n")
