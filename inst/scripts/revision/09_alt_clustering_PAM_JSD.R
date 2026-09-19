#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 09  Independent clustering: PAM on Jensen-Shannon distance (Reviewer 1, major 2)
#     The core / variable classification comes from ERP (DMM). Is it circular, i.e.
#     an artefact of DMM? Here the ERPM size-perturbation is repeated with a
#     clustering that shares nothing with DMM: partitioning around medoids on the
#     Jensen-Shannon distance between relative-abundance profiles (the Arumugam
#     2011 / Costea 2018 protocol, as in the enterotyper code of the mirror), K = 4,
#     clusters named by the same rank heuristic as the paper. For every sample and
#     fraction: PAM robustness (frequency of its consensus PAM label over the draws
#     in which it appeared). Tested: are ERP-core individuals more robust than
#     ERP-variable ones under PAM (Wilcoxon per enterotype and fraction, as in
#     Figure 3B), and how often do the PAM and DMM consensus labels agree.
#
# Inputs : METACARDIS_DATA (master table, clinical groups), out/cache/erpm_legacy.rds
#          (DMM ERPM, script 00), mirror figure objects (published core/variable lists)
# Outputs: out/09_pam_jsd/  pam_robustness.csv, core_vs_variable.csv, agreement_with_dmm.csv,
#          silhouette_full_cohort.csv, Figure_PAM_core_variable.pdf, summary.txt ;
#          out/cache/pam_assignments.rds
# Usage  : Rscript inst/scripts/revision/09_alt_clustering_PAM_JSD.R [n_cores]
# ------------------------------------------------------------------------------
suppressMessages({ pkgload::load_all(".", quiet = TRUE); library(ggplot2); library(cluster) })
source("inst/scripts/metacardis_data.R"); source("inst/scripts/cache_utils.R")
mirror <- Sys.getenv("ENTEROTYPES_MIRROR")
if (!nzchar(mirror) || !dir.exists(mirror)) stop("Set ENTEROTYPES_MIRROR to the directory holding the 2022 analyses (private, not distributed).")
args <- commandArgs(trailingOnly = TRUE); n_cores <- if (length(args)) as.integer(args[1]) else 12L
fig <- file.path(mirror, "analyses/9.analyses_mars2026")
out <- "out/09_pam_jsd"; dir.create(out, showWarnings = FALSE, recursive = TRUE); dir.create("out/cache", showWarnings = FALSE)
ent <- enterotype_names()
fractions <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.75, 0.90); n_iter <- 500; K <- 4
md <- load_metacardis_data(); A <- md$abundance; groups <- md$groups[rownames(A)]
stopifnot(!anyNA(groups))

# ---- 1. Jensen-Shannon distance between all samples (once) ------------------------------------
# sqrt(JSD) on relative abundances: zeros replaced by a pseudocount and each profile
# renormalised to sum 1 (the enterotyper's getDist.R adds the pseudocount without renormalising)
jsd_all <- function(P, pseudocount = 1e-6) {
  X <- t(P); X[X == 0] <- pseudocount                       # features x samples
  X <- sweep(X, 2, colSums(X), "/")
  KLD <- sapply(as.data.frame(X), function(x) crossprod(x, log(2 * x / (x + X))))
  D <- sqrt(0.5 * (KLD + t(KLD))); dimnames(D) <- list(colnames(X), colnames(X)); D
}
t0 <- Sys.time()
D <- jsd_all(A / rowSums(A))
cat(sprintf("JSD matrix %d x %d in %.1f s\n", nrow(D), ncol(D), as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# ---- 2. PAM + rank naming --------------------------------------------------------------------
name_by_rank <- function(cl, ids) {                          # Bact2 = highest mean Bacteroides, Bact1 = second, Prevotella = highest Prevotella of the rest
  bact <- tapply(A[ids, md$marker_features[["bacteroides"]]], cl, mean)
  prev <- tapply(A[ids, md$marker_features[["prevotella"]]], cl, mean)
  o <- order(bact, decreasing = TRUE); nm <- character(K); names(nm) <- names(bact)
  nm[o[1]] <- "Bact2"; nm[o[2]] <- "Bact1"
  rest <- o[3:4]; collision <- which.max(prev) %in% o[1:2]
  nm[rest[which.max(prev[rest])]] <- "Prevotella"; nm[nm == ""] <- "Ruminococcus"
  list(names = nm, collision = collision)
}
one_draw <- function(i, Q) {
  ids <- stratified_subsample(groups, Q, seed = i)
  cl <- cluster::pam(stats::as.dist(D[ids, ids]), k = K, diss = TRUE, cluster.only = TRUE, pamonce = 5)
  nb <- name_by_rank(cl, ids)
  list(ids = ids, labels = unname(nb$names[as.character(cl)]), collision = nb$collision)
}
cache_file <- "out/cache/pam_assignments.rds"
manifest <- cache_manifest(params = list(fractions = fractions, n_iter = n_iter, K = K, distance = "sqrt JSD, pseudocount 1e-6, renormalised", naming = "rank", pamonce = 5),
                           inputs = list(A = A, groups = groups))
pam <- cache_read(cache_file, manifest)
if (is.null(pam)) {
  pam <- list()
  for (Q in fractions) {
    t1 <- Sys.time()
    runs <- parallel::mclapply(seq_len(n_iter), one_draw, Q = Q, mc.cores = n_cores, mc.preschedule = TRUE)
    failed <- vapply(runs, inherits, logical(1), "try-error")
    if (any(failed)) stop("Q = ", Q, ": ", sum(failed), " PAM draws failed. First error: ", as.character(runs[[which(failed)[1]]]))
    m <- matrix(NA_character_, nrow(A), n_iter, dimnames = list(rownames(A), paste0("iter", seq_len(n_iter))))
    for (j in seq_len(n_iter)) m[runs[[j]]$ids, j] <- runs[[j]]$labels
    pam[[sprintf("Q%.2f", Q)]] <- list(fraction = Q, assignments = m, collisions = sum(vapply(runs, `[[`, logical(1), "collision")))
    cat(sprintf("Q = %.2f: %d draws of %d samples, %d Prevotella collisions (%.1f min)\n", Q, n_iter, length(runs[[1]]$ids),
                pam[[sprintf("Q%.2f", Q)]]$collisions, as.numeric(difftime(Sys.time(), t1, units = "mins"))))
  }
  cache_write(pam, cache_file, manifest)
}

# ---- 3. Robustness per sample and fraction; DMM counterpart from the cache ----------------------
cache <- load_erpm_cache(mirror)
freq_int <- function(m) { nd <- rowSums(!is.na(m)); f <- vapply(1:4, function(k) rowSums(m == k, na.rm = TRUE) / nd, numeric(nrow(m)))
  dimnames(f) <- list(rownames(m), ent); f[nd == 0, ] <- NA; attr(f, "n_drawn") <- nd; f }
rob <- do.call(rbind, lapply(pam, function(p) {
  r <- robustness(assignment_frequency(p$assignments)); r <- r[!is.na(r$enterotype), ]
  key <- cache$fractions$file[abs(cache$fractions$fraction - p$fraction) < 1e-9]
  d <- robustness(freq_int(cache$assignments[[key]])); d <- d[!is.na(d$enterotype), ]
  m <- match(r$sample, d$sample)
  data.frame(fraction = p$fraction, sample = r$sample, pam_enterotype = r$enterotype, pam_robustness = r$robustness, pam_n_drawn = r$n_drawn,
             dmm_enterotype = d$enterotype[m], dmm_robustness = d$robustness[m], stringsAsFactors = FALSE)
}))
write.csv(rob, file.path(out, "pam_robustness.csv"), row.names = FALSE)
agree <- do.call(rbind, lapply(split(rob, rob$fraction), function(x) data.frame(
  fraction = x$fraction[1], n = nrow(x), same_label = mean(x$pam_enterotype == x$dmm_enterotype, na.rm = TRUE),
  ari = adjusted_rand_index(x$pam_enterotype, x$dmm_enterotype),
  spearman_robustness = suppressWarnings(stats::cor(x$pam_robustness, x$dmm_robustness, method = "spearman", use = "complete.obs")),
  mean_pam = mean(x$pam_robustness), mean_dmm = mean(x$dmm_robustness, na.rm = TRUE))))
write.csv(agree, file.path(out, "agreement_with_dmm.csv"), row.names = FALSE)

# ---- 4. ERP-core vs ERP-variable under PAM (Figure 3B analogue) --------------------------------
st <- load_published_core_variable(fig)$table[, c("sample", "enterotype", "status")]
st <- st[!duplicated(st[, c("sample", "enterotype")]), ]
cv <- do.call(rbind, lapply(split(rob, rob$fraction), function(x) do.call(rbind, lapply(ent, function(e) {
  ids_c <- st$sample[st$enterotype == e & st$status == "core"]; ids_v <- st$sample[st$enterotype == e & st$status == "variable"]
  for (clus in c("pam", "dmm")) {
    lab <- x[[paste0(clus, "_enterotype")]]; rb <- x[[paste0(clus, "_robustness")]]
    rc <- rb[lab == e & x$sample %in% ids_c]; rv <- rb[lab == e & x$sample %in% ids_v]
    rc <- rc[!is.na(rc)]; rv <- rv[!is.na(rv)]
    p <- if (length(rc) >= 3 && length(rv) >= 3) stats::wilcox.test(rc, rv, exact = FALSE)$p.value else NA
    row <- data.frame(clustering = clus, fraction = x$fraction[1], enterotype = e, n_core = length(rc), n_variable = length(rv),
                      mean_core = mean(rc), mean_variable = mean(rv), difference = mean(rc) - mean(rv), p_wilcoxon = p)
    if (clus == "pam") out_rows <- row else out_rows <- rbind(out_rows, row)
  }
  out_rows
}))))
# one family of tests per clustering (fractions x enterotypes): BH within family, and the direction is checked
cv$p_adj <- stats::ave(cv$p_wilcoxon, cv$clustering, FUN = function(v) { a <- v; ok <- !is.na(v); a[ok] <- stats::p.adjust(v[ok], "BH"); a })
cv$core_higher_adj <- !is.na(cv$p_adj) & cv$p_adj < 0.05 & cv$difference > 0
write.csv(cv, file.path(out, "core_vs_variable.csv"), row.names = FALSE)
# curve-level test, one per clustering and enterotype: each sample's robustness averaged over the fractions, then core vs variable
curve_test <- do.call(rbind, lapply(c("pam", "dmm"), function(clus) do.call(rbind, lapply(ent, function(e) {
  ids_c <- st$sample[st$enterotype == e & st$status == "core"]; ids_v <- st$sample[st$enterotype == e & st$status == "variable"]
  lab <- rob[[paste0(clus, "_enterotype")]]; rb <- rob[[paste0(clus, "_robustness")]]
  x <- rob[!is.na(lab) & lab == e & !is.na(rb), ]; rbv <- rb[!is.na(lab) & lab == e & !is.na(rb)]
  per_sample <- tapply(rbv, x$sample, mean)
  rc <- per_sample[names(per_sample) %in% ids_c]; rv <- per_sample[names(per_sample) %in% ids_v]
  data.frame(clustering = clus, enterotype = e, n_core = length(rc), n_variable = length(rv), mean_core = mean(rc), mean_variable = mean(rv),
             difference = mean(rc) - mean(rv), p_wilcoxon = if (length(rc) >= 3 && length(rv) >= 3) stats::wilcox.test(rc, rv, exact = FALSE)$p.value else NA)
}))))
curve_test$p_adj <- stats::p.adjust(curve_test$p_wilcoxon, "BH")
write.csv(curve_test, file.path(out, "core_vs_variable_curve_level.csv"), row.names = FALSE)
# status taken from ERP under PAM regardless of the PAM label (a sample may be named differently by PAM)
cv_any <- do.call(rbind, lapply(split(rob, rob$fraction), function(x) {
  core_ids <- unique(st$sample[st$status == "core"]); var_ids <- setdiff(unique(st$sample[st$status == "variable"]), core_ids)
  rc <- x$pam_robustness[x$sample %in% core_ids]; rv <- x$pam_robustness[x$sample %in% var_ids]
  data.frame(fraction = x$fraction[1], n_core = length(rc), n_variable = length(rv), mean_core = mean(rc), mean_variable = mean(rv),
             difference = mean(rc) - mean(rv), p_wilcoxon = stats::wilcox.test(rc, rv, exact = FALSE)$p.value)
}))
cv_any$p_adj <- stats::p.adjust(cv_any$p_wilcoxon, "BH"); cv_any$core_higher_adj <- cv_any$p_adj < 0.05 & cv_any$difference > 0

# ---- 5. Silhouette of PAM on the full cohort, context ----------------------------------------
sil <- do.call(rbind, lapply(2:6, function(k) { cl <- cluster::pam(stats::as.dist(D), k = k, diss = TRUE, cluster.only = TRUE, pamonce = 5)
  data.frame(k = k, mean_silhouette = mean(cluster::silhouette(cl, stats::as.dist(D))[, 3])) }))
write.csv(sil, file.path(out, "silhouette_full_cohort.csv"), row.names = FALSE)

# ---- 6. Figure and summary -------------------------------------------------------------------
long <- rbind(data.frame(cv[, c("clustering", "fraction", "enterotype")], status = "core", mean = cv$mean_core),
              data.frame(cv[, c("clustering", "fraction", "enterotype")], status = "variable", mean = cv$mean_variable))
long$clustering <- factor(long$clustering, levels = c("dmm", "pam"), labels = c("DMM (published)", "PAM on Jensen-Shannon"))
p <- ggplot(long, aes(fraction, mean, colour = status)) + geom_line() + geom_point(size = 0.8) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey40") +
  facet_grid(enterotype ~ clustering) + scale_colour_manual(values = c(core = "#1b9e77", variable = "#d95f02")) +
  labs(x = "fraction of the cohort", y = "mean robustness", title = "ERP-core vs ERP-variable individuals under DMM and under an independent clustering (PAM / JSD)") +
  theme_bw()
ggsave(file.path(out, "Figure_PAM_core_variable.pdf"), p, width = 9, height = 9)
sink(file.path(out, "summary.txt"))
cat("09  PAM on Jensen-Shannon, run", format(Sys.time()), "\n")
cat(sprintf("%d fractions x %d draws, K = %d, rank naming; Prevotella collisions per fraction: %s\n\n", length(fractions), n_iter, K,
            paste(vapply(pam, function(p) sprintf("%.2f:%d", p$fraction, p$collisions), character(1)), collapse = " ")))
cat("Silhouette of PAM on the full cohort:\n"); print(sil, row.names = FALSE, digits = 3)
cat("\nAgreement between PAM and DMM consensus labels and robustness, per fraction:\n"); print(agree, row.names = FALSE, digits = 3)
cat("\nERP-core vs ERP-variable robustness under PAM, all enterotypes pooled (status from ERP, label-free):\n"); print(cv_any, row.names = FALSE, digits = 3)
cat("\nCurve-level test (one per clustering and enterotype; each sample's robustness averaged over the fractions; BH across the 8 tests):\n")
print(curve_test, row.names = FALSE, digits = 3)
cat("\nPer enterotype, fractions at which core > variable (difference > 0 and BH-adjusted p < 0.05 within clustering), PAM vs DMM:\n")
sig <- aggregate(cbind(core_higher_adj, tested = !is.na(p_wilcoxon)) ~ clustering + enterotype, cv, sum, na.action = stats::na.pass)
sig$count <- sprintf("%d/%d", sig$core_higher_adj, sig$tested); sig <- sig[, c("clustering", "enterotype", "count")]
print(stats::reshape(sig, idvar = "clustering", timevar = "enterotype", direction = "wide"), row.names = FALSE)
cat("\nMean difference core - variable at Q = 0.25 / 0.50 / 0.75:\n")
print(cv[cv$fraction %in% c(0.25, 0.5, 0.75), ], row.names = FALSE, digits = 3)
cat("\nSession:\n"); print(sessionInfo())
sink()
cat(readLines(file.path(out, "summary.txt"), n = 60), sep = "\n")
