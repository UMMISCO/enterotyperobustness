#!/usr/bin/env Rscript
# Run-to-run stability of erp() at full scale, after switching to one fork per
# iteration. Reference = sequential refits of the iterations that differed
# between two earlier prescheduled parallel runs; then two full parallel runs.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
prev <- readRDS("out/reproduce_check/objects.rds")$new
r32 <- readRDS("out/reproduce_check/run3_vs_run2.rds")
bad <- which(r32$post_diff > 1e-8)
mirror <- Sys.getenv("ENTEROTYPES_MIRROR", "/Users/jdz/PROGRAMMATION/R/ENTEROTYPES")
e <- new.env(); load(file.path(mirror, "analyses/2.Estimating_DMM/Article/Experiments/Healthy-Healthy_TotalEnvironment.RData"), envir = e)
A <- e$Remaining_Dataset; anchors <- prev$samples$anchors; comp <- prev$samples$complementary
mk <- c(bacteroides = "Bacteroides", prevotella = "Prevotella")
run <- function(iters, cores) {
  erp(A, anchors, comp, n_complementary = 86, n_iter = 500, marker_features = mk, diversity = e$GenRichness,
      iteration_seeds = 1:500, dmm_seed = 314, n_cores = cores, keep_fits = FALSE, verbose = FALSE)
}
# sequential reference on the differing iterations only: reuse erp() machinery via id_list
id_list <- lapply(prev$index[bad], function(ix) c(anchors, comp[ix]))
seq1 <- run_dmm_iterations(A, id_list, 4, 314, "trunc", mk, e$GenRichness, "rank", NULL, NULL, 1, FALSE, FALSE)
seq2 <- run_dmm_iterations(A, id_list, 4, 314, "trunc", mk, e$GenRichness, "rank", NULL, NULL, 1, FALSE, FALSE)
dseq <- vapply(seq_along(bad), function(k) max(abs(seq1[[k]]$posterior - seq2[[k]]$posterior)), numeric(1))
cat(sprintf("sequential x2 on the %d unstable iterations: max|dpost| = %.2e\n", length(bad), max(dseq)))
tA <- system.time(runA <- run(NULL, 10))[["elapsed"]]
runB <- run(NULL, 10)
cmp_post <- function(x, y) vapply(1:500, function(i) max(abs(x$posteriors[[i]] - y$posteriors[[i]][rownames(x$posteriors[[i]]), ])), numeric(1))
dAB <- cmp_post(runA, runB); dA2 <- cmp_post(runA, prev)
dAseq <- vapply(seq_along(bad), function(k) max(abs(runA$posteriors[[bad[k]]] - seq1[[k]]$posterior[rownames(runA$posteriors[[bad[k]]]), ])), numeric(1))
cat("\n=========== DETERMINISM WITH ONE FORK PER ITERATION ===========\n")
cat(sprintf("full run A time            : %.1f min on 10 cores\n", tA / 60))
cat(sprintf("run A vs run B (parallel)  : max|dpost| %.2e, iterations differing %d / 500\n", max(dAB), sum(dAB > 1e-8)))
cat(sprintf("run A vs sequential ref    : max|dpost| %.2e on the %d formerly unstable iterations\n", max(dAseq), length(bad)))
cat(sprintf("run A vs prescheduled run 2: iterations differing %d / 500\n", sum(dA2 > 1e-8)))
cat(sprintf("labels A vs B identical    : %d / 500 iterations\n", sum(colMeans(runA$assignments == runB$assignments, na.rm = TRUE) == 1)))
cat("=================================================================\n")
saveRDS(list(dAB = dAB, dA2 = dA2, dAseq = dAseq, runA_assign = runA$assignments), "out/reproduce_check/determinism_full.rds")
