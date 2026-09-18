#!/usr/bin/env Rscript
# Third full run of the Healthy-Healthy reproduction, compared iteration by
# iteration with the saved objects of the previous run (out/reproduce_check/objects.rds).
# Answers: is erp() bit-stable from one run to the next on this machine?
suppressMessages(pkgload::load_all(".", quiet = TRUE))
args <- commandArgs(trailingOnly = TRUE); n_cores <- if (length(args)) as.integer(args[1]) else 10L
obj <- readRDS("out/reproduce_check/objects.rds"); prev <- obj$new
mirror <- Sys.getenv("ENTEROTYPES_MIRROR", "/Users/jdz/PROGRAMMATION/R/ENTEROTYPES")
e <- new.env(); load(file.path(mirror, "analyses/2.Estimating_DMM/Article/Experiments/Healthy-Healthy_TotalEnvironment.RData"), envir = e)
abundance <- e$Remaining_Dataset
t0 <- Sys.time()
new <- erp(abundance, prev$samples$anchors, prev$samples$complementary, n_complementary = 86, n_iter = 500,
           marker_features = c(bacteroides = "Bacteroides", prevotella = "Prevotella"),
           diversity = e$GenRichness, iteration_seeds = 1:500, dmm_seed = 314,
           annotate = "rank", integerise = "trunc", n_cores = n_cores, keep_fits = FALSE, verbose = FALSE)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
same_draw <- all(mapply(identical, new$index, prev$index))
lab_same <- colMeans(new$assignments == prev$assignments, na.rm = TRUE)
post_diff <- vapply(1:500, function(i) max(abs(new$posteriors[[i]] - prev$posteriors[[i]][rownames(new$posteriors[[i]]), ])), numeric(1))
cat("\n================ RUN-TO-RUN STABILITY (run 3 vs run 2) ================\n")
cat(sprintf("run time                 : %.1f min on %d cores\n", elapsed, n_cores))
cat(sprintf("draws identical          : %s\n", same_draw))
cat(sprintf("iterations with identical labels : %d / 500\n", sum(lab_same == 1)))
cat(sprintf("posterior max |diff|     : %.2e (iterations > 1e-8: %d)\n", max(post_diff), sum(post_diff > 1e-8)))
if (any(post_diff > 1e-8)) cat("differing iterations     :", which(post_diff > 1e-8), "\n")
cat("=========================================================================\n")
saveRDS(list(lab_same = lab_same, post_diff = post_diff), "out/reproduce_check/run3_vs_run2.rds")
