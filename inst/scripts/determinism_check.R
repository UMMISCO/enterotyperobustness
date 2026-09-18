#!/usr/bin/env Rscript
# Is a DMM fit reproducible run to run? Refit a subset of Healthy-Healthy
# iterations twice sequentially and twice with forked workers, and compare.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
obj <- readRDS("out/reproduce_check/objects.rds")
rep <- readRDS("out/reproduce_check/report.rds")
new <- obj$new
mirror <- Sys.getenv("ENTEROTYPES_MIRROR", "/Users/jdz/PROGRAMMATION/R/ENTEROTYPES")
e <- new.env(); load(file.path(mirror, "analyses/2.Estimating_DMM/Article/Experiments/Healthy-Healthy_TotalEnvironment.RData"), envir = e)
abundance <- e$Remaining_Dataset; anchors <- new$samples$anchors; comp <- new$samples$complementary
iters <- c(head(rep$nle_cmp$iteration, 12), 1:8)          # 12 differing + 8 identical iterations
id_list <- lapply(new$index[iters], function(ix) c(anchors, comp[ix]))
markers <- c(bacteroides = "Bacteroides", prevotella = "Prevotella")
run <- function(cores) run_dmm_iterations(abundance, id_list, 4, 314, "trunc", markers, NULL,
                                          "rank", NULL, NULL, cores, TRUE, FALSE)
cmp <- function(a, b) vapply(seq_along(a), function(i) max(abs(dmm_posterior(a[[i]]$fit) - dmm_posterior(b[[i]]$fit))), numeric(1))
cat("BLAS:", extSoftVersion()[["BLAS"]], "\n")
s1 <- run(1); s2 <- run(1); p1 <- run(10); p2 <- run(10)
d_ss <- cmp(s1, s2); d_pp <- cmp(p1, p2); d_sp <- cmp(s1, p1)
d_ref <- vapply(seq_along(iters), function(k) max(abs(dmm_posterior(s1[[k]]$fit) - new$posteriors[[iters[k]]][rownames(dmm_posterior(s1[[k]]$fit)), new$label_map[iters[k], ]])), numeric(1))
cat(sprintf("iterations tested      : %s\n", paste(iters, collapse = " ")))
cat(sprintf("seq vs seq   max|dpost|: %.2e (n differing > 1e-8: %d)\n", max(d_ss), sum(d_ss > 1e-8)))
cat(sprintf("par vs par   max|dpost|: %.2e (n differing > 1e-8: %d)\n", max(d_pp), sum(d_pp > 1e-8)))
cat(sprintf("seq vs par   max|dpost|: %.2e (n differing > 1e-8: %d)\n", max(d_sp), sum(d_sp > 1e-8)))
cat(sprintf("seq vs 10-core run     : %.2e (n differing > 1e-8: %d)\n", max(d_ref), sum(d_ref > 1e-8)))
