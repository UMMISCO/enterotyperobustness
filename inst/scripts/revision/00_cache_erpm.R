#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 00  Read the 56 legacy ERPM environments once and cache what the revision
#     scripts need (02, 05, 06, 07). Each Conf*_TotalEnvironment.RData is 20 to
#     160 MB; reading them all takes about ten minutes and is never repeated.
#
# Cache: out/cache/erpm_legacy.rds, a list with
#   fractions   : data.frame(file, fraction, n_samples, n_iter)
#   assignments : list per fraction of integer matrix samples x iterations
#                 (1 = Bact1, 2 = Bact2, 3 = Prevotella, 4 = Ruminococcus, NA = not drawn)
#   prob_enterotype : the stored per-enterotype robustness lists (what Fig. 3/4 plotted)
#   flags       : list per fraction of the annotation flags (stored WARNING + re-derived)
#   gap         : list per fraction of the Bact1/Bact2 gap table (re-derived)
#   ties        : W_50 per fraction
#   samples     : the 2 022 sample IDs
# ------------------------------------------------------------------------------
suppressMessages(pkgload::load_all(".", quiet = TRUE))
mirror <- Sys.getenv("ENTEROTYPES_MIRROR", "/Users/jdz/PROGRAMMATION/R/ENTEROTYPES")
conf <- file.path(mirror, "analyses/2.Estimating_DMM/Article/Confidence")
dir.create("out/cache", showWarnings = FALSE, recursive = TRUE)
files <- sort(list.files(conf, pattern = "_TotalEnvironment\\.RData$", full.names = TRUE))
# the hand-maintained reference table maps file -> fraction (file names are irregular)
ref <- read.table(file.path(conf, "TotalEnvironmentTable_21052024.txt"), sep = "\t", header = FALSE,
                  quote = "", comment.char = "", stringsAsFactors = FALSE)
frac_of <- setNames(as.numeric(ref$V2) / 10, basename(ref$V1))
cat(length(files), "environments;", sum(basename(files) %in% names(frac_of)), "mapped in the reference table\n")
code <- setNames(1:4, enterotype_names())
cache <- list(fractions = NULL, assignments = list(), prob_enterotype = list(), flags = list(), gap = list(), ties = list())
t0 <- Sys.time()
for (f in files) {
  key <- basename(f)
  r <- tryCatch(read_legacy_erpm(f, keep_fits = FALSE), error = function(e) { message("skip ", key, ": ", conditionMessage(e)); NULL })
  if (is.null(r)) next
  a <- r$assignments
  m <- matrix(code[a], nrow(a), ncol(a), dimnames = dimnames(a)); storage.mode(m) <- "integer"
  frac <- if (key %in% names(frac_of)) frac_of[[key]] else r$fraction
  cache$fractions <- rbind(cache$fractions, data.frame(file = key, fraction = frac, n_samples = nrow(a),
                                                       n_iter = ncol(a), n_drawn_mean = mean(rowSums(!is.na(a)))))
  cache$assignments[[key]] <- m
  cache$prob_enterotype[[key]] <- r$legacy_prob_enterotype
  cache$flags[[key]] <- r$flags
  cache$gap[[key]] <- r$gap
  cache$ties[[key]] <- r$legacy_ties
  if (is.null(cache$samples)) cache$samples <- rownames(a)
  cat(sprintf("%-40s Q=%.2f  n=%d  %.1f min\n", key, frac, nrow(a), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  rm(r, a, m); invisible(gc())
}
cache$fractions <- cache$fractions[order(cache$fractions$fraction), ]
saveRDS(cache, "out/cache/erpm_legacy.rds")
cat("\ncached", nrow(cache$fractions), "fractions to out/cache/erpm_legacy.rds (",
    round(file.size("out/cache/erpm_legacy.rds") / 1e6), "MB ) in",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
print(cache$fractions[, c("fraction", "n_samples", "n_drawn_mean")], row.names = FALSE)
