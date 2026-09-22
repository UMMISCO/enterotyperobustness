#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 00b Read the ERPM sub-cohort environments once (Supplementary Figure S3 and
#     Table S3: ERPM repeated on sub-cohorts of the MetaCardis dataset, one
#     ConfidenceN=* folder per sub-cohort) and cache the integer-coded
#     assignments, as script 00 does for the full cohort.
#
# Cache: out/cache/erpm_subcohorts.rds, a list with
#   runs        : data.frame(subcohort, file, fraction, n_samples, n_iter, n_drawn_mean)
#   assignments : list per run (key = "<subcohort>/<file>") of integer matrix
#                 samples x iterations (1 = Bact1, 2 = Bact2, 3 = Prevotella, 4 = Ruminococcus)
#   samples     : list per subcohort of the sample IDs
# ------------------------------------------------------------------------------
suppressMessages(pkgload::load_all(".", quiet = TRUE)); source("inst/scripts/cache_utils.R")
mirror <- Sys.getenv("ENTEROTYPES_MIRROR")
if (!nzchar(mirror) || !dir.exists(mirror)) stop("Set ENTEROTYPES_MIRROR to the directory holding the 2022 analyses (private, not distributed).")
art <- file.path(mirror, "analyses/2.Estimating_DMM/Article")
dir.create("out/cache", showWarnings = FALSE, recursive = TRUE)
files <- subcohort_files(art)
cat(length(files), "sub-cohort environments in", length(subcohort_dirs(art)), "folders\n")
manifest <- subcohort_cache_manifest(art)
code <- setNames(1:4, enterotype_names())
cache <- list(runs = NULL, assignments = list(), samples = list())
t0 <- Sys.time()
for (f in files) {
  sub <- basename(dirname(f)); key <- file.path(sub, basename(f))
  r <- tryCatch(read_legacy_erpm(f, keep_fits = FALSE), error = function(e) { message("skip ", key, ": ", conditionMessage(e)); NULL })
  if (is.null(r)) next
  a <- r$assignments
  m <- matrix(code[a], nrow(a), ncol(a), dimnames = dimnames(a)); storage.mode(m) <- "integer"
  frac <- if (!is.na(r$fraction)) r$fraction else as.numeric(sub("^Conf_?([0-9.]+)_.*$", "\\1", basename(f)))
  cache$runs <- rbind(cache$runs, data.frame(subcohort = sub, file = basename(f), fraction = frac, n_samples = nrow(a),
                                             n_iter = ncol(a), n_drawn_mean = mean(rowSums(!is.na(a)))))
  cache$assignments[[key]] <- m
  if (is.null(cache$samples[[sub]])) cache$samples[[sub]] <- rownames(a)
  cat(sprintf("%-50s Q=%.2f  n=%d  %.1f min\n", key, frac, nrow(a), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  rm(r, a, m); invisible(gc())
}
cache$runs <- cache$runs[order(cache$runs$subcohort, cache$runs$fraction), ]
cache_write(cache, "out/cache/erpm_subcohorts.rds", manifest)
cat("\ncached", nrow(cache$runs), "runs to out/cache/erpm_subcohorts.rds in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
print(aggregate(fraction ~ subcohort + n_samples, cache$runs, function(v) sprintf("%d fractions, %.2f to %.2f", length(v), min(v), max(v))), row.names = FALSE)
