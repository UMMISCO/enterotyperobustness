#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 05  Annotation anomalies: Supplementary Tables S1 to S3 regenerated from the raw
#     files, compared cell by cell with the submitted tables, and the healthy-anchor
#     results re-aggregated with the flagged iterations excluded (Reviewer 1, major 1:
#     "state whether anomalous iterations were excluded or retained; if retained,
#     re-run with them excluded and show the trade-off survives").
#
# Inputs : mirror (Experiments/*_K=4_WARNING.csv, ConfidenceN=*/Conf*_K=4_WARNING.csv and
#          *_Samples50ProbOfBeingInTwoEnt.RData), out/cache/erpm_legacy.rds (script 00),
#          the manuscript's supplemental.tex for the comparison.
# Outputs: out/05_anomaly_tables/  tableS1.{csv,tex} tableS2.{csv,tex} tableS3.{csv,tex},
#          comparison_with_submitted.csv, excluded_iterations_healthy.csv, summary.txt
# ------------------------------------------------------------------------------
suppressMessages(pkgload::load_all(".", quiet = TRUE))
mirror <- Sys.getenv("ENTEROTYPES_MIRROR")
if (!nzchar(mirror) || !dir.exists(mirror)) stop("Set ENTEROTYPES_MIRROR to the directory holding the 2022 analyses (private, not distributed).")
art <- file.path(mirror, "analyses/2.Estimating_DMM/Article")
tex_file <- Sys.getenv("SUPPLEMENTAL_TEX")          # optional: the manuscript's supplemental.tex, for the comparison with the submitted tables
if (!nzchar(tex_file)) message("SUPPLEMENTAL_TEX not set: the comparison with the submitted tables is skipped.")
out <- "out/05_anomaly_tables"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
read_warning <- function(f) { w <- utils::read.table(f, header = TRUE, stringsAsFactors = FALSE); w[, paste0("V", 1:5)] != "0" }
count_flags <- function(w) c(BACT2 = sum(w[, 1]), RUM = sum(w[, 2]), n_iter = nrow(w))

# ---- 1. Table S1: ERP experiments ------------------------------------------------------
exp <- legacy_erp_experiments()
s1_label <- c("Healthy-Diabetes" = "Healthy--Diabetes", "Healthy-Healthy" = "Healthy--Healthy", "Healthy-Obese" = "Healthy--Obese",
              "Healthy-DiabetesT" = "Healthy--Diabetes T", "Healthy-DiabetesNT" = "Healthy--Diabetes NT",
              "Healthy-ObeseNT" = "Healthy--Obese NT", "Healthy-ObeseT" = "Healthy--Obese T",
              "DiabeticsT-Healthy" = "Diabetes T--Healthy", "DiabeticsT-Obese" = "Diabetes T--Obese",
              "DiabeticsT-DiabeticsT" = "Diabetes T--Diabetes T", "DiabeticsT-DiabeticsNT" = "Diabetes T--Diabetes NT",
              "Obese-Healthy" = "Obese--Healthy", "Obese-Obese" = "Obese--Obese", "Obese-DiabeticsT" = "Obese--Diabetes T",
              "Obese-DiabeticsNT" = "Obese--Diabetes NT", "Obese-Diabetics" = "Obese--Diabetes (not shown)")
s1 <- do.call(rbind, lapply(seq_len(nrow(exp)), function(i) {
  f <- file.path(art, "Experiments", exp$folder[i], paste0(exp$label[i], "_K=4_WARNING.csv"))
  f2 <- file.path(art, "Experiments/WarningsFolder_July2024", paste0(exp$label[i], "_K=4_WARNING.csv"))
  if (!file.exists(f)) return(NULL)
  cf <- count_flags(read_warning(f))
  data.frame(label = exp$label[i], row = s1_label[[exp$label[i]]], BACT2 = cf[["BACT2"]], RUM = cf[["RUM"]], n_iter = cf[["n_iter"]],
             same_as_july2024_copy = if (file.exists(f2)) unname(tools::md5sum(f) == tools::md5sum(f2)) else NA)
}))
s1 <- s1[match(names(s1_label), s1$label), ]; s1 <- s1[!is.na(s1$label), ]
write.csv(s1, file.path(out, "tableS1.csv"), row.names = FALSE)
writeLines(c("\\begin{tabular}{lcc}", "\\toprule", "Anchor--Complementary dataset & BACT2 & RUM \\\\", "\\midrule",
             sprintf("%-25s & %d & %d \\\\", s1$row[s1$label != "Obese-Diabetics"], s1$BACT2[s1$label != "Obese-Diabetics"], s1$RUM[s1$label != "Obese-Diabetics"]),
             "\\bottomrule", "\\end{tabular}"), file.path(out, "tableS1.tex"))

# ---- 2. Table S2: ERPM on the full cohort, from the cache -----------------------------------
cache <- readRDS("out/cache/erpm_legacy.rds")
s2 <- do.call(rbind, lapply(cache$fractions$file, function(k) {
  fl <- cache$flags[[k]]
  data.frame(fraction = cache$fractions$fraction[cache$fractions$file == k],
             BACT2 = sum(fl[["bact2_not_lowest_diversity"]]), RUM = sum(fl[["rum_not_highest_diversity"]]),
             W = if (is.null(cache$ties[[k]])) NA_integer_ else as.integer(cache$ties[[k]]),
             BACT2_rederived = sum(fl[["bact2_not_lowest_diversity_rederived"]]), RUM_rederived = sum(fl[["rum_not_highest_diversity_rederived"]]))
}))
s2 <- s2[order(s2$fraction), ]
write.csv(s2, file.path(out, "tableS2.csv"), row.names = FALSE)
four_col_table <- function(d, header, caption_note) {
  d <- d[order(d$fraction), ]; d$W <- ifelse(is.na(d$W), "NI", as.character(d$W))
  n <- nrow(d); per <- ceiling(n / 4)
  cols <- lapply(1:4, function(j) { idx <- ((j - 1) * per + 1):min(j * per, n); d[idx, ] })
  rows <- vapply(seq_len(per), function(r) paste(vapply(cols, function(cc) if (r <= nrow(cc)) sprintf("%.2f & %d & %d & %s", cc$fraction[r], cc$BACT2[r], cc$RUM[r], cc$W[r]) else " & & & ", character(1)), collapse = " & "), character(1))
  c("\\begin{tabular}{cccccccccccccccc}", "\\toprule", paste(rep(header, 4), collapse = " & "), "\\midrule",
    paste0(rows, " \\\\"), "\\bottomrule", "\\end{tabular}", paste0("% ", caption_note))
}
writeLines(four_col_table(s2, "\\% & BACT2 & RUM & W", "Q = 1.00 (no subsampling) has no run and is not listed"), file.path(out, "tableS2.tex"))

# ---- 3. Table S3: sub-cohorts ----------------------------------------------------------------
sub <- c("90" = "ConfidenceN=2000", "80" = "ConfidenceN=80P", "70" = "ConfidenceN=70P", "50" = "ConfidenceN=50P", "30" = "ConfidenceN=30P")
s3 <- do.call(rbind, lapply(names(sub), function(M) {
  d <- file.path(art, sub[[M]]); wf <- list.files(d, pattern = "_K=4_WARNING\\.csv$", full.names = TRUE)
  do.call(rbind, lapply(wf, function(f) {
    frac <- as.numeric(sub("^Conf_?", "", sub("_K=4_WARNING\\.csv$", "", basename(f))))
    tf <- file.path(d, sub("_K=4_WARNING\\.csv$", "_Samples50ProbOfBeingInTwoEnt.RData", basename(f)))
    W <- if (file.exists(tf)) { e <- new.env(); load(tf, envir = e); as.integer(e$W_50) } else NA_integer_
    cf <- count_flags(read_warning(f))
    data.frame(M = as.integer(M), fraction = frac, BACT2 = cf[["BACT2"]], RUM = cf[["RUM"]], W = W, n_iter = cf[["n_iter"]])
  }))
}))
s3 <- s3[order(s3$M, s3$fraction), ]
write.csv(s3, file.path(out, "tableS3.csv"), row.names = FALSE)
s3_shown <- s3[s3$fraction >= 0.5, ]
tex3 <- c("\\begin{tabular}{cccccccccccccccccccc}", "\\toprule",
          paste(sprintf("\\multicolumn{4}{c}{$M=%s$}", c(90, 80, 70, 50, 30)), collapse = " & "), "\\\\",
          paste(rep("\\% & BACT2 & RUM & W", 5), collapse = " & "), "\\\\", "\\midrule")
fr <- sort(unique(s3_shown$fraction))
for (q in fr) {
  cells <- vapply(c(90, 80, 70, 50, 30), function(M) { r <- s3_shown[s3_shown$M == M & s3_shown$fraction == q, ]
    if (nrow(r)) sprintf("%.2f & %d & %d & %s", q, r$BACT2, r$RUM, ifelse(is.na(r$W), "NI", r$W)) else sprintf("%.2f & -- & -- & --", q) }, character(1))
  tex3 <- c(tex3, paste0(paste(cells, collapse = " & "), " \\\\"))
}
writeLines(c(tex3, "\\bottomrule", "\\end{tabular}"), file.path(out, "tableS3.tex"))

# ---- 4. Comparison with the submitted tables (parsed from supplemental.tex) ------------------
cmp <- NULL
if (nzchar(tex_file) && file.exists(tex_file)) {
  tx <- readLines(tex_file)
  cells <- function(line) trimws(strsplit(sub("\\\\\\\\.*$", "", line), "&")[[1]])
  # S1
  i1 <- grep("^Healthy--Diabetes +&", tx)[1]; rows1 <- tx[i1:(i1 + 14)]
  sub1 <- do.call(rbind, lapply(rows1, function(l) { c <- cells(l); data.frame(row = c[1], BACT2_sub = as.integer(c[2]), RUM_sub = as.integer(c[3])) }))
  m1 <- merge(s1[, c("row", "BACT2", "RUM")], sub1, by = "row")
  cmp <- rbind(cmp, data.frame(table = "S1", key = m1$row, BACT2_new = m1$BACT2, BACT2_sub = m1$BACT2_sub, RUM_new = m1$RUM, RUM_sub = m1$RUM_sub, W_new = NA, W_sub = NA))
  # S2: rows of 4 groups x (frac, BACT2, RUM, W)
  i2 <- grep("^0\\.05 & 13 & 142", tx)[1]; rows2 <- tx[i2:(i2 + 13)]
  sub2 <- do.call(rbind, lapply(rows2, function(l) { c <- cells(l); do.call(rbind, lapply(0:3, function(g) data.frame(fraction = as.numeric(c[4 * g + 1]), BACT2_sub = as.integer(c[4 * g + 2]), RUM_sub = as.integer(c[4 * g + 3]), W_sub = c[4 * g + 4]))) }))
  m2 <- merge(s2[, c("fraction", "BACT2", "RUM", "W")], sub2, by = "fraction")
  cmp <- rbind(cmp, data.frame(table = "S2", key = sprintf("Q=%.2f", m2$fraction), BACT2_new = m2$BACT2, BACT2_sub = m2$BACT2_sub, RUM_new = m2$RUM, RUM_sub = m2$RUM_sub, W_new = m2$W, W_sub = suppressWarnings(as.integer(m2$W_sub))))
  # S3: rows of 5 groups
  i3 <- grep("^0\\.50 & 0 & 2 & 0 & 0\\.50", tx)[1]; rows3 <- tx[i3:(i3 + 10)]
  Ms <- c(90, 80, 70, 50, 30)
  sub3 <- do.call(rbind, lapply(rows3, function(l) { c <- cells(l); do.call(rbind, lapply(0:4, function(g) data.frame(M = Ms[g + 1], fraction = as.numeric(c[4 * g + 1]), BACT2_sub = as.integer(c[4 * g + 2]), RUM_sub = as.integer(c[4 * g + 3]), W_sub = c[4 * g + 4]))) }))
  m3 <- merge(s3[, c("M", "fraction", "BACT2", "RUM", "W")], sub3, by = c("M", "fraction"))
  cmp <- rbind(cmp, data.frame(table = "S3", key = sprintf("M=%d Q=%.2f", m3$M, m3$fraction), BACT2_new = m3$BACT2, BACT2_sub = m3$BACT2_sub, RUM_new = m3$RUM, RUM_sub = m3$RUM_sub, W_new = m3$W, W_sub = suppressWarnings(as.integer(m3$W_sub))))
  cmp$BACT2_ok <- cmp$BACT2_new == cmp$BACT2_sub; cmp$RUM_ok <- cmp$RUM_new == cmp$RUM_sub
  cmp$W_ok <- ifelse(is.na(cmp$W_sub) & is.na(cmp$W_new), TRUE, cmp$W_new == cmp$W_sub)
  write.csv(cmp, file.path(out, "comparison_with_submitted.csv"), row.names = FALSE)
}

# ---- 5. Healthy-anchor experiments with flagged iterations excluded (A1d) -----------------
exp_dir <- file.path(art, "Experiments")
healthy <- c("Healthy-Healthy", "Healthy-Obese", "Healthy-Diabetes", "Healthy-DiabetesT", "Healthy-DiabetesNT", "Healthy-ObeseNT", "Healthy-ObeseT")
excl <- list(); prev_tab <- list()
for (lab in healthy) {
  r <- read_legacy_erp(exp_dir, lab, keep_fits = FALSE)
  anchors <- r$samples$anchors
  flagged <- r$flags$bact2_not_lowest_diversity | r$flags$rum_not_highest_diversity
  a_all <- r$assignments[anchors, ]; a_ex <- r$assignments[anchors, !flagged, drop = FALSE]
  c_all <- robustness(assignment_frequency(a_all)); c_ex <- robustness(assignment_frequency(a_ex))
  prev_tab[[lab]] <- data.frame(experiment = lab, n_flagged = sum(flagged),
    do.call(cbind, lapply(enterotype_names(), function(e) stats::setNames(data.frame(100 * mean(c_all$enterotype == e), 100 * mean(c_ex$enterotype == e)), paste0(e, c("_all", "_excl"))))),
    consensus_changed = sum(c_all$enterotype != c_ex$enterotype),
    max_abs_dfreq = max(abs(assignment_frequency(a_all) - assignment_frequency(a_ex))))
  excl[[lab]] <- c_all$enterotype; names(excl[[lab]]) <- anchors
}
prev_tab <- do.call(rbind, prev_tab); rownames(prev_tab) <- NULL
write.csv(prev_tab, file.path(out, "excluded_iterations_healthy.csv"), row.names = FALSE)

# ---- 6. Summary --------------------------------------------------------------------------------
sink(file.path(out, "summary.txt"))
cat("05  Annotation anomaly tables, run", format(Sys.time()), "\n\n")
cat("Table S1 (ERP), BACT2 / RUM flags over 500 iterations:\n"); print(s1[, c("row", "BACT2", "RUM", "n_iter", "same_as_july2024_copy")], row.names = FALSE)
cat(sprintf("\nTable S2 (ERPM full cohort): %d fractions; stored vs re-derived flags identical in %d / %d cells\n",
            nrow(s2), sum(s2$BACT2 == s2$BACT2_rederived & s2$RUM == s2$RUM_rederived), nrow(s2)))
cat(sprintf("Table S3 (sub-cohorts): %d rows over M = %s\n", nrow(s3), paste(unique(s3$M), collapse = ", ")))
if (!is.null(cmp)) {
  cat("\nComparison with the submitted tables (cells equal / total):\n")
  print(aggregate(cbind(BACT2_ok, RUM_ok, W_ok) ~ table, cmp, function(x) sprintf("%d/%d", sum(x), length(x))), row.names = FALSE)
  bad <- cmp[!(cmp$BACT2_ok & cmp$RUM_ok & cmp$W_ok), ]
  if (nrow(bad)) { cat("\nDiscrepant cells:\n"); print(bad[, c("table", "key", "BACT2_new", "BACT2_sub", "RUM_new", "RUM_sub", "W_new", "W_sub")], row.names = FALSE) } else cat("No discrepant cell.\n")
}
cat("\nHealthy-anchor experiments, all iterations vs flagged iterations excluded (prevalence %, 142 anchors):\n")
print(prev_tab, row.names = FALSE, digits = 3)
cat("\nSession:\n"); print(sessionInfo())
sink()
cat(readLines(file.path(out, "summary.txt"), n = 45), sep = "\n")
