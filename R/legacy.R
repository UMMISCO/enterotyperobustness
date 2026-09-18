#' The ERP experiments of the paper, as stored by the 2022-2023 scripts
#'
#' Labels and folders under `.../2.Estimating_DMM/Article/Experiments/`, and
#' which paper panel each feeds.
#'
#' @return A data frame with `label`, `folder`, `anchor`, `pool`, `panel`.
#' @export
legacy_erp_experiments <- function() {
  data.frame(
    label = c("Healthy-Healthy", "Healthy-Obese", "Healthy-Diabetes",
              "Healthy-DiabetesT", "Healthy-DiabetesNT", "Healthy-ObeseT", "Healthy-ObeseNT",
              "DiabeticsT-Healthy", "DiabeticsT-Obese", "DiabeticsT-DiabeticsT", "DiabeticsT-DiabeticsNT",
              "Obese-Healthy", "Obese-Obese", "Obese-Diabetics", "Obese-DiabeticsT", "Obese-DiabeticsNT"),
    folder = c(".", "Obese_T_All2a2b", ".", "Diabetes_T_NT", "Diabetes_T_NT", "Obese_T_T", "Obese_T_NT",
               rep("Experiment_CoreFixedDiabetics", 4), rep("Experiment_CoreFixedObese", 5)),
    anchor = c(rep("healthy (142)", 7), rep("T2D metformin+ (226)", 4), rep("obese (209)", 5)),
    pool = c("other healthy", "obese 2a+2b", "T2D", "T2D metformin+", "T2D metformin-",
             "obese statin+ (33, with replacement)", "obese statin-",
             "healthy", "obese", "T2D metformin+", "T2D metformin-",
             "healthy", "other obese", "T2D", "T2D metformin+", "T2D metformin-"),
    panel = c("Fig1 A-C", "Fig1 A-C", "Fig1 A-C", "Fig1 D-F", "Fig1 D-F", "Fig1 D-F", "Fig1 D-F",
              rep("Fig2 A-C", 4), "Fig2 D-F", "Fig2 D-F", "not shown", "Fig2 D-F", "Fig2 D-F"),
    stringsAsFactors = FALSE)
}

#' Read one ERP experiment stored by the 2022-2023 scripts
#'
#' Rebuilds an `erp_result` from the files the original pipeline wrote:
#' the 500 `DMN` fits (`<label>_EachIterationDMMestimation.RData`, or, when that
#' file is absent as for `Healthy-Diabetes`, the `result` object inside
#' `<label>_TotalEnvironment.RData`), the per-iteration name-to-component maps,
#' the draws, the anchor count matrix and the `WARNING` flags. Labels are
#' rebuilt from each fit's MAP component and the stored map, and checked
#' against the stored `EnterotypeAssigment` vectors.
#'
#' @param root Path to the `Experiments/` directory of the mirror.
#' @param label Experiment label, see [legacy_erp_experiments()].
#' @param folder Sub-folder; looked up in [legacy_erp_experiments()] if `NULL`.
#' @param abundance Optional full abundance matrix (samples x features). If
#'   given, every iteration is re-annotated with [annotate_clusters()] so that
#'   `gap` is available, and the re-derived names are compared with the stored
#'   ones.
#' @param marker_features,diversity,extra_checks Passed to [annotate_clusters()]
#'   when `abundance` is given.
#' @param keep_fits Keep the `DMN` objects in the result.
#' @return An `erp_result` with two extra elements: `legacy_agreement`, the
#'   fraction of stored anchor labels reproduced from the fits, and
#'   `legacy_files`, the files read.
#' @export
read_legacy_erp <- function(root, label, folder = NULL, abundance = NULL,
                            marker_features = c(bacteroides = "Bacteroides", prevotella = "Prevotella"),
                            diversity = NULL, extra_checks = NULL, keep_fits = TRUE) {
  if (is.null(folder)) {
    tab <- legacy_erp_experiments()
    if (!label %in% tab$label) stop("Unknown label; see legacy_erp_experiments().")
    folder <- tab$folder[tab$label == label]
  }
  dir <- file.path(root, folder)
  f <- function(suffix) file.path(dir, paste0(label, suffix))
  files <- c(fits = f("_EachIterationDMMestimation.RData"),
             maps = f("_EachIterationFinalEnterotypes.RData"),
             index = f("index.RData"),
             fixed = f("_K=4_Fixed_Dataset.RData"),
             warning = f("_K=4_WARNING.csv"),
             assignment = f("_EnterotypeAssigment.RData"),
             environment = f("_TotalEnvironment.RData"))
  env <- new.env()
  if (file.exists(files[["fits"]])) {
    load(files[["fits"]], envir = env)
  } else if (file.exists(files[["environment"]])) {
    message("No fits file for ", label, "; reading `result` from the TotalEnvironment (slow).")
    load(files[["environment"]], envir = env)
  } else stop("No DMM fits found for ", label, " in ", dir)
  for (k in c("maps", "index", "fixed", "assignment")) {
    if (!file.exists(files[[k]])) stop("Missing legacy file: ", files[[k]])
    load(files[[k]], envir = env)
  }
  fits <- env$result
  maps <- env$EachIterationFinalEnterotypes
  index <- env$index
  anchors <- rownames(env$Fixed_Dataset)
  n_iter <- length(fits)
  if (length(maps) != n_iter || length(index) != n_iter) stop("Inconsistent iteration counts in ", label)

  runs <- lapply(seq_len(n_iter), function(i) {
    fit <- fits[[i]]
    map <- unlist(maps[[i]])[enterotype_names()]
    component_names <- character(4); component_names[map] <- names(map)
    post <- dmm_posterior(fit)
    ids <- rownames(post)
    colnames(post) <- component_names
    post <- post[, enterotype_names(), drop = FALSE]
    labels <- component_names[dmm_assign(fit)]
    keep <- !duplicated(ids)
    ann <- NULL
    if (!is.null(abundance)) {
      ann <- annotate_clusters(fit, abundance[ids, , drop = FALSE], marker_features, diversity,
                               method = "rank", extra_checks = extra_checks)
    }
    list(ids = ids[keep], labels = labels[keep], posterior = post[keep, , drop = FALSE],
         map = map,
         flags = if (is.null(ann)) data.frame(row.names = 1) else ann$flags,
         gap = if (is.null(ann)) data.frame(row.names = 1) else ann$gap,
         reannotated_map = if (is.null(ann)) NULL else ann$map,
         fit = if (keep_fits) fit else NULL)
  })
  all_ids <- unique(c(anchors, unlist(lapply(runs, `[[`, "ids"))))
  out <- assemble_result(runs, all_samples = all_ids, keep_fits = keep_fits)

  # Stored WARNING flags (500 x 5, "0" = no anomaly)
  if (file.exists(files[["warning"]])) {
    w <- utils::read.table(files[["warning"]], header = TRUE, stringsAsFactors = FALSE)
    w <- w[, paste0("V", 1:5)] != "0"
    colnames(w) <- c("bact2_not_lowest_diversity", "rum_not_highest_diversity",
                     "rum_not_highest_FS", "rum_not_highest_FSR", "rum_not_highest_FSRE")
    stored <- cbind(iteration = seq_len(n_iter), as.data.frame(w))
    extra <- out$flags[, setdiff(names(out$flags), "iteration"), drop = FALSE]
    out$flags <- if (ncol(extra)) cbind(stored, extra) else stored   # cbind(NULL) would fail
  }
  # Agreement between rebuilt labels and the stored EnterotypeAssigment vectors (anchors only)
  stored <- env$EnterotypeAssigment
  agreement <- vapply(anchors, function(s) {
    v <- stored[[s]]; v <- v[v != "NotSelected"]
    if (length(v) != n_iter) return(NA_real_)
    mean(v == out$assignments[s, ])
  }, numeric(1))
  if (!is.null(abundance)) {
    same_map <- vapply(runs, function(r) identical(unname(r$map), unname(r$reannotated_map)), logical(1))
    out$reannotation_agreement <- mean(same_map)
  }
  out$index <- index
  out$samples <- list(anchors = anchors, complementary = setdiff(all_ids, anchors))
  out$legacy_agreement <- mean(agreement, na.rm = TRUE)
  out$legacy_files <- files[file.exists(files)]
  out$params <- list(procedure = "ERP (legacy)", label = label, n_iter = n_iter, K = 4,
                     n_complementary = length(index[[1]]), annotate = "rank (stored)")
  class(out) <- c("erp_result", "dmm_iterations")
  out
}

#' Read one ERPM fraction stored by the 2022 scripts
#'
#' Loads a `Conf<Q>_TotalEnvironment.RData` file into a private environment
#' and rebuilds an `erpm_fraction` object from `result`, the name maps, the
#' draws, `WARNING` and `W_50`. The environment also contains the abundance
#' table used (`Data_Used`), so the annotation is re-derived and `gap` is
#' filled without any extra input.
#'
#' @param file Path to a `*_TotalEnvironment.RData` file.
#' @param marker_features,diversity,extra_checks Passed to [annotate_clusters()].
#' @param keep_fits Keep the `DMN` objects.
#' @return An `erpm_fraction` (see [erpm()]) with extra elements
#'   `legacy_prob_enterotype` (the stored per-enterotype robustness lists, i.e.
#'   what Figure 3 plotted), `legacy_ties` (`W_50`) and `abundance` (the table).
#' @export
read_legacy_erpm <- function(file, marker_features = c(bacteroides = "Bacteroides", prevotella = "Prevotella"),
                             diversity = NULL, extra_checks = NULL, keep_fits = FALSE) {
  env <- new.env()
  load(file, envir = env)
  need <- c("result", "EachIterationFinalEnterotypes", "index", "Data_Used")
  miss <- setdiff(need, ls(env))
  if (length(miss)) stop("Not an ERPM environment (missing ", paste(miss, collapse = ", "), "): ", file)
  fits <- env$result; maps <- env$EachIterationFinalEnterotypes; index <- env$index
  abundance <- env$Data_Used
  n_iter <- length(fits)
  runs <- lapply(seq_len(n_iter), function(i) {
    fit <- fits[[i]]
    map <- unlist(maps[[i]])[enterotype_names()]
    component_names <- character(4); component_names[map] <- names(map)
    post <- dmm_posterior(fit)
    ids <- rownames(post)
    colnames(post) <- component_names
    post <- post[, enterotype_names(), drop = FALSE]
    ann <- annotate_clusters(fit, abundance[ids, , drop = FALSE], marker_features, diversity,
                             method = "rank", extra_checks = extra_checks)
    list(ids = ids, labels = component_names[dmm_assign(fit)], posterior = post,
         map = map, flags = ann$flags, gap = ann$gap, fit = if (keep_fits) fit else NULL)
  })
  out <- assemble_result(runs, all_samples = rownames(abundance), keep_fits = keep_fits)
  if (!is.null(env$WARNING)) {
    w <- env$WARNING[, 1:2] != "0"
    out$flags$stored_bact2_not_lowest_diversity <- w[, 1]
    out$flags$stored_rum_not_highest_diversity <- w[, 2]
  }
  out$index <- index
  out$fraction <- if (!is.null(env$Percen)) env$Percen else NA_real_
  out$samples <- list(cohort = rownames(abundance))
  out$legacy_prob_enterotype <- env$Prob_Enterotype
  out$legacy_ties <- env$W_50
  out$abundance <- abundance
  out$params <- list(procedure = "ERPM (legacy)", file = basename(file), n_iter = n_iter, K = 4,
                     fraction = out$fraction)
  class(out) <- c("erpm_fraction", "dmm_iterations")
  out
}
