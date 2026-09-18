# Mixed-effects analyses of the ERP designs: the same anchor individuals are
# observed in every condition, so the individual enters as a random intercept.
# lme4 / lmerTest / emmeans / mclogit / nnet are Suggests.

#' Linear mixed model of entropy across ERP conditions
#'
#' Fits `value ~ condition * label + (1 | sample)` (or `value ~ condition +
#' (1 | sample)` without `label`) with \pkg{lme4}, so that the paper's framing
#' "entropy among Bact2 individuals" can keep the enterotype assigned **in each
#' condition** as a covariate while the random intercept accounts for the
#' repeated individual. Omnibus tests are type III F-tests with Satterthwaite
#' degrees of freedom from \pkg{lmerTest} (whose type III construction does
#' not depend on the contrast coding); pairwise contrasts between conditions
#' within each label come from \pkg{emmeans}, adjusted with `p_adjust` within
#' label. The same model on rank-transformed values is fitted as a robustness
#' check, since entropy is bounded and right-skewed.
#'
#' @inheritParams paired_entropy_test
#' @param label Optional column: the enterotype assigned to the sample in that
#'   condition (varies within sample). `NULL` for a condition-only model.
#' @param rank_check Also fit the model on `rank(value)` and report its omnibus table.
#' @return A list with `model`, `formula`, `omnibus` (term, F, NumDF, DenDF, p),
#'   `pairwise` (label, cond1, cond2, estimate, se, df, t, p, p_adj), `variance`
#'   (random-intercept and residual variances, ICC), `converged`, `messages`,
#'   `n_samples`, and `omnibus_rank` when `rank_check` is `TRUE`.
#' @export
mixed_entropy_model <- function(data, sample = "sample", condition = "condition",
                                value = "entropy", label = NULL, p_adjust = "BH",
                                rank_check = TRUE) {
  for (pkg in c("lme4", "lmerTest", "emmeans")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Package ", pkg, " is needed.")
  }
  d <- data.frame(.y = data[[value]], .cond = factor(data[[condition]]),
                  .id = factor(data[[sample]]), stringsAsFactors = FALSE)
  if (!is.null(label)) d$.lab <- factor(data[[label]])
  f <- if (is.null(label)) .y ~ .cond + (1 | .id) else .y ~ .cond * .lab + (1 | .id)
  msgs <- character(0)
  m <- withCallingHandlers(
    lmerTest::lmer(f, data = d, REML = TRUE),
    warning = function(w) { msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning") })
  conv <- m@optinfo$conv$lme4$messages
  converged <- is.null(conv) || length(conv) == 0
  omni <- as.data.frame(stats::anova(m, type = 3, ddf = "Satterthwaite"))
  omni <- data.frame(term = gsub("\\.cond", "condition", gsub("\\.lab", "label", rownames(omni))),
                     F = omni[["F value"]], NumDF = omni$NumDF, DenDF = omni$DenDF,
                     p = omni[["Pr(>F)"]], row.names = NULL)
  emm <- if (is.null(label)) emmeans::emmeans(m, ~ .cond) else emmeans::emmeans(m, ~ .cond | .lab)
  ct <- as.data.frame(emmeans::contrast(emm, method = "pairwise", adjust = "none"))
  split_pair <- function(x) {
    x <- gsub("[()]", "", as.character(x))
    parts <- strsplit(x, " - ", fixed = TRUE)
    list(cond1 = vapply(parts, `[`, character(1), 1), cond2 = vapply(parts, `[`, character(1), 2))
  }
  sp <- split_pair(ct$contrast)
  pw <- data.frame(stratum = if (is.null(label)) "all" else as.character(ct$.lab),
                   cond1 = sp$cond1, cond2 = sp$cond2,
                   estimate = ct$estimate, se = ct$SE, df = ct$df, t = ct$t.ratio, p = ct$p.value,
                   stringsAsFactors = FALSE)
  pw$p_adj <- stats::ave(pw$p, pw$stratum, FUN = function(v) stats::p.adjust(v, method = p_adjust))
  vc <- as.data.frame(lme4::VarCorr(m))
  v_id <- vc$vcov[vc$grp == ".id"]; v_res <- vc$vcov[vc$grp == "Residual"]
  out <- list(model = m,
              formula = gsub("\\.y", value, gsub("\\.cond", condition, gsub("\\.lab", if (is.null(label)) "" else label,
                             gsub("\\.id", sample, deparse(f))))),
              omnibus = omni, pairwise = pw,
              variance = data.frame(random_intercept = v_id, residual = v_res, icc = v_id / (v_id + v_res)),
              converged = converged, messages = c(msgs, conv),
              n_samples = nlevels(d$.id), p_adjust = p_adjust)
  if (rank_check) {
    d$.y <- rank(d$.y)
    mr <- suppressWarnings(lmerTest::lmer(f, data = d, REML = TRUE))
    or <- as.data.frame(stats::anova(mr, type = 3, ddf = "Satterthwaite"))
    out$omnibus_rank <- data.frame(term = gsub("\\.cond", "condition", gsub("\\.lab", "label", rownames(or))),
                                   F = or[["F value"]], NumDF = or$NumDF, DenDF = or$DenDF, p = or[["Pr(>F)"]],
                                   row.names = NULL)
  }
  class(out) <- "mixed_entropy_model"
  out
}

#' @export
print.mixed_entropy_model <- function(x, ...) {
  cat("Linear mixed model:", x$formula, "\n")
  cat(sprintf("  %d individuals; random-intercept variance %.4f, residual %.4f, ICC %.2f; converged: %s\n",
              x$n_samples, x$variance$random_intercept, x$variance$residual, x$variance$icc, x$converged))
  print(x$omnibus, row.names = FALSE, digits = 4)
  cat(sprintf("  %d pairwise contrasts, %d with adjusted p < 0.05 (%s)\n",
              nrow(x$pairwise), sum(x$pairwise$p_adj < 0.05, na.rm = TRUE), x$p_adjust))
  invisible(x)
}

#' Multinomial mixed logit of enterotype labels across ERP conditions
#'
#' `mclogit::mblogit(label ~ condition, random = ~ 1 | sample)`: a baseline-
#' category logit with a random intercept per individual, estimated by PQL
#' (or MQL). Replaces the Chi-square on the condition x label table, which
#' treated the repeated anchors as independent. The effect of `condition` is
#' tested with a joint Wald test on all its coefficients (the deviance-based
#' `anova` is not valid for a PQL fit), and per enterotype logit with BH.
#'
#' @param labels Character matrix, samples x conditions.
#' @param reference_condition Baseline condition (the anchors clustered with
#'   their own group in the paper). Default: the first column.
#' @param reference_enterotype Baseline category of the logit.
#' @param method `"PQL"` (default) or `"MQL"`.
#' @param p_adjust Adjustment across the per-enterotype tests.
#' @return A list with `converged`, `messages`, `wald` (statistic, df, p),
#'   `per_enterotype` (enterotype, statistic, df, p, p_adj), `coefficients`
#'   (enterotype, condition, estimate, se, z, p; log-odds vs the reference
#'   condition), `model`. When the fit fails, `converged` is `FALSE`, the tables
#'   are `NULL` and `messages` says why.
#' @export
multinomial_mixed_test <- function(labels, reference_condition = colnames(labels)[1],
                                   reference_enterotype = "Bact1", method = c("PQL", "MQL"),
                                   p_adjust = "BH") {
  method <- match.arg(method)
  if (!requireNamespace("mclogit", quietly = TRUE)) stop("Package mclogit is needed.")
  d <- long_labels(labels, reference_condition, reference_enterotype)
  msgs <- character(0)
  fit <- tryCatch(withCallingHandlers(
    mclogit::mblogit(label ~ condition, random = ~ 1 | sample, data = d, method = method,
                     control = mclogit::mmclogit.control(trace = FALSE)),
    warning = function(w) { msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(e) { msgs <<- c(msgs, conditionMessage(e)); NULL })
  fail <- function(why) list(converged = FALSE, messages = unique(c(msgs, why)), wald = NULL,
                             per_enterotype = NULL, coefficients = NULL, model = fit, method = method)
  if (is.null(fit)) return(fail("mblogit failed"))
  b <- stats::coef(fit); V <- tryCatch(stats::vcov(fit), error = function(e) NULL)
  if (is.null(V) || any(!is.finite(V)) || any(!is.finite(b))) return(fail("non-finite coefficients or covariance"))
  if (isFALSE(fit$converged)) msgs <- c(msgs, "mblogit reports non-convergence")
  # PQL/MQL fits whose inner (random-effects) iterations did not converge are not
  # trustworthy: with anchors whose label never changes the random intercept is
  # weakly identified. Treat any convergence warning as failure so callers fall back.
  inner_fail <- any(grepl("conve?rge", msgs, ignore.case = TRUE))
  idx <- which(grepl("~condition", names(b), fixed = TRUE))   # coefficient names are "<enterotype>~condition<level>"
  wald <- function(i) {
    bb <- b[i]; VV <- V[i, i, drop = FALSE]
    W <- tryCatch(as.numeric(t(bb) %*% solve(VV) %*% bb), error = function(e) NA_real_)
    data.frame(statistic = W, df = length(bb), p = stats::pchisq(W, length(bb), lower.tail = FALSE))
  }
  joint <- wald(idx)
  ent <- sub("~.*", "", names(b)[idx])
  per <- do.call(rbind, lapply(unique(ent), function(e) cbind(enterotype = e, wald(idx[ent == e]))))
  per$p_adj <- stats::p.adjust(per$p, method = p_adjust)
  se <- sqrt(diag(V))[idx]
  coefs <- data.frame(enterotype = ent, condition = sub(".*~condition", "", names(b)[idx]),
                      estimate = unname(b[idx]), se = unname(se), z = unname(b[idx] / se),
                      p = 2 * stats::pnorm(-abs(unname(b[idx] / se))), stringsAsFactors = FALSE)
  list(converged = !isFALSE(fit$converged) && !inner_fail && is.finite(joint$statistic), messages = unique(msgs),
       wald = joint, per_enterotype = per, coefficients = coefs, model = fit, method = method,
       reference = c(condition = reference_condition, enterotype = reference_enterotype))
}

#' Permutation test of a multinomial Wald statistic with individuals as blocks
#'
#' Fallback for [multinomial_mixed_test()] when the mixed logit does not
#' converge. Fits `nnet::multinom(label ~ condition)` on the long data, takes
#' the joint Wald statistic of the `condition` coefficients, and compares it
#' with its distribution under within-individual permutation of the conditions
#' (each individual's labels shuffled across its own conditions). No random
#' effect is needed because the permutation scheme itself respects the blocks.
#'
#' @inheritParams multinomial_mixed_test
#' @param n_perm Number of permutations.
#' @param seed Seed.
#' @return A list with `statistic`, `df`, `p`, `n_perm`, `coefficients`.
#' @export
multinomial_permutation_test <- function(labels, reference_condition = colnames(labels)[1],
                                         reference_enterotype = "Bact1", n_perm = 999, seed = 314) {
  if (!requireNamespace("nnet", quietly = TRUE)) stop("Package nnet is needed.")
  d <- long_labels(labels, reference_condition, reference_enterotype)
  stat <- function(dd) {
    fit <- nnet::multinom(label ~ condition, data = dd, trace = FALSE, maxit = 500)
    b <- as.vector(t(stats::coef(fit))); nm <- as.vector(t(outer(rownames(stats::coef(fit)), colnames(stats::coef(fit)), paste, sep = ":")))
    V <- stats::vcov(fit)
    sel <- grepl("condition", nm)
    W <- tryCatch(as.numeric(t(b[sel]) %*% solve(V[sel, sel]) %*% b[sel]), error = function(e) NA_real_)
    list(W = W, df = sum(sel), coef = stats::setNames(b, nm))
  }
  obs <- stat(d)
  set.seed(seed)
  k <- ncol(labels); n <- nrow(labels)
  perm <- vapply(seq_len(n_perm), function(b) {
    m <- labels
    for (i in seq_len(n)) m[i, ] <- labels[i, sample.int(k)]
    stat(long_labels(m, reference_condition, reference_enterotype))$W
  }, numeric(1))
  list(statistic = obs$W, df = obs$df, p = (1 + sum(perm >= obs$W, na.rm = TRUE)) / (1 + n_perm),
       n_perm = n_perm, coefficients = obs$coef)
}

#' Long data frame from a samples x conditions label matrix (internal)
#' @keywords internal
long_labels <- function(labels, reference_condition, reference_enterotype) {
  if (!reference_condition %in% colnames(labels)) stop("`reference_condition` is not a column of `labels`.")
  d <- data.frame(sample = rep(rownames(labels), times = ncol(labels)),
                  condition = rep(colnames(labels), each = nrow(labels)),
                  label = as.vector(labels), stringsAsFactors = FALSE)
  d$condition <- stats::relevel(factor(d$condition), ref = reference_condition)
  lev <- unique(c(reference_enterotype, sort(unique(d$label))))
  d$label <- factor(d$label, levels = lev[lev %in% d$label])
  d$sample <- factor(d$sample)
  d
}
