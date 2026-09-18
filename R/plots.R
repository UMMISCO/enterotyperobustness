# Figure panels. ggplot2 and ggalluvial are Suggests: the package works without
# them, these functions need them.

#' Default enterotype colours (the palette of the submitted figures)
#' @export
enterotype_colours <- function() {
  c(Bact2 = "#F2C57C", Bact1 = "#E08E45", Prevotella = "#FF6B6B", Ruminococcus = "#6E2C00")
}

#' Bar plot of mean entropy per condition with significance brackets
#'
#' Reproduces panels C and F of Figures 1 and 2: mean +/- standard error of the
#' entropy per condition, one facet per stratum, with brackets drawn from a
#' `pairwise` table (e.g. `paired_entropy_test()$pairwise`).
#'
#' @inheritParams paired_entropy_test
#' @param pairwise Data frame with `stratum`, `cond1`, `cond2`, `p_adj` (and
#'   optionally `p`). Brackets are drawn for pairs with `p_adj < alpha`,
#'   labelled `"**"`; if `also_raw` is `TRUE`, pairs with only `p < alpha` get `"*"`.
#' @param alpha Significance level.
#' @param also_raw Draw single-star brackets for raw p < alpha.
#' @param colours Named vector of fill colours per stratum.
#' @param condition_levels Order of the conditions on the x axis.
#' @return A ggplot object.
#' @export
plot_entropy <- function(data, sample = "sample", condition = "condition", value = "entropy",
                         stratum = NULL, pairwise = NULL, alpha = 0.05, also_raw = FALSE,
                         colours = enterotype_colours(), condition_levels = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package ggplot2 is needed.")
  d <- data
  d$.stratum <- if (is.null(stratum)) "all" else d[[stratum]]
  d$.cond <- if (is.null(condition_levels)) factor(d[[condition]]) else factor(d[[condition]], levels = condition_levels)
  agg <- stats::aggregate(d[[value]], by = list(stratum = d$.stratum, condition = d$.cond),
                          FUN = function(v) c(mean = mean(v), se = stats::sd(v) / sqrt(length(v)), n = length(v)))
  agg <- data.frame(stratum = agg$stratum, condition = agg$condition, agg$x)
  p <- ggplot2::ggplot(agg, ggplot2::aes(x = condition, y = mean, fill = stratum)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - se, ymax = mean + se), width = 0.25) +
    ggplot2::scale_fill_manual(values = colours, guide = "none") +
    ggplot2::facet_wrap(~stratum, nrow = 1) +
    ggplot2::labs(x = "Experiment", y = "mean entropy ± std. err") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
  if (!is.null(pairwise) && nrow(pairwise)) {
    pw <- pairwise
    pw$label <- ifelse(pw$p_adj < alpha, "**", ifelse(also_raw & !is.null(pw$p) & pw$p < alpha, "*", ""))
    pw <- pw[pw$label != "", ]
    if (nrow(pw)) {
      lev <- levels(d$.cond)
      pw$x1 <- match(pw$cond1, lev); pw$x2 <- match(pw$cond2, lev)
      top <- stats::aggregate(agg$mean + agg$se, by = list(stratum = agg$stratum), FUN = max)
      pw$base <- top$x[match(pw$stratum, top$stratum)]
      pw <- pw[order(pw$stratum, abs(pw$x2 - pw$x1)), ]
      pw$y <- stats::ave(pw$base, pw$stratum, FUN = function(b) b * (1 + 0.08 * seq_along(b)))
      p <- p +
        ggplot2::geom_segment(data = pw, ggplot2::aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE) +
        ggplot2::geom_text(data = pw, ggplot2::aes(x = (x1 + x2) / 2, y = y, label = label), vjust = -0.2,
                           size = 5, inherit.aes = FALSE)
    }
  }
  p
}

#' Alluvial plot of enterotype assignments across conditions
#'
#' Reproduces panels B and E of Figures 1 and 2: one stratum bar per condition,
#' flows following each anchor sample. Significance marks are **not** drawn
#' here; put them in the caption or add them from a test result, so that they
#' remain traceable.
#'
#' @inheritParams paired_entropy_test
#' @param label Column with the enterotype assigned in that condition.
#' @param colours Named fill colours.
#' @param condition_levels Order of the conditions.
#' @param ylab Y-axis label.
#' @return A ggplot object.
#' @export
plot_alluvial <- function(data, sample = "sample", condition = "condition", label = "enterotype",
                          colours = enterotype_colours(), condition_levels = NULL,
                          ylab = paste(length(unique(data[[sample]])), "anchor samples")) {
  if (!requireNamespace("ggalluvial", quietly = TRUE)) stop("Package ggalluvial is needed.")
  d <- data
  d$.cond <- if (is.null(condition_levels)) factor(d[[condition]]) else factor(d[[condition]], levels = condition_levels)
  d$.lab <- factor(d[[label]], levels = names(colours))
  ggplot2::ggplot(d, ggplot2::aes(x = .cond, stratum = .lab, alluvium = .data[[sample]], fill = .lab)) +
    ggalluvial::geom_flow(stat = "alluvium") +
    ggalluvial::geom_stratum() +
    ggplot2::scale_fill_manual(values = colours, name = "enterotype") +
    ggplot2::labs(x = "Experiment", y = ylab) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
}
