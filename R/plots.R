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
    ggplot2::labs(x = "Experiment", y = "mean entropy \u00b1 std. err") +
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

#' Assemble the B, C, E, F panels of Figures 1 and 2, leaving A and D empty
#'
#' Lays out two rows (alluvial + entropy panels) with \pkg{patchwork} at the
#' proportions of the published page, with an empty slot at the left of each
#' row where the nVennR "flower" schematic (panels A and D) is added in
#' Inkscape afterwards. Panel letters are set with `labs(tag =)` on each plot.
#'
#' @param alluvial_top,entropy_top,alluvial_bottom,entropy_bottom ggplot objects.
#' @param letters Panel letters for the four plots, in that order.
#' @param widths Relative widths of the flower slot, the alluvial and the entropy panel.
#' @return A patchwork object; save with `ggplot2::ggsave(width = 13.6, height = 11)`
#'   to match the 978 x 795 pt page of the published figure.
#' @export
plot_figure_panels <- function(alluvial_top, entropy_top, alluvial_bottom, entropy_bottom,
                               letters = c("B", "C", "E", "F"), widths = c(1.1, 1.1, 2.6)) {
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("Package patchwork is needed.")
  tag <- function(p, l) p + ggplot2::labs(tag = l) + ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 16))
  row <- function(a, e) patchwork::plot_spacer() + a + e + patchwork::plot_layout(widths = widths)
  top <- row(tag(alluvial_top, letters[1]), tag(entropy_top, letters[2]))
  bottom <- row(tag(alluvial_bottom, letters[3]), tag(entropy_bottom, letters[4]))
  patchwork::wrap_plots(top, bottom, ncol = 1)
}

#' Robustness against cohort fraction with a confidence ribbon
#'
#' Reproduces Figures 3 and 4: mean robustness per group (enterotype or pair)
#' along the cohort fraction, with a ribbon between `lower` and `upper`.
#'
#' @param ci Data frame with `fraction`, `group`, `mean`, `lower`, `upper` and
#'   optionally `facet`.
#' @param colours Named colours per group; groups absent from it get defaults.
#' @param ylab,xlab Axis labels.
#' @return A ggplot object.
#' @export
plot_robustness <- function(ci, colours = enterotype_colours(), ylab = "mean robustness",
                            xlab = "fraction of the cohort") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package ggplot2 is needed.")
  groups <- unique(ci$group)
  cols <- colours[intersect(names(colours), groups)]
  missing <- setdiff(groups, names(cols))
  if (length(missing)) cols <- c(cols, stats::setNames(grDevices::hcl.colors(length(missing), "Dark 3"), missing))
  p <- ggplot2::ggplot(ci, ggplot2::aes(x = fraction, y = mean, colour = group, fill = group)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.25, colour = NA) +
    ggplot2::geom_line() + ggplot2::geom_point(size = 1) +
    ggplot2::scale_colour_manual(values = cols) + ggplot2::scale_fill_manual(values = cols) +
    ggplot2::geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey40") +
    ggplot2::scale_y_continuous(limits = c(NA, 1)) +
    ggplot2::labs(x = xlab, y = ylab, colour = NULL, fill = NULL) +
    ggplot2::theme_bw()
  if ("facet" %in% names(ci)) p <- p + ggplot2::facet_wrap(~facet, nrow = 1)
  p
}
