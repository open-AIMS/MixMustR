#' @importFrom ggplot2 ggplot geom_errorbarh aes geom_point scale_shape_manual
#' @importFrom ggplot2 geom_abline labs xlim ylim facet_wrap theme_bw
#' @importFrom ggsci scale_fill_jco
#' @importFrom rlang .data
#' @export
plot_multiple_faceted_scatter_avg <- function (data, ...) {
  ggplot(data = data) +
    geom_errorbarh(
      mapping = aes(y = .data$Observed, xmin = .data$ymin, xmax = .data$ymax),
      height = 0, linewidth = 0.2, colour = "grey60", alpha = 0.8
    ) +
    geom_point(
      mapping = aes(
        x = .data$Predicted, y = .data$Observed, shape = .data$`Variant:`,
        fill = .data$`Variant:`
      ), size = 2, alpha = 0.8
    ) +
    scale_shape_manual(values = 21:24) +
    scale_fill_jco() +
    geom_abline(slope = 1, linetype = 2) +
    labs(x = "Expected (mean posterior fit +/- 95% HDI)",
         title = "Mixing proportions", ...) +
    xlim(c(0, 1)) +
    ylim(c(0, 1)) +
    facet_wrap(~ .data$source) +
    theme_bw()
}
