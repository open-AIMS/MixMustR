#' @export
plot_multiple_faceted_scatter_avg <- function (data, ...) {
  ggplot(data = data) +
    geom_errorbarh(
      mapping = aes(y = Observed, xmin = ymin, xmax = ymax),
      height = 0, linewidth = 0.2, colour = "grey60", alpha = 0.8
    ) +
    geom_point(
      mapping = aes(
        x = Predicted, y = Observed, shape = `Variant:`, fill = `Variant:`
      ), size = 2, alpha = 0.8
    ) +
    scale_shape_manual(values = 21:24) +
    ggsci::scale_fill_jco() +
    geom_abline(slope = 1, linetype = 2) +
    labs(x = "Expected (mean posterior fit +/- 95% HDI)",
         title = "Mixing proportions", ...) +
    xlim(c(0, 1)) +
    ylim(c(0, 1)) +
    facet_wrap(~source) + 
    theme_bw()
}
