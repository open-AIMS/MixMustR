#' Plot Multiple Faceted Scatter with Averages
#'
#' This function creates a faceted scatter plot to compare observed and predicted values, including error bars and visual grouping by variants. It is designed for visualising mixing proportions or similar data with multiple sources. It's mostly used to plot outputs of make_post_prop_long.
#'
#' @param data A data frame containing the data to be plotted. It must include the following columns:
#'   \describe{
#'     \item{Observed}{Numeric. The observed values for each source.}
#'     \item{Predicted}{Numeric. The predicted values for each source.}
#'     \item{ymin}{Numeric. The lower bound of the error bar for the predicted values.}
#'     \item{ymax}{Numeric. The upper bound of the error bar for the predicted values.}
#'     \item{Variant:}{Factor or character. The grouping variable for different variants.}
#'     \item{source}{Factor or character. The source variable used for faceting.}
#'   }
#' @param ... Additional arguments passed to the `labs` function for customizing plot labels.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object representing the faceted
# scatter plot.
#'
#' @details
#' The function generates a scatter plot with error bars for predicted values, points representing 
#' observed vs. predicted values, and a diagonal reference line (slope = 1) for visual comparison. 
#' The plot is faceted by the `source` variable and grouped by the `Variant:` variable, with 
#' customizable shapes and colors for each variant.
#'
#' The `scale_fill_jco` function from the `ggsci` package is used to apply a color palette, 
#' and the `scale_shape_manual` function is used to define point shapes.
#'
#' @importFrom ggplot2 ggplot geom_errorbarh aes geom_point scale_shape_manual
#' @importFrom ggplot2 geom_abline labs xlim ylim facet_wrap theme_bw
#' @importFrom ggsci scale_fill_jco
#' @importFrom rlang .data
#'
#' @examples
#' data <- data.frame(
#'   Observed = c(0.2, 0.4, 0.6, 0.8),
#'   Predicted = c(0.25, 0.35, 0.65, 0.75),
#'   ymin = c(0.2, 0.3, 0.6, 0.7),
#'   ymax = c(0.3, 0.5, 0.7, 0.9),
#'   `Variant:` = c("A", "B", "A", "B"),
#'   source = c("Source1", "Source1", "Source2", "Source2"),
#'   check.names = FALSE
#' )
#' plot_multiple_faceted_scatter_avg(data)
#'
#' @seealso
#'   \code{\link{make_post_prop_long}}
#' 
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
