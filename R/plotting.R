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
#' library(MixMustR)
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

#' Visually inspect uncertainty for input mixing proportion
#'
#' This function plots the expected uncertainty around each input mixing
#' proportion given its uncertainty on the log scale.
#'
#' @param pi A numeric vector containing mixing proportions
#' (needs to sum to 1). If named, the function uses the names for plotting
#' purposes.
#' @param ln_sigma_rho A numeric vector containing the uncertainty 
#' around `pi` on the log scale.
#' @param iter Integer. The number of iterations to run the simulation for.
#' Defaults to 10,000.
#' @param seed Integer. Seed for simulation reproducibility. Defaults to 10.
#' Defaults to 10.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object representing the faceted
# scatter plot.
#'
#' @details
#' We recognise that choosing the appropriate value for the uncertainty around
#' data stream 2, \eqn{\sigma_{\text{ln}\rho}}, is not trivial. Therefore, the
#' MixMustR package offers a simple tool which plots the expected variability
#' in mixing proportion as a function of an input vector of mixing proportions
#' from data stream 2 (one from each source), and their respective degrees of
#' confidence or measurement error. Although the function does not inform the
#' user what the expected uncertainty for each source and observation should be
#' (what has to come from expect knowledge), it does allow the user to fine
#' tune \eqn{\sigma_{\text{ln}\rho}} until the appropriate variability in input
#' mixing proportions is achieved.
#'
#' @importFrom ggplot2 ggplot geom_density aes labs scale_x_continuous
#' @importFrom ggplot2 facet_wrap theme_bw
#' @importFrom tidyselect everything
#' @importFrom tidyr pivot_longer
#' @importFrom stats rnorm
#' @importFrom rlang .data
#'
#' @examples
#' x = c(0.06, 0.6, 0.04, 0.3)
#' ln_sigma_x <- rep(1, length(x))
#' MixMustR::evaluate_uncertainty(x, ln_sigma_x)$plot
#' 
#' @export
evaluate_uncertainty <- function(pi, ln_sigma_rho, iter = 1e4, seed = 10) {
  softmax <- function(x) {
    exp(x) / sum(exp(x))
  }
  if (length(pi) != length(ln_sigma_rho)) {
    stop("Arguments `pi` and `ln_sigma_rho` must have the same length.")
  }
  if (!is.numeric(pi) | !is.numeric(ln_sigma_rho)) {
    stop("Arguments `pi` and `ln_sigma_rho` must be numeric vectors.")
  }
  ln_x <- log(pi)
  out <- matrix(0, iter, length(ln_x))
  set.seed(seed)
  for (i in seq_len(iter)) {
    out[i, ] <- rnorm(length(ln_x), ln_x, ln_sigma_rho) |>
      softmax()
  }
  if (is.null(names(pi))) {
    colnames(out) <- paste0("Source ", seq_len(ncol(out)))
  } else {
    colnames(out) <- names(pi)
  }
  out <- data.frame(out, check.names = FALSE) |>
    pivot_longer(everything(), names_to = "Sources", values_to = "Values")
  plot_out <- ggplot(data = out) +
    geom_density(mapping = aes(x = .data$Values), trim = TRUE, adjust = 2,
                 fill = "grey30", alpha = 0.8) +
    labs(x = "Mixing proportions", y = "Density") +
    scale_x_continuous(limits = c(0, 1)) +
    facet_wrap(~.data$Sources, scales = "free_y") +
    theme_bw()
  list(data = out, plot = plot_out)
}
