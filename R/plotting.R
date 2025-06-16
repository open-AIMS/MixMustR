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

#' @export
compare_mixing_proportions <- function(synth_df_d, synth_df_c, mu_tab) {
  rbind(
    dplyr::left_join(
      reshape_ref_data(
        synth_df_d, target = "stream_1_df", order_ref = mu_tab$Group
      ) |>
        data.frame(check.names = FALSE) |>
        dplyr::mutate(N = seq_len(n())) |>
        tidyr::pivot_longer(!N, names_to = "source", values_to = "Tracers"),
      reshape_ref_data(
        synth_df_d, target = "stream_2_df", order_ref = mu_tab$Group
      ) |>
        data.frame(check.names = FALSE) |>
        dplyr::mutate(N = seq_len(n())) |>
        tidyr::pivot_longer(!N, names_to = "source", values_to = "eDNA"),
      by = join_by(N, source)
    ) |>
      dplyr::mutate(`dataset` = "Disagreement (Dataset 2)"),
    dplyr::left_join(
      reshape_ref_data(
        synth_df_c, target = "stream_1_df", order_ref = mu_tab$Group
      ) |>
        data.frame(check.names = FALSE) |>
        dplyr::mutate(N = seq_len(n())) |>
        tidyr::pivot_longer(!N, names_to = "source", values_to = "Tracers"),
      reshape_ref_data(
        synth_df_c, target = "stream_2_df", order_ref = mu_tab$Group
      ) |>
        data.frame(check.names = FALSE) |>
        dplyr::mutate(N = seq_len(n())) |>
        tidyr::pivot_longer(!N, names_to = "source", values_to = "eDNA"),
      by = join_by(N, source)
    ) |>
      dplyr::mutate(`dataset` = "Agreement (Dataset 1)")
  ) |>
    ggplot(data = _) +
      geom_point(
        mapping = aes(x = Tracers, y = eDNA, fill = dataset, shape =  dataset),
        size = 2, alpha = 0.5
      ) +
      geom_abline(slope = 1, linetype = 2) +
      scale_fill_manual(values = c("dodgerblue3", "tomato3")) +
      scale_shape_manual(values = 21:22) +
      labs(x = "From chemical tracers (data stream 1)",
           y = "From eDNA (data stream 2)",
           title = "Simulated mixing proportions",
           fill = "Matrices in:", shape = "Matrices in:") +
      xlim(c(0, 1)) +
      ylim(c(0, 1)) +
      facet_wrap(~source) +
      theme_bw() +
      theme(
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "inside",
        legend.position.inside = c(0.55, 0.2)
      ) +
      guides(
        shape = guide_legend(override.aes = list(size = 3.5)),
        fill = guide_legend(override.aes = list(alpha = 1))
      )
}
