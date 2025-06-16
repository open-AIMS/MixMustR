#' @export
make_post_prop_long <- function(modfit, mu_tab, synth_df, n, ...) {
  array_of_draws <- rstan::extract(modfit)$p
  ordered_sources <- c(mu_tab$Group, "Unsampled")
  names(ordered_sources) <- as.character(seq_along(ordered_sources))
  pred_mean_hdci <- apply(array_of_draws, 3, function(x) {
    apply(x, 2, ggdist::mean_hdci) |>
      dplyr::bind_rows(.id = "N") |>
      dplyr::select(!(`.width`:`.interval`))
  }) |>
    dplyr::bind_rows(.id = "source") |>
    dplyr::mutate(N = as.integer(N), source = .env$ordered_sources[source]) |>
    dplyr::rename(Predicted = y)
  dplyr::left_join(
    pred_mean_hdci,
    reshape_ref_data(synth_df, order_ref = mu_tab$Group, ...) |>
      data.frame(check.names = FALSE) |>
      dplyr::mutate(N = seq_len(n())) |>
      tidyr::pivot_longer(
        cols = !N, names_to = "source", values_to = "Observed"
      ), by = dplyr::join_by(N, source)
  ) |>
    dplyr::mutate(`Variant:` = as.character(.env$n))
}

#' @export
make_jags_post_prop_long <- function(modfit, mu_tab, synth_df, n, ...) {
  ordered_sources <- sort(mu_tab$Group)
  names(ordered_sources) <- as.character(seq_along(ordered_sources))
  dplyr::left_join(
    modfit$BUGSoutput$sims.matrix |>
      apply(2, ggdist::mean_hdci) |>
      dplyr::bind_rows(.id = "Variable") |>
      dplyr::filter(grepl("p\\.ind\\[", Variable)) |>
      tidyr::separate_wider_regex(
        Variable, c("p\\.ind\\[", N = "\\d+", "\\,", source = "\\d+", "\\]")
      ) |>
      dplyr::mutate(dplyr::across(c(N, source), as.numeric)) |>
      dplyr::mutate(source = .env$ordered_sources[source]) |>
      dplyr::select(source, N, Predicted = y, ymin, ymax),
    reshape_ref_data(synth_df, order_ref = mu_tab$Group, ...) |>
      data.frame(check.names = FALSE) |>
      dplyr::select(-Unsampled) |>
      dplyr::mutate(N = seq_len(n())) |>
      tidyr::pivot_longer(
        cols = !N, names_to = "source", values_to = "Observed"
      ), by = dplyr::join_by(N, source)
  ) |>
    dplyr::mutate(`Variant:` = as.character(.env$n))
}
