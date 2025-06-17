#' @importFrom rstan extract
#' @importFrom ggdist mean_hdci
#' @importFrom dplyr bind_rows select mutate rename left_join n
#' @importFrom tidyr pivot_longer
#' @importFrom rlang .env .data
#' @export
make_post_prop_long <- function(modfit, mu_tab, synth_df, n, ...) {
  array_of_draws <- extract(modfit)$p
  ordered_sources <- c(mu_tab$Group, "Unsampled")
  names(ordered_sources) <- as.character(seq_along(ordered_sources))
  pred_mean_hdci <- apply(array_of_draws, 3, function(x) {
    apply(x, 2, mean_hdci) |>
      bind_rows(.id = "N") |>
      select(!(.data$`.width`:.data$`.interval`))
  }) |>
    bind_rows(.id = "source") |>
    mutate(N = as.integer(.data$N),
           source = .env$ordered_sources[.data$source]) |>
    rename(Predicted = .data$y)
  left_join(
    pred_mean_hdci,
    reshape_ref_data(synth_df, order_ref = mu_tab$Group, ...) |>
      data.frame(check.names = FALSE) |>
      mutate(N = seq_len(n())) |>
      pivot_longer(
        cols = !.data$N, names_to = "source", values_to = "Observed"
      ), by = c("N", "source")
  ) |>
    mutate(`Variant:` = as.character(.env$n))
}

#' @importFrom ggdist mean_hdci
#' @importFrom dplyr bind_rows select mutate left_join filter n
#' @importFrom dplyr across
#' @importFrom tidyr pivot_longer separate_wider_regex
#' @importFrom rlang .env .data
#' @export
make_jags_post_prop_long <- function(modfit, mu_tab, synth_df, n, ...) {
  ordered_sources <- sort(mu_tab$Group)
  names(ordered_sources) <- as.character(seq_along(ordered_sources))
  left_join(
    modfit$BUGSoutput$sims.matrix |>
      apply(2, mean_hdci) |>
      bind_rows(.id = "Variable") |>
      filter(grepl("p\\.ind\\[", .data$Variable)) |>
      separate_wider_regex(
        .data$Variable,
        c("p\\.ind\\[", N = "\\d+", "\\,", source = "\\d+", "\\]")
      ) |>
      mutate(across(c(.data$N, .data$source), as.numeric)) |>
      mutate(source = .env$ordered_sources[.data$source]) |>
      select(
        .data$source, .data$N, Predicted = .data$y, .data$ymin, .data$ymax
      ),
    reshape_ref_data(synth_df, order_ref = mu_tab$Group, ...) |>
      data.frame(check.names = FALSE) |>
      select(-.data$Unsampled) |>
      mutate(N = seq_len(n())) |>
      pivot_longer(
        cols = !.data$N, names_to = "source", values_to = "Observed"
      ), by = c("N", "source")
  ) |>
    mutate(`Variant:` = as.character(.env$n))
}
