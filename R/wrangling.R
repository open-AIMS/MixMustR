#' @export
wrangle_tracer_pars <- function(raw_data_si, raw_data_fa) {
  mu_tab <- dplyr::left_join(
    reshape_isotope_df(raw_data_si) |>
      dplyr::select(!c(Study, sd, tracer_family)) |>
      tidyr::pivot_wider(names_from = "marker", values_from = "mean"),
    reshape_fattyacids_df(raw_data_fa) |>
      dplyr::select(!c(Study, sd, tracer_family)) |>
      tidyr::pivot_wider(names_from = "marker", values_from = "mean"),
    by = dplyr::join_by(Group)
  )
  sig_tab <- dplyr::left_join(
    reshape_isotope_df(raw_data_si) |>
      dplyr::select(!c(Study, mean, tracer_family)) |>
      tidyr::pivot_wider(names_from = "marker", values_from = "sd"),
    reshape_fattyacids_df(raw_data_fa) |>
      dplyr::select(!c(Study, mean, tracer_family)) |>
      tidyr::pivot_wider(names_from = "marker", values_from = "sd"),
    by = dplyr::join_by(Group)
  ) |>
    dplyr::mutate(across(where(is.numeric), function(x) {
      x[x == 0] <- 0.01 # NB: check this. Model can't deal with 0 variance
      x
    }))
  list(mu_tab = mu_tab, sig_tab = sig_tab)
}

#' @importFrom dplyr %>% mutate
#' @export
reshape_ref_data <- function(x, target = "stream_1_df", order_ref) {
  (x[[target]][, order_ref]) %>%
    mutate(Unsampled = 1 - rowSums(.)) |>
    as.matrix() |>
    abs()
}

#' @export
reshape_isotope_df <- function(x) {
  x |>
    dplyr::select(-ends_with("sd")) |>
    tidyr::pivot_longer(
      cols = starts_with("d("), values_to = "mean", names_to = "marker"
    ) |>
    dplyr::left_join(
      x |>
        dplyr::select(-starts_with("d(")) |>
        tidyr::pivot_longer(
          cols = ends_with("sd"), values_to = "sd", names_to = "marker"
        ) |>
        dplyr::mutate(
          marker = dplyr::case_match(marker,
            "d13C sd" ~ "d(13C/12C)", "d15N sd" ~ "d(15N/14N)",
            .default = NA
          )
        ),
      by = dplyr::join_by(Group, Study, marker)
    ) |>
    dplyr::mutate(tracer_family = "si")
}

#' @export
reshape_fattyacids_df <- function(x) {
  x |>
    dplyr::select(-Taxa, -ends_with("(SD)")) |>
    tidyr::pivot_longer(
      cols = `24:0`:`20:5w3`, values_to = "mean", names_to = "marker"
    ) |>
    dplyr::left_join(
      x |>
        dplyr::select(Group, -Taxa, ends_with("(SD)"), Study) |>
        tidyr::pivot_longer(
          cols = ends_with("(SD)"), values_to = "sd", names_to = "marker"
        ) |>
        dplyr::mutate(marker = gsub(" (SD)", "", marker, fixed = TRUE)),
      by = dplyr::join_by(Group, Study, marker)
    ) |>
    dplyr::mutate(tracer_family = "fa")
}
