#' @noRd
vary_x <- function(x, yes, ...) {
  if (yes) {
    unsampled_diff <- make_fractionation_mat(x, ...)
    out <- x + unsampled_diff
    out[, "Unsampled"] <- rowSums(unsampled_diff)
    out <- t(apply(out, 1, function(x) {
      x <- linear_rescale(x, c(0, 1))
      abs(x) / sum(abs(x))
    }))
    out[, -which(colnames(out) == "Unsampled")]
  } else {
    x
  }
}

#' @noRd
make_fractionation_mat <- function(x, delta = -0.05) {
  delta <- abs(delta)
  matrix(runif(nrow(x) * ncol(x), -1 * delta, delta), nrow(x), ncol(x))
}

#' @noRd
produce_mix_props <- function(data, sources, seed_1, seed_2,
                              make_unsampled = TRUE, ...) {
  group_n <- data |>
    dplyr::group_by(group) |>
    dplyr::count() |>
    dplyr::ungroup()
  n_groups <- nrow(group_n)
  alphas <- c(1, 1.5, 2, 2.5, 0.5, 2.5) # each value represents a source
  if (length(sources) != length(alphas)) {
    stop("Number of Dirichlet concentration parameters not the same as number",
         " of sources.")
  }
  set.seed(seed_1)
  x <- gtools::rdirichlet(n = n_groups, alpha = alphas) |>
    round(3) |>
    fix_sum_to_one()
  set.seed(seed_2)
  out <- list()
  for (i in seq_len(nrow(group_n))) {
    out[[i]] <- gtools::rdirichlet(n = group_n$n[i], alpha = x[i, ] * 10) |>
      round(3) |>
      fix_sum_to_one() |>
      data.frame() |>
      tibble::as_tibble() |>
      dplyr::mutate(group = .env$group_n$group[i])
  }
  out <- dplyr::bind_rows(out)
  names(out)[seq_along(alphas)] <- sources
  if (make_unsampled) {
    out[, sources] <- vary_x(out[, sources], yes = TRUE, ...)
  }
  out[, order(names(out))] # ensure alphabetical order of names
}

#' @noRd
make_mixture_data <- function(si_df, fa_df, stream_1_props, stream_2_props,
                              truth_stream = 1, template_df, ...) {
  if (truth_stream == 1) {
    mixing_props_df <- as.matrix(stream_1_props)
  } else if (truth_stream == 2) {
    mixing_props_df <- as.matrix(stream_2_props)
  } else {
    stop("truth_stream can only have value 1 or 2.")
  }
  out <- dplyr::bind_rows(
    reshape_isotope_df(si_df), reshape_fattyacids_df(fa_df)
  ) |>
    dplyr::mutate(
      a = dplyr::case_match(marker, "d(13C/12C)" ~ -Inf, .default = 0), b = Inf
    ) |>
    calc_marker_estimate(...) |>
    dplyr::arrange(Group, tracer_family, marker) |> # ensure alphabetical order
    split(f = ~ tracer_family + marker, drop = TRUE) |>
    purrr::map(function(tracer_df, props_df) {
      (t(props_df[, tracer_df$Group]) * tracer_df$estimate) |>
        colSums() |>
        data.frame(check.names = FALSE) |>
        assign_new_names(unique(tracer_df$marker))
    }, props_df = mixing_props_df)
  # check if bind_cols is working as intended
  df_stream_1 <- cbind(template_df, dplyr::bind_cols(out))
  df_stream_2 <- cbind(template_df, stream_2_props)
  stream_1_props <- cbind(template_df, stream_1_props)
  list(df_stream_1 = df_stream_1, df_stream_2 = df_stream_2,
       stream_1_props = stream_1_props)
}

#' @noRd
calc_marker_estimate <- function(x, rand_gen = FALSE, sd_ = NULL, seed = 10) {
  if (rand_gen) {
    set.seed(seed)
    x$sd <- if (!is.null(sd_)) sd_ else x$sd
    dplyr::group_by(x, Group, tracer_family, marker) |>
      dplyr::summarise(
        estimate = trun_na_zr(a = a, b = b, mean = mean, sd = sd)
      ) |>
      dplyr::ungroup()
  } else {
    dplyr::group_by(x, Group, tracer_family, marker) |>
      dplyr::summarise(estimate = mean) |>
      dplyr::ungroup()
  }
}

#' @noRd
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
  list(mus = mu_tab, sigmas = sig_tab)
}

#' @noRd
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

#' @noRd
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

#' @noRd
compare_mixing_proportions <- function(synth_df_d, synth_df_c, mu_tab) {
  rbind(
    dplyr::left_join(
      reshape_ref_data(
        synth_df_d, target = "stream_1_props", order_ref = mu_tab$Group
      ) |>
        data.frame(check.names = FALSE) |>
        dplyr::mutate(N = seq_len(n())) |>
        tidyr::pivot_longer(!N, names_to = "source", values_to = "Tracers"),
      reshape_ref_data(
        synth_df_d, target = "df_stream_2", order_ref = mu_tab$Group
      ) |>
        data.frame(check.names = FALSE) |>
        dplyr::mutate(N = seq_len(n())) |>
        tidyr::pivot_longer(!N, names_to = "source", values_to = "eDNA"),
      by = join_by(N, source)
    ) |>
      dplyr::mutate(`dataset` = "Disagreement (Dataset 2)"),
    dplyr::left_join(
      reshape_ref_data(
        synth_df_c, target = "stream_1_props", order_ref = mu_tab$Group
      ) |>
        data.frame(check.names = FALSE) |>
        dplyr::mutate(N = seq_len(n())) |>
        tidyr::pivot_longer(!N, names_to = "source", values_to = "Tracers"),
      reshape_ref_data(
        synth_df_c, target = "df_stream_2", order_ref = mu_tab$Group
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
