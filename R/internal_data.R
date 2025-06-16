#' @export
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

#' @export
make_fractionation_mat <- function(x, delta = -0.05) {
  delta <- abs(delta)
  matrix(runif(nrow(x) * ncol(x), -1 * delta, delta), nrow(x), ncol(x))
}

#' @export
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

#' @export
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
  sim_data <- cbind(template_df, dplyr::bind_cols(out))
  stream_1_df <- cbind(template_df, stream_1_props)
  stream_2_df <- cbind(template_df, stream_2_props)
  list(sim_data = sim_data, stream_1_df = stream_1_df,
       stream_2_df = stream_2_df)
}

#' @export
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
  list(mus = mu_tab, sigmas = sig_tab)
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
