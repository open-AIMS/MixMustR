#' Simulate a hierarchical logistic-normal source composition
#'
#' `simulate_mvn_mixture()` simulates an `N = G * M` composition of `J` source
#' categories using a two-level hierarchical logistic-normal process:
#' group-level centres in unconstrained (multivariate normal) space, and
#' observation-level draws around each group centre, both mapped to the
#' simplex with `softmax()`. This mirrors the composition-generating process
#' used for the MixMustR simulation study (see the accompanying manuscript).
#'
#' @param G Integer. Number of groups.
#' @param M Integer. Number of observations simulated within each group.
#' @param J Integer. Number of source categories (including any pooled
#'   unsampled category), i.e. the dimension of the composition.
#' @param noise Numeric scalar. Within-group variance in unconstrained
#'   (pre-softmax) space; larger values yield more dispersed, less
#'   group-consistent compositions. Defaults to `0.5`.
#' @param seed Integer. Seed used for both the group-centre draws and the
#'   within-group observation draws. Defaults to `1`.
#'
#' @return A list with two elements: `dataset`, a data frame with columns
#'   `group` and `Source_1, ..., Source_J` containing the `N x J` simulated
#'   proportions (rows sum to 1), and `group_means`, the group-level average
#'   of `dataset` for each source.
#'
#' @details
#' Group centres are drawn as `rnorm(G * J, mean = 0, sd = 2)` in
#' unconstrained space. Observations within a group are then drawn from a
#' multivariate normal distribution centred on that group's centre with
#' covariance `diag(J) * noise`, and mapped to the simplex with `softmax()`.
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats rnorm aggregate
#' @keywords internal
#' @noRd
simulate_mvn_mixture <- function(G, M, J, noise = 0.5, seed = 1) {
  softmax <- function(x) {
    exp(x) / sum(exp(x))
  }
  set.seed(seed)
  group_means_unconstrained <- matrix(
    rnorm(G * J, mean = 0, sd = 2), nrow = G, ncol = J
  )
  sigma <- diag(J) * noise
  results_list <- vector(mode = "list", length = G)
  for (g in seq_len(G)) {
    unconstrained_data <- mvrnorm(
      n = M, mu = group_means_unconstrained[g, ], Sigma = sigma
    )
    proportions_data <- t(apply(unconstrained_data, 1, softmax))
    results_list[[g]] <- cbind(data.frame(group = rep(g, M)), proportions_data)
  }
  full_dataset <- do.call(rbind, results_list)
  colnames(full_dataset)[2:(J + 1)] <- paste0("Source_", seq_len(J))
  full_dataset$group <- factor(full_dataset$group)
  group_means_table <- aggregate(. ~ group, data = full_dataset, FUN = mean)
  list(dataset = full_dataset, group_means = group_means_table)
}

#' Rename generic source columns to real source names
#'
#' `rename_sources()` renames the `Source_1, ..., Source_J` columns produced
#' by `simulate_mvn_mixture()` to the source names supplied in
#' `source_names`, matched by the numeric suffix of each `Source_N` column.
#'
#' @param df A data frame containing one or more `Source_N` columns, typically
#'   the `dataset` element returned by `simulate_mvn_mixture()`.
#' @param source_names A character vector of replacement names, one per
#'   `Source_N` column, in the same order (`Source_1` first, and so on).
#'
#' @return `df` with its `Source_N` columns renamed to `source_names`.
#'
#' @importFrom dplyr rename all_of
#' @importFrom stats setNames
#' @keywords internal
#' @noRd
rename_sources <- function(df, source_names) {
  source_cols <- grep("^Source_[0-9]+$", names(df), value = TRUE)
  source_cols <- source_cols[
    order(as.integer(sub("^Source_", "", source_cols)))
  ]
  if (length(source_names) != length(source_cols)) {
    stop(
      "Length of `source_names` must match the number of Source_N columns.\n",
      "Got ", length(source_names), " source names but found ",
      length(source_cols), " Source_N columns: ",
      paste(source_cols, collapse = ", ")
    )
  }
  rename_map <- setNames(source_cols, source_names)
  rename(df, all_of(rename_map))
}

#' Strip the pooled-unsampled column from simulated proportion data
#'
#' `rm_unsampled()` removes the `Unsampled` proportion column from the
#' `df_stream_2` and `stream_1_props` elements of a simulated mixture-data
#' list, mimicking the realistic situation in which users do not supply the
#' unsampled-source proportion directly; \code{\link{mixmustr_wrangle_input}}
#' reconstructs it as `1 - rowSums(sampled proportions)`.
#'
#' @param x A list with `df_stream_1`, `df_stream_2`, and `stream_1_props`
#'   elements, typically produced by `make_mixture_data()`.
#'
#' @return `x`, with the `Unsampled` column removed from `df_stream_2` and
#'   `stream_1_props`.
#'
#' @importFrom dplyr select
#' @keywords internal
#' @noRd
rm_unsampled <- function(x) {
  x$df_stream_2 <- select(x$df_stream_2, -"Unsampled")
  x$stream_1_props <- select(x$stream_1_props, -"Unsampled")
  x
}

#' Construct an out-of-hull unsampled-source tracer signature
#'
#' `make_unsampled_signature()` builds a synthetic tracer signature for a
#' pooled unsampled source, positioned `delta` scaled units away from the
#' sampled-source centroid along the leading principal-component direction of
#' the standardised sampled-source signatures. This reproduces the simulation
#' design used to test MixMustR's ability to recover a source whose signature
#' does not overlap the sampled-source multivariate space (see the
#' accompanying manuscript).
#'
#' @param mus_tab A data frame of sampled-source mean tracer signatures, with
#'   a `source` column followed by one numeric column per tracer, e.g.
#'   `tracer_parameters$mus`.
#' @param delta Numeric scalar. Number of scaled units to displace the
#'   unsampled signature from the sampled-source centroid along the first
#'   principal component. Defaults to `3`.
#'
#' @return A one-row tibble with a `source` column (`"Unsampled"`) followed by
#'   the same tracer columns as `mus_tab`, on the original (unstandardised)
#'   tracer scale.
#'
#' @importFrom tibble column_to_rownames as_tibble_row
#' @importFrom dplyr mutate
#' @importFrom stats prcomp
#' @keywords internal
#' @noRd
make_unsampled_signature <- function(mus_tab, delta = 3) {
  s <- column_to_rownames(mus_tab, "source") |> as.matrix()
  s_scaled <- scale(s)
  center_tracer <- attr(s_scaled, "scaled:center")
  scale_tracer <- attr(s_scaled, "scaled:scale")
  # Direction of strongest source separation in scaled tracer space.
  pc <- prcomp(s_scaled, center = FALSE, scale. = FALSE)
  u_scaled <- colMeans(s_scaled) + delta * pc$rotation[, 1]
  u_original <- u_scaled * scale_tracer + center_tracer
  as_tibble_row(as.list(u_original)) |>
    mutate(source = "Unsampled", .before = 1)
}

#' Append a simulated unsampled-source row to the tracer signature tables
#'
#' `augment_bcs_with_unsampled()` combines `bcs_si` and `bcs_fa`, homogenises
#' all standard-deviation and sample-size columns to the values used in the
#' MixMustR simulation study (`1` and `10`, respectively), and appends a
#' pooled-unsampled row built by `make_unsampled_signature()`. The result is
#' split back into `bcs_si`- and `bcs_fa`-shaped tables ready to be passed to
#' `make_mixture_data()`.
#'
#' @param bcs_si A data frame of sampled-source stable-isotope signatures,
#'   structured as \code{\link{bcs_si}}.
#' @param bcs_fa A data frame of sampled-source fatty-acid signatures,
#'   structured as \code{\link{bcs_fa}}.
#' @param mus_tab A data frame of sampled-source mean tracer signatures used
#'   to position the unsampled signature, e.g. `tracer_parameters$mus`.
#' @param delta Numeric scalar passed to `make_unsampled_signature()`.
#'   Defaults to `3`.
#'
#' @return A list with two elements, `bcs_si` and `bcs_fa`, matching the
#'   structure of the package's `bcs_si`/`bcs_fa` datasets but with an
#'   additional `"Unsampled"` row and homogenised SD/n columns throughout.
#'
#' @details
#' Homogenising the process SD (`1`) and sample size (`10`) across all
#' sources and tracers isolates the effect of the other simulation-design
#' factors (number of tracers, number of sources, stream agreement, etc.) on
#' recovery of the mixing proportions, matching the factorial simulation
#' experiment described in the accompanying manuscript.
#'
#' @importFrom dplyr left_join mutate across bind_rows select relocate
#' @importFrom tidyselect contains
#' @keywords internal
#' @noRd
augment_bcs_with_unsampled <- function(bcs_si, bcs_fa, mus_tab, delta = 3) {
  bcs <- left_join(
    bcs_si, select(bcs_fa, !c("Study", "Taxa")), by = "source"
  ) |>
    mutate(across(contains("SD"), ~1), across(contains("(n)"), ~10L))
  unsmpd <- make_unsampled_signature(mus_tab, delta = delta) |>
    mutate(Study = "Simulated")
  # SD/n are set directly to the homogenised values above; deriving them from
  # `bcs` (as in the original simulation script) would return the same
  # constants because `bcs` was already homogenised.
  for (trc in setdiff(names(mus_tab), "source")) {
    unsmpd[[paste0(trc, " (SD)")]] <- 1
    unsmpd[[paste0(trc, " (n)")]] <- 10L
  }
  bcs <- bind_rows(bcs, unsmpd)
  list(
    bcs_si = select(bcs, "source":"Study"),
    bcs_fa = select(bcs, "source", "24:0":"20:5w3 (n)", "Study") |>
      mutate(Taxa = NA) |>
      relocate("Taxa", .after = "source")
  )
}

#' @importFrom dplyr bind_rows bind_cols mutate recode_values arrange
#' @importFrom purrr map
#' @importFrom rlang .data
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
  out <- bind_rows(
    reshape_isotope_df(si_df), reshape_fattyacids_df(fa_df)
  ) |>
    mutate(
      a = recode_values(.data$tracer, "d(13C/12C)" ~ -Inf, default = 0), b = Inf
    ) |>
    calc_tracer_estimate(...) |>
    arrange(.data$source, .data$tracer_family, .data$tracer) |> # ensure alphabetical order
    split(f = ~ tracer_family + tracer, drop = TRUE) |>
    map(function(tracer_df, props_df) {
      (t(props_df[, tracer_df$source]) * tracer_df$estimate) |>
        colSums() |>
        data.frame(check.names = FALSE) |>
        assign_new_names(unique(tracer_df$tracer))
    }, props_df = mixing_props_df)
  # check if bind_cols is working as intended
  df_stream_1 <- cbind(template_df, bind_cols(out))
  df_stream_2 <- cbind(template_df, stream_2_props)
  stream_1_props <- cbind(template_df, stream_1_props)
  list(df_stream_1 = df_stream_1, df_stream_2 = df_stream_2,
       stream_1_props = stream_1_props)
}

#' @importFrom dplyr group_by summarise ungroup
#' @importFrom rlang .data
#' @noRd
calc_tracer_estimate <- function(x, rand_gen = FALSE, sd_ = NULL, seed = 10) {
  if (rand_gen) {
    set.seed(seed)
    x$sd <- if (!is.null(sd_)) sd_ else x$sd
    group_by(x, .data$source, .data$tracer_family, .data$tracer) |>
      summarise(
        estimate = trun_na_zr(a = .data$a, b = .data$b, mean = .data$mean,
                              sd = .data$sd)
      ) |>
      ungroup()
  } else {
    group_by(x, .data$source, .data$tracer_family, .data$tracer) |>
      summarise(estimate = .data$mean) |>
      ungroup()
  }
}

#' @importFrom dplyr left_join select mutate across
#' @importFrom tidyr pivot_wider
#' @importFrom tidyselect where
#' @noRd
wrangle_tracer_pars <- function(raw_data_si, raw_data_fa) {
  mu_tab <- left_join(
    reshape_isotope_df(raw_data_si) |>
      select(!c("Study", "sd", "n", "tracer_family")) |>
      pivot_wider(names_from = "tracer", values_from = "mean"),
    reshape_fattyacids_df(raw_data_fa) |>
      select(!c("Study", "sd", "n", "tracer_family")) |>
      pivot_wider(names_from = "tracer", values_from = "mean"),
    by = "source"
  )
  sig_tab <- left_join(
    reshape_isotope_df(raw_data_si) |>
      select(!c("Study", "mean", "n", "tracer_family")) |>
      pivot_wider(names_from = "tracer", values_from = "sd"),
    reshape_fattyacids_df(raw_data_fa) |>
      select(!c("Study", "mean", "n", "tracer_family")) |>
      pivot_wider(names_from = "tracer", values_from = "sd"),
    by = "source"
  )
  n_tab <- left_join(
    reshape_isotope_df(raw_data_si) |>
      select(!c("Study", "mean", "sd", "tracer_family")) |>
      pivot_wider(names_from = "tracer", values_from = "n"),
    reshape_fattyacids_df(raw_data_fa) |>
      select(!c("Study", "mean", "sd", "tracer_family")) |>
      pivot_wider(names_from = "tracer", values_from = "n"),
    by = "source"
  )
  list(mus = mu_tab, sigmas = sig_tab, ns = n_tab)
}

#' @importFrom dplyr select left_join mutate case_match
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect starts_with ends_with
#' @importFrom rlang .data
#' @noRd
reshape_isotope_df <- function(x) {
  x |>
    select(-ends_with("(SD)"), -ends_with("(n)")) |>
    pivot_longer(
      cols = starts_with("d("), values_to = "mean", names_to = "tracer"
    ) |>
    left_join(
      x |>
        select("source", ends_with("(SD)"), "Study") |>
        pivot_longer(
          cols = ends_with("(SD)"), values_to = "sd", names_to = "tracer"
        ) |>
        mutate(tracer = gsub(" (SD)", "", .data$tracer, fixed = TRUE)),
      by = c("source", "Study", "tracer")
    ) |>
    left_join(
      x |>
        select("source", ends_with("(n)"), "Study") |>
        pivot_longer(
          cols = ends_with("(n)"), values_to = "n", names_to = "tracer"
        ) |>
        mutate(tracer = gsub(" (n)", "", .data$tracer, fixed = TRUE)),
      by = c("source", "Study", "tracer")
    ) |>
    mutate(tracer_family = "si")
}

#' @importFrom dplyr select left_join mutate
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect ends_with
#' @importFrom rlang .data
#' @noRd
reshape_fattyacids_df <- function(x) {
  x |>
    select(-"Taxa", -ends_with("(SD)"), -ends_with("(n)")) |>
    pivot_longer(
      cols = "24:0":"20:5w3", values_to = "mean",
      names_to = "tracer"
    ) |>
    left_join(
      x |>
        select("source", ends_with("(SD)"), "Study") |>
        pivot_longer(
          cols = ends_with("(SD)"), values_to = "sd", names_to = "tracer"
        ) |>
        mutate(tracer = gsub(" (SD)", "", .data$tracer, fixed = TRUE)),
      by = c("source", "Study", "tracer")
    ) |>
    left_join(
      x |>
        select("source", ends_with("(n)"), "Study") |>
        pivot_longer(
          cols = ends_with("(n)"), values_to = "n", names_to = "tracer"
        ) |>
        mutate(tracer = gsub(" (n)", "", .data$tracer, fixed = TRUE)),
      by = c("source", "Study", "tracer")
    ) |>
    mutate(tracer_family = "fa")
}

#' Compare Mixing Proportions Between Data Streams
#'
#' This function generates a plot to compare the mixing proportions between two data streams 
#' (e.g., chemical tracers and eDNA) for synthetic datasets with agreement and disagreement.
#'
#' @param synth_list_d A list representing the synthetic dataset with divergent mixing proportions.
#' @param synth_list_c A list representing the synthetic dataset with convergent mixing proportions.
#' @param mu_tab A data frame containing the mean tracer
#' signatures for each source. The first column should be named `source`, and
#' the remaining columns should contain numeric tracer values.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object visualizing the comparison of mixing proportions between the two data streams.
#'
#' @details
#' The function compares the mixing proportions from two data streams (`stream_1_props` and `df_stream_2`) 
#' for both divergent and convergent synthetic datasets. It reshapes and aligns the data based on the 
#' reference order provided in `mu_tab`, and then creates a scatter plot with points representing 
#' the proportions from each data stream. The plot includes facets for each source and highlights 
#' agreement or disagreement between the datasets.
#'
#' @importFrom dplyr left_join mutate n
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot geom_point aes geom_abline scale_fill_manual
#' @importFrom ggplot2 scale_shape_manual labs xlim ylim facet_wrap theme_bw
#' @importFrom ggplot2 theme element_text guides guide_legend
#' @importFrom rlang .data
#'
#' @examples
#' library(MixMustR)
#' compare_mixing_proportions(
#'   synth_list_d = synthetic_df_divergent,
#'   synth_list_c = synthetic_df_convergent,
#'   mu_tab = tracer_parameters$mus
#' )
#'
#' @seealso
#' \code{\link{synthetic_df_divergent}}, \code{\link{synthetic_df_convergent}}
#'
#' @export
compare_mixing_proportions <- function(synth_list_d, synth_list_c, mu_tab) {
  rbind(
    left_join(
      reshape_ref_data(
        synth_list_d, target = "stream_1_props", order_ref = mu_tab$source
      ) |>
        data.frame(check.names = FALSE) |>
        mutate(N = seq_len(n())) |>
        pivot_longer(!"N", names_to = "source", values_to = "Tracers"),
      reshape_ref_data(
        synth_list_d, target = "df_stream_2", order_ref = mu_tab$source
      ) |>
        data.frame(check.names = FALSE) |>
        mutate(N = seq_len(n())) |>
        pivot_longer(!"N", names_to = "source", values_to = "eDNA"),
      by = c("N", "source")
    ) |>
      mutate(`dataset` = "Disagreement (Dataset 2)"),
    left_join(
      reshape_ref_data(
        synth_list_c, target = "stream_1_props", order_ref = mu_tab$source
      ) |>
        data.frame(check.names = FALSE) |>
        mutate(N = seq_len(n())) |>
        pivot_longer(!"N", names_to = "source", values_to = "Tracers"),
      reshape_ref_data(
        synth_list_c, target = "df_stream_2", order_ref = mu_tab$source
      ) |>
        data.frame(check.names = FALSE) |>
        mutate(N = seq_len(n())) |>
        pivot_longer(!"N", names_to = "source", values_to = "eDNA"),
      by = c("N", "source")
    ) |>
      mutate(`dataset` = "Agreement (Dataset 1)")
  ) |>
    ggplot(data = _) +
      geom_point(
        mapping = aes(x = .data$Tracers, y = .data$eDNA, fill = .data$dataset,
                      shape = .data$dataset), size = 2, alpha = 0.5
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
      facet_wrap(~ .data$source) +
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
