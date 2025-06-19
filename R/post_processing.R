#' Create Long-Format Posterior Proportions from mixmustr model fits
#'
#' This function processes posterior draws from a Stan model to compute mean and highest density credible intervals (HDI) for mixing proportions. It reshapes the data into a long format for comparison with observed values.
#'
#' @inheritParams mixmustr_wrangle_input
#' @param modfit A fitted Stan model object containing posterior draws.
#' Typically produced by function \code{\link{run_mixmustr_models}}.
#' @param mu_tab A data frame containing the mean tracer
#' signatures for each source. The first column should be named `source`, and
#' the remaining columns should contain numeric tracer values.
#' @param n A numeric or character value representing the variant identifier for
#' the model.
#' @param ... Additional arguments passed to internal function.
#'
#' @return A data frame in long format with the following columns:
#' \describe{
#'   \item{N}{Integer. The observation index.}
#'   \item{source}{Character. The source identifier (e.g., sampled sources or "Unsampled").}
#'   \item{Predicted}{Numeric. The mean posterior predicted value for the mixing proportion.}
#'   \item{ymin}{Numeric. The lower bound of the 95% HDI for the predicted value.}
#'   \item{ymax}{Numeric. The upper bound of the 95% HDI for the predicted value.}
#'   \item{Observed}{Numeric. The observed mixing proportion for the corresponding source.}
#'   \item{Variant:}{Character. The variant identifier for the dataset.}
#' }
#'
#' @details
#' The function extracts posterior draws of mixing proportions from the fitted Stan model (`modfit`) and computes the mean and HDI for each source. It then aligns these predictions with observed values from the input data streams (`data_streams_list`) and reshapes the data into a long format for further analysis or visualization.
#'
#' The `Variant:` column is used to distinguish between different models.
#'
#' @importFrom rstan extract
#' @importFrom ggdist mean_hdci
#' @importFrom dplyr bind_rows select mutate rename left_join n
#' @importFrom tidyr pivot_longer
#' @importFrom rlang .env .data
#'
#' @examples
#' \dontrun{
#' library(mixmustr)
#' # mixmustr_models[6, ] runs quickest
#' model_fits <- run_mixmustr_models(
#'   mixmustr_models[6, ], synthetic_df_convergent, tracer_parameters,
#'   sigma_ln_rho = 0.1, iter = 1e4, warmup = 5e3, chains = 4, cores = 4
#' )
#' make_post_prop_long(model_fits[[1]]$model, tracer_parameters$mus,
#'                     synthetic_df_convergent, target = "df_stream_2", n = 1)
#' }
#'
#' @export
make_post_prop_long <- function(modfit, mu_tab, data_streams_list, n, ...) {
  array_of_draws <- extract(modfit)$p
  ordered_sources <- c(mu_tab$source, "Unsampled")
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
    reshape_ref_data(data_streams_list, order_ref = mu_tab$source, ...) |>
      data.frame(check.names = FALSE) |>
      mutate(N = seq_len(n())) |>
      pivot_longer(
        cols = !.data$N, names_to = "source", values_to = "Observed"
      ), by = c("N", "source")
  ) |>
    mutate(`Variant:` = as.character(.env$n))
}
