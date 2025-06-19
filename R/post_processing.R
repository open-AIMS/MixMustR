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
#' @param ... Further arguments passed to \code{\link{reshape_ref_data}}.
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

#' Compute a Bayesian version of R-squared for mixture models
#' 
#' @inheritParams mixmustr_wrangle_input
#' @inheritParams make_post_prop_long
#' @param summary Should summary statistics be returned
#' instead of the raw values? Default is `TRUE`.
#' 
#' @return If \code{summary = TRUE}, a (J + 1) x 4 matrix is returned
#'  (J + 1 = number of sampled + 1 unsampled source). Columns contain source
#'  names and summary statistics (mean +/- 95% HDI) of the Bayesian R-squared
#'  values. If \code{summary = FALSE}, the posterior draws of the Bayesian
#'  R-squared values are returned in an S x (J + 1) matrix (S is the number of
#' draws).
#'
#' @details For an introduction to the approach, see Gelman et al. (2019)
#'  and \url{https://github.com/jgabry/bayes_R2/}.
#'
#' @references Andrew Gelman, Ben Goodrich, Jonah Gabry & Aki Vehtari. (2019).
#'   R-squared for Bayesian regression models, \emph{The American Statistician},
#'   73(3):307-309. \code{10.1080/00031305.2018.1549100}
#' 
#' @importFrom rstan extract
#' @importFrom ggdist mean_hdci
#' @importFrom dplyr bind_rows select rename mutate across
#' @importFrom tidyselect where
#' @importFrom stats var
#' @importFrom rlang .data
#' 
#' @examples
#' \dontrun{
#' library(mixmustr)
#' # mixmustr_models[6, ] runs quickest
#' model_fits <- run_mixmustr_models(
#'   mixmustr_models[6, ], synthetic_df_convergent, tracer_parameters,
#'   sigma_ln_rho = 0.1, iter = 1e4, warmup = 5e3, chains = 4, cores = 4
#' )
#' # Now compute R2, using data stream 2 as point of comparison
#' mixmustr_bayes_R2(
#'   model_fits[[1]]$model, data_streams_list = synthetic_df_convergent,
#'   target = "df_stream_2", order_ref = tracer_parameters$mus$source
#' )
#' }
#'
#' @export
mixmustr_bayes_R2 <- function(modfit, summary = TRUE, ...) {
  ref_data <- reshape_ref_data(...)
  all_posts <- rstan::extract(modfit)$p
  new_dims <- c(ncol(all_posts), dim(all_posts)[3], nrow(all_posts))
  pred_y <- res_y <- array(0, dim = new_dims)
  for (i in seq_len(dim(all_posts)[3])) {
    for (j in seq_len(nrow(all_posts))) {
      pred_y[, i, j] <- all_posts[j, , i]
    }
  }
  for (i in seq_len(dim(res_y)[3])) {
    res_y[, , i] <- ref_data - pred_y[, , i]
  }
  out_r2 <- matrix(0, dim(res_y)[3], ncol(res_y))
  colnames(out_r2) <- colnames(ref_data)
  for (i in seq_len(nrow(out_r2))) {
    pred_i <- pred_y[, , i]
    res_i <- res_y[, , i]
    var_pred <- apply(pred_i, 2, var)
    var_res <- apply(res_i, 2, var)
    out_r2[i, ] <- var_pred / (var_pred + var_res)
  }
  if (summary) {
    apply(out_r2, 2, mean_hdci) |>
      bind_rows(.id = "Source") |>
      select(!(.data$`.width`:.data$`.interval`)) |>
      rename(mean = .data$y, `2.5%HDI` = .data$ymin, `97.5%HDI` = .data$ymax) |>
      mutate(across(where(is.numeric), ~round(.x, 2)))
  } else {
    out_r2
  }
}
