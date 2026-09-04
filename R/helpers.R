#' Reshape a data stream into a full source-proportion matrix
#'
#' Extracts the sampled-source proportions for `target` from
#' `data_streams_list`, in the column order given by `order_ref`, and appends
#' the reconstructed unsampled-source proportion
#' (`1 - rowSums(sampled proportions)`).
#'
#' @inheritParams mixmustr_wrangle_input
#' @param target The target element of `data_streams_list`, typically
#' `"df_stream_2"`.
#' @param order_ref A character vector specifying the correct order of the
#' source names.
#' 
#' @importFrom dplyr mutate rowwise c_across
#' @importFrom tidyselect everything
#' 
#' @return A matrix of N x (J + 1), with the + 1 representing the
#' mixing proportions from the unsampled source(s).
#' 
#' @examples
#' library(MixMustR)
#' reshape_ref_data(synthetic_df_convergent, target = "df_stream_2",
#'                  order_ref = tracer_parameters$mus$source)
#' 
#' @export
reshape_ref_data <- function(data_streams_list, target = "df_stream_2",
                             order_ref) {
  (data_streams_list[[target]][, order_ref]) |>
    rowwise() |>
    mutate(Unsampled = 1 - sum(c_across(everything()))) |>
    as.matrix() |>
    abs()
}

#' @importFrom truncnorm rtruncnorm
#' @importFrom tidyr replace_na
#' @noRd
trun_na_zr <- function(...) {
  rtruncnorm(1, ...) |>
    replace_na(0)
}

#' @noRd
assign_new_names <- function(df_, new_names) {
  names(df_) <- new_names
  df_
}
