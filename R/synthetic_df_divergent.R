#' Synthetic Convergent Dataset
#'
#' This dataset is a synthetic dataset generated to simulate mixture data for testing and validation purposes in the `mixmustr` package. It is anchored to empirical values of stable isotopes and fatty acids for a range of plant carbon sources in marine soils stable isotope data (`bcs_si`), fatty acid data (`bcs_fa`), and synthetic mixture proportions. `synthetic_df_divergent` exhibits great difference in the underlying mixing proportions between data streams 1 and 2.
#'
#' @format A list containing three data frames:
#' \describe{
#'   \item{df_stream_1}{A data frame containing the simulated mixture data for the first stream, including tracer estimates calculated from the stable isotope and fatty acid data.}
#'   \item{df_stream_2}{A data frame containing the simulated synthetic
#' proportions for the second stream, a column for each source.}
#'   \item{stream_1_props}{A data frame containing the synthetic proportions for
#' the first stream, used as input for generating the mixture data. This is only
#' used for testing purposes}
#' }
#' @details
#' The dataset was generated using the non-exported `make_mixture_data` function, which combines stable isotope data (`bcs_si`), fatty acid data (`bcs_fa`), and synthetic proportions (`stream_1_props` and `stream_2_props`) to produce simulated mixture data. The `truth_stream` parameter determines which stream's proportions are used as the "true" source contributions.
#'
#' @source
#' The synthetic dataset was generated programmatically using the `mixmustr` package. The input data sources are:
#'
#' - Stable isotope data: See `bcs_si` documentation.
#' - Fatty acid data: See `bcs_fa` documentation.
#'
#' @seealso
#'   \code{\link{synthetic_df_convergent}}
#' 
#' @examples
#' data(synthetic_df_divergent)
#' str(synthetic_df_divergent)
"synthetic_df_divergent"
