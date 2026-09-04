#' Synthetic Divergent Dataset
#'
#' This dataset is a synthetic dataset generated to simulate mixture data for testing and validation purposes in the `MixMustR` package. It is anchored to empirical values of stable isotopes and fatty acids for a range of plant carbon sources in marine soils stable isotope data (`bcs_si`), fatty acid data (`bcs_fa`), and synthetic mixture proportions. `synthetic_df_divergent` exhibits great difference in the underlying mixing proportions between data streams 1 and 2.
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
#' The non-exported `augment_bcs_with_unsampled` function first combines stable isotope data (`bcs_si`) and fatty acid data (`bcs_fa`), homogenises their sampled-source SD/sample-size columns, and appends a synthetic pooled-unsampled tracer signature (built by the non-exported `make_unsampled_signature` function) positioned away from the sampled-source multivariate space. The non-exported `simulate_mvn_mixture` function then generates the two composition streams (`stream_1_props` and `df_stream_2`) using a hierarchical logistic-normal process across 10 groups of 10 observations. The non-exported `make_mixture_data` function combines these signatures and proportions to produce the simulated tracer measurements in `df_stream_1`; the `truth_stream` parameter determines which stream's proportions are used as the "true" source contributions. Finally, the non-exported `rm_unsampled` function removes the explicit `Unsampled` proportion column from `df_stream_2` and `stream_1_props`, mimicking the realistic situation in which users do not supply it directly; \code{\link{mixmustr_wrangle_input}} reconstructs it internally.
#'
#' @source
#' The synthetic dataset was generated programmatically using the `MixMustR` package, following the simulation design described in the accompanying manuscript. The input data sources are:
#'
#' - Stable isotope data: See `bcs_si` documentation.
#' - Fatty acid data: See `bcs_fa` documentation.
#' - Pooled-unsampled tracer signature: synthetically constructed, offset from the sampled-source centroid (see `make_unsampled_signature`).
#'
#' @seealso
#'   \code{\link{synthetic_df_convergent}}
#' 
#' @examples
#' library(MixMustR)
#' data(synthetic_df_divergent)
#' str(synthetic_df_divergent)
"synthetic_df_divergent"
