#' Tracer Parameters
#'
#' This built-in list contains standardised mean (\eqn{\mu}), standard deviation
#' (\eqn{\sigma}) and sample size values for stable isotope and fatty acid
#' tracers, derived from the `bcs_si` and `bcs_fa` datasets. It is used
#' internally by the `mixmustr` package for modelling and analysis.
#'
#' @format A list with three elements:
#' \describe{
#'   \item{mus}{A data frame containing the mean (\eqn{\mu}) values for each
#' tracer, grouped by source. Columns represent tracers, and rows represent
#' sources.}
#'   \item{sigmas}{A data frame containing the standard deviation
#' (\eqn{\sigma}) values for each tracer, grouped by source. Columns represent
#' tracers, and rows represent sources.}
#'   \item{ns}{A data frame containing the sample size values for each tracer,
#' grouped by source. Columns represent tracers, and rows represent sources.}
#' }
#' @details
#' The `tracer_parameters` list is generated using the
#' `mixmustr:::wrangle_tracer_pars` function, which processes the `bcs_si` and
#' `bcs_fa` datasets to reshape and combine their tracer data. The `mus` element
#' contains the mean values for each tracer, the `sigmas` element contains
#' the corresponding standard deviations and the `ns` element contains the
#' sample sizes. These values are used as inputs for mixture modelling in the
#' `mixmustr` package.
#'
#' @source
#' The data is derived from the `bcs_si` and `bcs_fa` datasets. See their
#' respective documentation for more details.
#'
#' @examples
#' data(tracer_parameters)
#' str(tracer_parameters)
#'
#' # Access mean values
#' tracer_parameters$mus
#'
#' # Access standard deviation values
#' tracer_parameters$sigmas
#' 
#' # Access sample sizes
#' tracer_parameters$ns
"tracer_parameters"
