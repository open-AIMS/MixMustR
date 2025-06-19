#' BCS Stable Isotope Data
#'
#' This dataset contains stable isotope ratios
#' (\eqn{\delta^{13}C} and \eqn{\delta^{15}N}), their standard deviations and
#' sample sizes for various sources in the BCS study. Notice that the sample
#' sizes for Terrestrial grass were not provided in the original source
#' (Goni & Thomas, 2000) and were therefore arbitrarily set to 15.
#'
#' @format A data frame with 6 rows and 6 variables:
#' \describe{
#'   \item{source}{The source of the sample (e.g., Mangrove, Seagrass, etc.).}
#'   \item{d(13C/12C)}{The stable isotope ratio of carbon (\eqn{\delta^{13}C})
#' in parts per thousand (‰).}
#'   \item{d(13C/12C) (SD)}{The standard deviation of the \eqn{\delta^{13}C}
#'  measurements.}
#'   \item{d(13C/12C) (n)}{The sample size of the \eqn{\delta^{13}C}
#' measurements.}
#'   \item{d(15N/14N)}{The stable isotope ratio of nitrogen
#' (\eqn{\delta^{15}N}) in parts per thousand (‰).}
#'   \item{d(15N/14N) (SD)}{The standard deviation of the \eqn{\delta^{15}N}
#' measurements.}
#'   \item{d(15N/14N) (n)}{The sample size of the \eqn{\delta^{15}N}
#' measurements.}
#'   \item{Study}{The study or reference from which the data were obtained.}
#' }
#' @source
#' The data and methodology are described in the following references:
#'
#' - Nyunja, J., Ntiba, M., & Onyari, J. (2009). Carbon and nitrogen isotopic variation in mangrove ecosystems along the Kenyan coast. *Estuarine, Coastal and Shelf Science*, 84(3), 377-385. DOI:10.1016/j.ecss.2009.06.001
#' - Bulmer, R. H., Lundquist, C. J., & Schwendenmann, L. (2020). Sediment carbon and nitrogen isotopes reveal variability in estuarine food web sources. *Marine Ecology Progress Series*, 641, 13-26. DOI:10.3354/meps13312
#' - Goni, M. A., & Thomas, K. A. (2000). Sources and transformations of organic matter in surface soils and sediments from a Spartina alterniflora marsh. *Estuarine, Coastal and Shelf Science*, 50(5), 657-679. DOI:10.1006/ecss.2000.0597
#'
#' @examples
#' data(bcs_si)
#' head(bcs_si)
"bcs_si"
