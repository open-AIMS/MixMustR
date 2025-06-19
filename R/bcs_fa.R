#' Fatty Acid Composition of Various Sources
#'
#' This dataset contains the fatty acid composition of various sources,
#' including macroalgae, mangroves, saltmarsh, seagrass, SPOM (suspended
#' particulate organic matter), and
#' terrestrial grass. The data includes mean values, standard deviations and
#' sample sizes for different fatty acids, as well as the study from which the
#' data was sourced. Values of 0 standard deviation have been arbitrarily
#' transformed to 0.01.
#'
#' @format A data frame with 6 rows and 15 variables:
#' \describe{
#'   \item{source}{The source of the sample (e.g., Macroalgae, Mangrove, etc.).}
#'   \item{Taxa}{The taxa or species associated with the source.}
#'   \item{24:0}{Mean value of the fatty acid 24:0.}
#'   \item{24:0 (SD)}{Standard deviation of the fatty acid 24:0.}
#'   \item{24:0 (n)}{Sample size of the fatty acid 24:0.}
#'   \item{18:1w9}{Mean value of the fatty acid 18:1w9.}
#'   \item{18:1w9 (SD)}{Standard deviation of the fatty acid 18:1w9.}
#'   \item{18:1w9 (n)}{Sample size of the fatty acid 18:1w9.}
#'   \item{18:2w6}{Mean value of the fatty acid 18:2w6.}
#'   \item{18:2w6 (SD)}{Standard deviation of the fatty acid 18:2w6.}
#'   \item{18:2w6 (n)}{Sample size of the fatty acid 18:2w6.}
#'   \item{18:3w3}{Mean value of the fatty acid 18:3w3.}
#'   \item{18:3w3 (SD)}{Standard deviation of the fatty acid 18:3w3.}
#'   \item{18:3w3 (n)}{Sample size of the fatty acid 18:3w3.}
#'   \item{20:4w6}{Mean value of the fatty acid 20:4w6.}
#'   \item{20:4w6 (SD)}{Standard deviation of the fatty acid 20:4w6.}
#'   \item{20:4w6 (n)}{Sample size of the fatty acid 20:4w6.}
#'   \item{20:5w3}{Mean value of the fatty acid 20:5w3.}
#'   \item{20:5w3 (SD)}{Standard deviation of the fatty acid 20:5w3.}
#'   \item{20:5w3 (n)}{Sample size of the fatty acid 20:5w3.}
#'   \item{Study}{The study from which the data was sourced.}
#' }
#' @source The data and methodology are described in the following references:
#' - Khotimchenko, S. V. (1991). Fatty acid composition of seven Sargassum species. *Phytochemistry*, 30(8), 2639-2641.
#' - Meziane, T., Lee, S. Y., Mfilinge, P. L., Shin, P. K. S., Lam, M. H. W., & Tsuchiya, M. (2007). Inter-specific and geographical variations in the fatty acid composition of mangrove leaves: implications for using fatty acids as a taxonomic tool and tracers of organic matter. *Marine Biology*, 150, 1103-1113. DOI:10.1007/s00227-006-0424-z
#' - Weete, J. D., Rivers, W. G., & Webert, D. J. (1970). Hydrocarbon and fatty acid distribution in the halophyte, *Salicornia bigelovii.* *Phytochemistry*, 9, 2041-2045.
#' - Belicka, L. L., Burkholder, D., Fourqurean, J. W., Heithaus, M. R., Macko, S. A., & Jaffé, R. (2012). Stable isotope and fatty acid biomarkers of seagrass, epiphytic, and algal organic matter to consumers in a pristine seagrass ecosystem. *Marine and Freshwater Research*, 63(3), 1085-1097. DOI:10.1071/MF12027
#' - Burns, K. A., Volkman, J. K., Cavanagh, J.-A., & Brinkman, D. (2003). Lipids as biomarkers for carbon cycling on the Northwest Shelf of Australia: results from a sediment trap study. *Marine Chemistry*, 80, 103-128.
#' - Wang, S., Chu, T., Huang, D., Li, B., & Wu, J. (2014). Incorporation of Exotic *Spartina alterniflora* into Diet of Deposit-Feeding Snails in the Yangtze River Estuary Salt Marsh: Stable Isotope and Fatty Acid Analyses. *Ecosystems*, 17, 567-577. DOI:10.1007/s10021-013-9743-3
#' 
#' @examples
#' data(bcs_fa)
#' head(bcs_fa)
"bcs_fa"
