#' mixmustr Models Configuration
#'
#' This data frame defines the configuration of models used in the `mixmustr` package. It specifies combinations of key parameters for generating and running Stan models, which are fundamental to the functions `mixmustr_wrangle_input` and `run_mixmustr_models`.
#'
#' @format A data frame with 8 rows and 4 variables:
#' \describe{
#'   \item{sample_tracer}{A logical vector, defaults to FALSE.
#' Should the model estimate uncertainty around sampled sources signatures?}
#'   \item{fix_unsampled}{A logical vector, defaults to FALSE.
#' Should the model estimate uncertainty around the unsampled source
#' signatures or should it be fixed to mean across all sampled sources?}
#'   \item{hierarchical}{A logical vector, defaults to FALSE.
#' Should all observations be treated as independent or should the model include
#' a hierarchical grouping structure?}
#'   \item{code_path}{A character vector indicating the file path to the
#' corresponding Stan model file, based on the parameter combination.}
#' }
#' @details
#' The `mixmustr_models` data frame is programmatically generated using all possible combinations of the `sample_tracer`, `fix_unsampled`, and `hierarchical` parameters. Each combination corresponds to a specific Stan model file, whose path is stored in the `code_path` column.
#'
#' This data frame is used internally by the `mixmustr` package to determine which Stan model to use for a given analysis, based on user-specified options.
#'
#' @source
#' The data frame is generated in the `data-raw/mixmustr_models.R` script.
#'
#' @examples
#' data(mixmustr_models)
#' head(mixmustr_models)
#' 
#' # Access the Stan model path for a specific configuration
#' mixmustr_models$code_path[mixmustr_models$sample_tracer & mixmustr_models$hierarchical]
"mixmustr_models"
