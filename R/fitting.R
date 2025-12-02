#' Builds underlying Stan code for MixMustR model of choice
#'
#' @param sample_tracer A logical vector, defaults to FALSE.
# Should the model estimate uncertainty around sampled sources signatures?
#' @param fix_unsampled A logical vector, defaults to FALSE.
# Should the model estimate uncertainty around the unsampled source
# signatures or should it be fixed to mean across all sampled sources?
#' @param hierarchical A logical vector, defaults to FALSE.
#' Should all observations be treated as independent or should the model include
#' a hierarchical grouping structure?
#' @param code_path A character vector indicating the file path to the
# corresponding Stan model file, based on the parameter combination.
#' 
#' @details `MixMustR` currently allows for eight model variants
#' which result from three user-driven binary choices: 1) should the model only
#' ingest the mean sampled-source tracer signatures (equivalent to
#' "residual-only error" structure of MixSIAR; Stock et al. 2018) or should it
#' incorporate their uncertainty based on user-provided mean, variance and
#' sample size information (equivalent to "process error" structure of MixSIAR);
#' 2) should the unsampled-source tracer signatures be fixed at the mean across
#' all sampled sources, or should they rather be estimated based on a prior
#' informed by the mean and variance across the sampled tracer signatures? and
#' 3) should all observations be treated as independent or should the model
#' include a hierarchical grouping structure?
#' 
#' Given these three choices, \code{\link{build_stancode}} saves a stan code
#' to an output directory of the users' choosing. Nothing is returned.
#' 
#' @references
#' Stock BC, Jackson AL, Ward EJ, Parnell AC, Phillips DL, Semmens BX (2018)
#' Analyzing mixing systems using a new generation of Bayesian tracer mixing
#' models. PeerJ, 6:e5096. doi:10.7717/peerj.5096.
#'
#' @return Nothing.
#' 
#' @seealso
#'   \code{\link{run_mixmustr_models}}
#' 
#' @examples
#' library(MixMustR)
#' data(mixmustr_models)
#'
#' build_stancode(mixmustr_models$sample_tracer[1],
#'                mixmustr_models$fix_unsampled[1],
#'                mixmustr_models$hierarchical[1],
#'                mixmustr_models$code_path[1])
#' 
#' @export
build_stancode <- function(sample_tracer = FALSE, fix_unsampled = FALSE,
                           hierarchical = FALSE, code_path) {
  if (any(sapply(c(sample_tracer, fix_unsampled, hierarchical),
                 function(x)!is.logical(x)))) {
    stop("Arguments 'sample_tracer', 'fix_unsampled' and 'hierarchical' must ",
         "be logical.")
  }
  if (missing(code_path)) {
    stop("Please provide a filename for the output stan code.")
  }
  script <- c(
    "// Code built by MixMustR on ", as.character(Sys.time()), "\n",
    "// User options:\n",
    paste0("\t// sample_tracer = ", sample_tracer, "\n"),
    paste0("\t// fix_unsampled = ", fix_unsampled, "\n"),
    paste0("\t// hierarchical = ", hierarchical, "\n"),
    "data {\n",
    "\tint<lower=1> N; // number of observations\n",
    "\tint<lower=1> J; // number of sampled (a.k.a considered) sources\n",
    "\tint<lower=1> K; // number of tracers\n",
    "\tmatrix[N, K] Y; // data stream 1\n",
    "\tmatrix[N, J + 1] ln_rho; // Log mixing proportions, data stream 2\n",
    "\tmatrix[N, J + 1] sigma_ln_rho; // Confidence around `ln_rho`\n",
    "\tmatrix[J, K] x; // sample means for each source and tracer\n"
  )
  if (sample_tracer) {
    script <- c(script,
      "\tmatrix<lower=0>[J, K] s; // sample SDs for each source and tracer\n",
      "\tmatrix<lower=1>[J, K] m; // source- and tracer- specific sample size\n"
    )
  }
  if (!fix_unsampled & !sample_tracer) {
    script <- c(script,
      "\tmatrix<lower=0>[J, K] s; // sample SDs for each source and tracer\n"
    )
  }
  if (hierarchical) {
    script <- c(script,
      "\tint<lower=1> R; // number of levels in the hierarchical structure\n",
      "\tarray[N] int<lower=1> YR; // dummy vector of random levels\n"
    )
  }
  script <- c(script, "}\n", "parameters {\n",
    "\tarray[N] vector[J + 1] zeta; // zeta from softmax\n"
  )
  if (hierarchical) {
    script <- c(script,
      "\tarray[R] vector[J + 1] randvec; // random effects\n",
      "\tarray[J + 1] real<lower=0> randsd; // random effects hyper parameter\n",
      "\tarray[J + 1] real<lower=0> ext_sigma; // residual error, likelihood 2\n"
    )
  }
  if (sample_tracer & fix_unsampled) {
    script <- c(script,
      "\tmatrix[J, K] mu; // mean for each source and tracer (J x K)\n",
      "\tmatrix<lower=0>[J, K] nu; // degrees of freedom\n"
    )
  }
  if (sample_tracer & !fix_unsampled) {
    script <- c(script,
      "\tmatrix[J + 1, K] mu; // mean for each source and tracer (J x K)\n",
      "\tmatrix<lower=0>[J, K] nu; // degrees of freedom\n"
    )
  }
  if (!sample_tracer) {
    script <- c(script,
      "\tarray[K] real<lower=0> sigma; // residual SDs for each tracer\n"
    )
  }
  if (!sample_tracer & !fix_unsampled) {
    script <- c(script,
      "\tarray[K] real sampled_mean_x; // mean unsampled estimate for each tracer\n"
    )
  }
  script <- c(script, "}\n", "transformed parameters {\n",
    "\tarray[N] simplex[J + 1] p; // probability vector for each site\n"
  )
  if (sample_tracer) {
    script <- c(script,
      "\tmatrix[J, K] omega; // SD for each source and tracer\n"
    )
  }
  script <- c(script,
    "\tfor (n in 1:N) {\n\t\tp[n] = softmax(zeta[n]);\n\t}\n"
  )
  if (sample_tracer) {
    script <- c(script,
      "\tfor (j in 1:J) {\n",
      "\t\tfor (k in 1:K) {\n",
      "\t\t\tomega[j, k] = sqrt((s[j, k]^2 * (m[j, k] - 1)) / nu[j, k]);\n",
      "\t\t}\n", "\t}\n"
    )
  }
  script <- c(script, "}\n", "model {\n")
  if (sample_tracer) {
    script <- c(script,
      "\t// Sampling sources\n",
      "\tfor (j in 1:J) {\n",
      "\t\tfor (k in 1:K) {\n",
      "\t\t\ttarget += chi_square_lpdf(m[j, k] | nu[j, k]);\n",
      "\t\t\ttarget += normal_lpdf(x[j, k] | mu[j, k], s[j, k] / sqrt(m[j, k]));\n",
      "\t\t}\n",
      "\t}\n"
    )
  }
  sd_term <- ifelse(sample_tracer, "sqrt(p[n][j]^2 * omega[j, k]^2)",
                    "sigma[k]")
  script <- c(script,
    "\t// Mixture model\n",
    "\tfor (n in 1:N) {\n",
    "\t\tfor (k in 1:K) {\n",
    "\t\t\tarray[J + 1] real phi;\n",
    "\t\t\t// Contribution from sampled sources\n",
    "\t\t\tfor (j in 1:J) {\n"
  )
  if (sample_tracer) {
    script <- c(script,
      "\t\t\t\tphi[j] = log(p[n][j]) + normal_lpdf(Y[n, k] | mu[j, k], "
    )
  } else {
    script <- c(script,
      "\t\t\t\tphi[j] = log(p[n][j]) + normal_lpdf(Y[n, k] | x[j, k], "
    )
  }
  script <- c(script, sd_term, ");\n",
    "\t\t\t}\n\t\t\t// Contribution from unsampled sources\n"
  )
  if (fix_unsampled & sample_tracer) {
    script <- c(script,
      "\t\t\treal mean_mu_k = mean(mu[, k]);\n",
      "\t\t\treal mean_omega_k = sqrt(sum(omega[, k]^2) / J);\n",
      "\t\t\tphi[J + 1] = log(p[n][J + 1]) + normal_lpdf(Y[n, k] | mean_mu_k, ",
      "sqrt(p[n][J + 1]^2 * mean_omega_k^2));\n"
    )
  } else if (fix_unsampled & !sample_tracer) {
    script <- c(script,
      "\t\t\treal mean_x_k = mean(x[, k]);\n",
      "\t\t\tphi[J + 1] = log(p[n][J + 1]) + normal_lpdf(Y[n, k] | mean_x_k, ",
      "sigma[k]);\n"
    )
  } else if (!fix_unsampled & sample_tracer) {
    script <- c(script,
      "\t\t\treal mean_mu_k = mean(mu[1:J, k]);\n",
      "\t\t\treal mean_omega_k = sqrt(sum(omega[, k]^2) / J);\n",
      "\t\t\tmu[J + 1, k] ~ normal(mean_mu_k, mean_omega_k);\n",
      "\t\t\tphi[J + 1] = log(p[n][J + 1]) + normal_lpdf(Y[n, k] | ",
      "mu[J + 1, k], sqrt(p[n][J + 1]^2 * mean_omega_k^2));\n"
    )
  } else if (!fix_unsampled & !sample_tracer) {
    script <- c(script,
      "\t\t\treal mean_x_k = mean(x[, k]);\n",
      "\t\t\treal mean_s_k = sqrt(sum(s[, k]^2)) / J;\n",
      "\t\t\tsampled_mean_x[k] ~ normal(mean_x_k, mean_s_k);\n",
      "\t\t\tphi[J + 1] = log(p[n][J + 1]) + normal_lpdf(Y[n, k] | ",
      "sampled_mean_x[k], sqrt(p[n][J + 1]^2 * mean_s_k^2));\n"
    )
  }
  script <- c(script,
    "\t\t\t// Log-sum-exp for numerical stability\n",
    "\t\t\ttarget += log_sum_exp(phi);\n\t\t}\n\t}\n"
  )
  script <- c(script,
    "\t// Independent evaluation of mixing proportion (zeta scale)\n",
    "\tfor (n in 1:N) {\n",
    "\t\tln_rho[n, ] ~ normal(zeta[n], sigma_ln_rho[n, ]);\n"
  )
  if (hierarchical) {
    script <- c(script,
      "\t\tzeta[n] ~ normal(randvec[YR[n]], ext_sigma);\n"
    )
  }
  script <- c(script, "\t}\n")
  if (hierarchical) {
    script <- c(script,
      "\tfor (j in 1:(J+1)) {\n",
      "\t\trandvec[][j] ~ normal(0, randsd[j]);\n",
      "\t}\n",
      "\trandsd ~ normal(0, 1);\n",
      "\text_sigma ~ normal(0, 1);\n"
    )
  }
  script <- c(script, "}")
  dir.create(dirname(code_path), recursive = TRUE, showWarnings = FALSE)
  message("Exporting Stan code to ", code_path)
  writeLines(paste0(script, collapse = ""), code_path)
}

#' @noRd
check_sigma_ln_rho <- function(x, ref) {
  er <- paste0("Parameter `sigma_ln_rho` needs to be either a single numeric",
               " value or a matrix of ", nrow(ref), " rows and ", ncol(ref),
                " columns, and cannot contain NA values.")
  if (is.matrix(x) && is.numeric(x)) {
    if (!all(dim(x) == dim(ref))) {
      stop(er)
    }
    if (sum(is.na(x)) > 0) {
      stop(er)
    }
  } else if (is.vector(x) && is.numeric(x)) {
    if (length(x) != 1) {
       stop(er)
    }
  } else {
    stop(er)
  }
}

#' @noRd
check_sd_tabs <- function(to_eval, mus, param = "SDs") {
  fct_eval <- if (param == "SDs") is.numeric else if (param == "ns") is.integer
  if (!is.data.frame(to_eval) | !is.data.frame(mus)) {
    stop("You need valid data.frames of tracers signature mean and ", param,
         ". See `tracer_parameters` for examples of each.")
  }
  if (!all(dim(to_eval) == dim(mus)) |
        !all(names(to_eval) == names(mus)) |
        !all(rownames(to_eval) == rownames(mus))
      ) {
    stop("Data frames of tracers signature mean and ", param, " do not match",
         " in structure.")
  }
  if (names(to_eval)[1] != "source") {
    stop("First column name in tracers signature mean and ", param, " should",
         " be `source`.")
  }
  if (!all(apply(to_eval[, -1], 2, fct_eval)) |
        !all(apply(to_eval[, -1], 2, function(x)sum(is.na(x)) == 0))) {
    type <- if (param == "SDs") "numeric" else if (param == "ns") "integer"
    stop("Tracers signature columns should all be ", type, " and cannot",
         " contain NAs.")
  }
}

#' Creates the input list for the Stan model
#'
#' @inheritParams build_stancode
#' @param data_streams_list A list containing the input
#' data streams for the model. It should include two
#' data frames named `df_stream_1` and`df_stream_2`. See
#' details for the expected structure of these.
#' @param tracer_list A named list containing 1--3 data frames of same size,
#' one for the mean signatures, one for their standard deviations, and another
#' one for sample size. The second data frame is mandatory if `sample_tracer`
#' is `TRUE`, or if `fix_unsampled` is `FALSE` and `sample_tracer` is `FALSE`.
#' The third data frame is mandatory if `sample_tracer` is `TRUE`. See details
#' for exact structure of these data frames.
#' @param model_path A character string specifying the path to the Stan model
#' file.
#' @param sigma_ln_rho A numeric value or matrix specifying the confidence
#' around the log mixing proportions. If a matrix, it must be of dimensions
#' N X J, with N being the number of observations and J being the number of
#' sources.
#' @param ... Additional unused arguments.
#'
#' @importFrom dplyr select mutate
#' @importFrom rlang .data
#' 
#' @details
#' The `mixmustr_wrangle_input` function prepares the input data for the Stan
#' model:
#'
#' - `data_streams_list`: A list containing two data frames:
#'   - `df_stream_1`: Contains the observed data. If this is a hierarchical
#' dataset, the first column should be named `group`, representing the group
#' labels, and the remaining columns should contain numeric tracer values.
#'   - `df_stream_2`: Contains the log mixing proportions. It must have the same
#' number of rows (`N`) as `df_stream_1`. If the design is hierarchical, the
#' first column should be named `group`, matching the `group` column in
#' `df_stream_1`. The remaining columns (e.g., `source1`, `source2`, etc.) must
#' contain proportions (i.e., row sums up to 1, and no negative values).
#'
#' - `tracer_list`: A list containing up to three data frames:
#' - `mus`: A data frame containing the mean tracer signatures for each
#' source. The first column should be named `source`, with source names
#' matching those from data_streams_list$df_stream_2. The remaining columns
#' should contain numeric tracer values, with names matching those in
#' data_streams_list$df_stream_1. This data frame is required across all models
#' in `MixMustR`.
#'
#' - `sigmas`: A data frame containing the standard deviations of tracer
#' signatures for each source. The structure must match that of `mus`. It
#' cannot contain NAs.
#' 
#' - `ns`: A data frame containing the sample size used to calculate mus and
#' sigmas The structure must match that of `mus`. It cannot contain NAs.
#'
#' - `model_path`: A character string specifying the path to the Stan model
#' file.
#'
#' - `sigma_ln_rho`: A numeric value or matrix specifying the confidence around
#' the log mixing proportions. If a matrix, it must have dimensions `N x (J +
#' 1)`, with N being the number of observations and J being the number of
#' sources. The additional J + 1 column represents the unsampled source, making
#' the final row sums across the J + 1 columns total 1.
#'
#' - `sample_tracer`: Logical. If `TRUE`, the model estimates uncertainty
#' around sampled source signatures.
#'
#' - `fix_unsampled`: Logical. If `TRUE`, the model fixes unsampled source signatures to the mean across all sampled sources.
#'
#' - `hierarchical`: Logical. If `TRUE`, the model includes a hierarchical grouping structure.
#' 
#' @return A list containing all the data elements necessary to run the Stan
#' code contained within `model_path`.
#' 
#' @seealso
#'   \code{\link{run_mixmustr_models}}
#' 
#' @examples
#' library(MixMustR)
#' 
#' # Example input data
#' synthetic_streams_list <- list(
#'   df_stream_1 = data.frame(
#'     group = c("A", "A", "B", "B"),
#'     tracer1 = c(1.2, 1.3, 2.1, 2.2),
#'     tracer2 = c(3.4, 3.5, 4.6, 4.7)
#'   ),
#'   df_stream_2 = data.frame(
#'     group = c("A", "A", "B", "B"),
#'     source1 = c(0.6, 0.4, 0.7, 0.3),
#'     source2 = c(0.4, 0.6, 0.3, 0.7)
#'   )
#' )
#' 
#' synthetic_tracer_list <- list(
#'   mus = data.frame(
#'     source = c("source1", "source2"),
#'     tracer1 = c(1.25, 2.15),
#'     tracer2 = c(3.45, 4.65)
#'   ),
#'   sigmas = data.frame(
#'     source = c("source1", "source2"),
#'     tracer1 = c(0.05, 0.1),
#'     tracer2 = c(0.15, 0.2)
#'   ),
#'   ns = data.frame(
#'     source = c("source1", "source2"),
#'     tracer1 = c(5L, 10L), # make sure these are integers
#'     tracer2 = c(7L, 9L)
#'   )
#' )
#' 
#' # Example parameters
#' model_path <- tempfile(fileext = ".stan")
#' sigma_ln_rho <- 0.01
#' sample_tracer <- TRUE
#' fix_unsampled <- FALSE
#' hierarchical <- TRUE
#' 
#' # Call the function
#' input_list <- mixmustr_wrangle_input(
#'   data_streams_list = synthetic_streams_list,
#'   tracer_list = synthetic_tracer_list,
#'   model_path = model_path,
#'   sigma_ln_rho = sigma_ln_rho,
#'   sample_tracer = sample_tracer,
#'   fix_unsampled = fix_unsampled,
#'   hierarchical = hierarchical
#' )
#' 
#' @export
mixmustr_wrangle_input <- function(data_streams_list, tracer_list,
                                   model_path, sigma_ln_rho, sample_tracer,
                                   fix_unsampled, hierarchical, ...) {
  mu_tab <- tracer_list$mus
  sig_tab <- tracer_list$sigmas
  n_tab <- tracer_list$ns
  yobs <- data_streams_list$df_stream_1 |>
    select(select(mu_tab, -.data$source) |> names())
  out <- list(
    N = nrow(yobs),
    J = nrow(mu_tab), # sources
    K = ncol(yobs), # tracers
    Y = yobs,
    x = select(mu_tab, -.data$source),
    ln_rho = reshape_ref_data(
      data_streams_list, target = "df_stream_2", order_ref = mu_tab$source
    ) |>
      as.matrix() |>
      abs_log(adjust = 0.0001)
  )
  check_sigma_ln_rho(sigma_ln_rho, yobs)
  if (is.vector(sigma_ln_rho)) {
    out[["sigma_ln_rho"]] <- matrix(
      sigma_ln_rho, nrow(out$ln_rho), ncol(out$ln_rho)
    )
  } else {
    out[["sigma_ln_rho"]] <- sigma_ln_rho
  }
  sig_tab_er <- "You need a valid data.frame of tracers signature SDs."
  ns_er <- "You need a valid data.frame of tracers signature sample sizes."
  if (sample_tracer) {
    if (is.null(sig_tab)) stop(sig_tab_er) else check_sd_tabs(sig_tab, mu_tab)
    out[["s"]] <- select(sig_tab, -.data$source)
    if (is.null(n_tab)) stop(ns_er) else check_sd_tabs(n_tab, mu_tab, "ns")
    out[["m"]] <- select(n_tab, -.data$source)
  }
  if (!fix_unsampled & !sample_tracer) {
    if (is.null(sig_tab)) stop(sig_tab_er) else check_sd_tabs(sig_tab, mu_tab)
    out[["s"]] <- select(sig_tab, -.data$source)
  }
  if (hierarchical) {
    if (!"group" %in% names(data_streams_list$df_stream_1)) {
      stop("You requested a hierarchical model. Your input data streams must",
           " contain a column named \"group\" indicating the grouping-level",
           " variable.")
    }
    reff_df <- data_streams_list$df_stream_1 |>
      mutate(YR = as.numeric(as.factor(.data$group)))
    out[["R"]] <- max(reff_df$YR)
    out[["YR"]] <- reff_df$YR
  }
  out
}

#' Runs a set of MixMustR models in Stan
#' 
#' @inheritParams mixmustr_wrangle_input
#' @param model_choices A data frame specifying the model configurations to
#' run. Each row should represent a model, with columns indicating the values
#' for `sample_tracer`, `fix_unsampled`, `hierarchical`, and `code_path`.
#' @param ... Additional arguments passed to \code{\link[rstan]{stan}}, such as 
#' `iter`, `chains`, or `control`.
#'
#' @details
#' The `run_mixmustr_models` function automates the process of running 
#' multiple `MixMustR` models based on user-specified configurations. It builds 
#' the Stan code for each model, prepares the input data, and runs the models 
#' sequentially. The results are returned as a list, with each element 
#' containing the run time and model output for a specific configuration. For
#' convenience, `MixMustR` provides the user with a built-in data frame
#' containing all of the potential allowed models \code{\link{mixmustr_models}}.
#'
#' @return A list where each element corresponds to a model configuration, 
#' containing the run time and model output.
#'
#' @seealso
#'   \code{\link{mixmustr_models}},
#'   \code{\link{mixmustr_wrangle_input}}
#'
#' @importFrom tools file_path_sans_ext
#' 
#' @examples
#' \dontrun{
#' library(MixMustR)
#' # mixmustr_models[6, ] runs quickest
#' model_fits <- run_mixmustr_models(
#'   mixmustr_models[6, ], synthetic_df_convergent, tracer_parameters,
#'   sigma_ln_rho = 0.1, iter = 1e4, warmup = 5e3, chains = 4, cores = 4
#' )
#' }
#' 
#' @export
run_mixmustr_models <- function(model_choices, data_streams_list, tracer_list,
                                sigma_ln_rho, ...) {
  mu_tab <- tracer_list$mus
  sig_tab <- tracer_list$sigmas
  n_tab <- tracer_list$ns
  sig_tab_required <- any(model_choices$sample_tracer | 
    (!model_choices$fix_unsampled & !model_choices$sample_tracer)
  )
  if (sig_tab_required) {
    sig_tab_er <- "You need a valid data.frame of tracers signature SDs."
    if (is.null(sig_tab)) stop(sig_tab_er) else check_sd_tabs(sig_tab, mu_tab)
  }
  models <- vector(mode = "list", length = nrow(model_choices))
  for (i in seq_len(nrow(model_choices))) {
    build_stancode(model_choices$sample_tracer[i],
                   model_choices$fix_unsampled[i],
                   model_choices$hierarchical[i],
                   model_choices$code_path[i])
    model_path <- model_choices$code_path[i]
    sdata <- mixmustr_wrangle_input(
      data_streams_list, tracer_list, model_path, sigma_ln_rho = sigma_ln_rho,
      sample_tracer = model_choices$sample_tracer[i],
      fix_unsampled = model_choices$fix_unsampled[i],
      hierarchical = model_choices$hierarchical[i]
    )
    mod_name_suffix_i <- file_path_sans_ext(basename(model_path))
    timing <- system.time({
      model <- run_mixmod(model_path, mod_name_suffix_i, data = sdata, ...)
    })
    models[[i]] <- list(timing = timing, model = model)
    names(models)[i] <- mod_name_suffix_i
  }
  models
}

#' Wrapper to run mixture model in Stan
#' 
#' @inheritParams mixmustr_wrangle_input
#' @param mod_name_suffix A character string used to create a unique model name.
#' @param ... Additional arguments passed to \code{\link[rstan]{stan}}.
#'
#' @importFrom tools file_path_sans_ext
#' @importFrom rstan stan
#' 
#' @seealso
#'   \code{\link{run_mixmustr_models}}
#' 
#' @export
run_mixmod <- function(model_path, mod_name_suffix, ...) {
  stan(file = model_path, model_name = mod_name_suffix, ...)
}
