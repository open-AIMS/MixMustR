#' Build Stan code for a MixMustR model variant
#'
#' `build_stancode()` generates the Stan program for one MixMustR model
#' variant. The variant is defined by three logical choices: whether to estimate
#' sampled-source tracer signatures, whether to fix or estimate the
#' unsampled-source tracer signature, and whether to impose a hierarchical
#' structure on mixing proportions.
#'
#' @description
#' In this generated Stan program, mixing proportions are represented with a
# baseline-logit simplex in
#' which the unsampled source is the reference category. The first data stream
#' is modelled with a heteroscedastic tracer likelihood that propagates
#' sampled-source process variance and estimates an additional residual tracer
#' scale. The second data stream is modelled on the log-composition scale using
#' `log_softmax()` of the latent logits. Tracer observations and source
#' signatures are expected to be standardised by
#' \code{\link{mixmustr_wrangle_input}} prior to Stan fitting.
#'
#' @param sample_tracer Logical. Should sampled-source tracer signatures be
#'   estimated as latent source means? If `TRUE`, the generated Stan code
#'   includes non-centred sampled-source signature parameters with prior scale
#'   `s / sqrt(m)`. If `FALSE`, sampled-source signatures are fixed at the
#'   user-supplied means.
#' @param fix_unsampled Logical. Should the unsampled-source tracer signature be
#'   fixed? If `TRUE`, the generated Stan code fixes the unsampled signature to
#'   the mean sampled-source signature on the standardised tracer scale. If
#'   `FALSE`, the unsampled signature is estimated with a covariance-informed
#'   prior derived from sampled-source signatures.
#' @param hierarchical Logical. Should the model include a group-level
#'   hierarchical structure on mixing proportions? If `TRUE`, the generated Stan
#'   program expects `R` and `YR` in the Stan data list. If `FALSE`, observations
#'   are treated as independent while retaining observation-level variation in
#'   logits.
#' @param code_path Optional character scalar. If supplied and `write_file` is
#'   `TRUE`, the generated Stan program is written to this path.
#' @param write_file Logical. Should the generated Stan program be written to
#'   `code_path`? Defaults to `TRUE` when `code_path` is supplied and `FALSE`
#'   otherwise.
#'
#' @return A character scalar containing the generated Stan program. If
#'   `write_file = TRUE`, the same program is also written to `code_path` as a
#'   side effect.
#'
#' @details
#' The 2 x 2 x 2 model grid is defined by `sample_tracer`, `fix_unsampled`, and
#' `hierarchical`.
#'
#' For all variants, the generated model uses `J` free logits for `J` sampled
#' sources and fixes the unsampled-source logit to zero. This removes the
#' additive non-identifiability of a full `J + 1` softmax representation while
#' retaining a full `J + 1` simplex after applying `softmax()`.
#'
#' The stream-2 likelihood is applied to log probabilities, not to raw logits.
#' Internally, the model computes `log_softmax(eta)` and compares it with
#' `ln_rho`, the log-composition matrix produced by
#' \code{\link{mixmustr_wrangle_input}}.
#' The stream-2 likelihood uses a Student-t distribution with estimated degrees
#' of freedom `nu_rho` to provide robustness to stream-2 deviations.
#'
#' The stream-1 likelihood is heteroscedastic. For observation `n` and tracer
#' `k`, the model combines an estimated residual scale, `sigma_y[k]`, with the
#' sampled-source process-variance contribution implied by source proportions
#' and sampled-source standard deviations. The unsampled-source process variance
#' is not supplied by users and is therefore absorbed into `sigma_y[k]` together
#' with residual observation and model error.
#'
#' When `sample_tracer = TRUE`, sampled-source mean signatures are estimated
#' with a non-centred parametrisation:
#' `mu_source[j, k] = x[j, k] + (s[j, k] / sqrt(m[j, k])) * z_source[j, k]`.
#' When `fix_unsampled = FALSE`, the unsampled-source mean signature is 
#' estimated as `x_bar_unsampled_prior + tau_unsampled * L_x_unsampled_prior *
#' z_unsampled`, where the prior mean and Cholesky factor are constructed from
#' sampled-source signatures by \code{\link{mixmustr_wrangle_input}}.
#'
#' @references
#' Stock BC, Jackson AL, Ward EJ, Parnell AC, Phillips DL, Semmens BX (2018)
#' Analyzing mixing systems using a new generation of Bayesian tracer mixing
#' models. PeerJ, 6:e5096. doi:10.7717/peerj.5096.
#' 
#' @examples
#' # Build code for the model variant with fixed sampled-source signatures,
#' # estimated unsampled-source signature, and hierarchical mixing proportions.
#' library(MixMustR)
#' code <- build_stancode(
#'   sample_tracer = FALSE,
#'   fix_unsampled = FALSE,
#'   hierarchical = TRUE
#' )
#' cat(substr(code, 1, 80))
#'
#' # Build and write a non-hierarchical model with sampled-source signatures
#' # estimated and the unsampled-source signature fixed.
#' stan_file <- tempfile(fileext = ".stan")
#' code <- build_stancode(
#'   sample_tracer = TRUE,
#'   fix_unsampled = TRUE,
#'   hierarchical = FALSE,
#'   code_path = stan_file
#' )
#' file.exists(stan_file)
#'
#' # Generate the complete 2 x 2 x 2 grid of Stan programs without writing
#' # files.
#' model_grid <- expand.grid(
#'   sample_tracer = c(FALSE, TRUE),
#'   fix_unsampled = c(FALSE, TRUE),
#'   hierarchical = c(FALSE, TRUE)
#' )
#' stan_codes <- Map(
#'   build_stancode,
#'   sample_tracer = model_grid$sample_tracer,
#'   fix_unsampled = model_grid$fix_unsampled,
#'   hierarchical = model_grid$hierarchical
#' )
#'
#' @seealso
#'   \code{\link{mixmustr_wrangle_input}},
#'   \code{\link{run_mixmod}},
#'   \code{\link{run_mixmustr_models}},
#'   \code{\link{mixmustr_models}}
#' @export
build_stancode <- function(sample_tracer, fix_unsampled, hierarchical,
                           code_path = NULL, write_file = !is.null(code_path)) {
  sample_tracer <- isTRUE(sample_tracer)
  fix_unsampled <- isTRUE(fix_unsampled)
  hierarchical <- isTRUE(hierarchical)
  out <- character()
  add <- function(...) {
    out <<- c(out, ...)
  }
  variant_label <- paste0(
    "sample_tracer_", toupper(as.character(sample_tracer)),
    "_fix_unsampled_", toupper(as.character(fix_unsampled)),
    "_hierarchical_", toupper(as.character(hierarchical))
  )
  add(
    paste0("// ", variant_label),
    paste0("// sample_tracer = ", toupper(as.character(sample_tracer)),
           ", fix_unsampled = ", toupper(as.character(fix_unsampled)),
           ", hierarchical = ", toupper(as.character(hierarchical))),
    paste0("// Generated by `MixMustR::build_stancode` on ",
           as.character(Sys.time()), "."),
    paste0("// Expected input scale: tracer observations and source signatures",
           " standardised in `MixMustR::mixmustr_wrangle_input`."),
    "",
    "data {",
    "  int<lower=1> N;",
    "  int<lower=1> J;                 // sampled sources",
    "  int<lower=1> K;                 // tracers"
  )
  if (hierarchical) {
    add(
      "  int<lower=1> R;                 // groups",
      "  array[N] int<lower=1, upper=R> YR;"
    )
  }
  add(
    "  array[N] vector[K] Y;           // stream-1 tracer observations, standardised",
    "",
    "  matrix[J, K] x;                 // sampled-source mean signatures, standardised",
    "  matrix<lower=0>[J, K] s;        // sampled-source SDs, standardised",
    "  array[J, K] int<lower=1> m;     // sampled-source sample sizes",
    "",
    "  matrix[N, J + 1] ln_rho;        // log stream-2 proportions, sampled + unsampled",
    "  matrix<lower=0>[N, J + 1] sigma_ln_rho;"
  )
  if (!fix_unsampled) {
    add(
      "",
      "  // Prior for estimated unsampled signature on the standardised tracer scale.",
      "  vector[K] x_bar_unsampled_prior;",
      "  matrix[K, K] L_x_unsampled_prior;"
    )
  }
  add("}", "")
  if (fix_unsampled) {
    add(
      "transformed data {",
      "  vector[K] x_bar;",
      "",
      "  // Fixed unsampled-source signature: mean sampled-source signature.",
      "  for (k in 1:K) {",
      "    real acc = 0;",
      "    for (j in 1:J) {",
      "      acc += x[j, k];",
      "    }",
      "    x_bar[k] = acc / J;",
      "  }",
      "}",
      ""
    )
  }
  add("parameters {")
  if (hierarchical) {
    add(
      "  // Baseline-logit hierarchy for mixing proportions.",
      "  // Only J logits are estimated; the unsampled source is the reference category.",
      "  matrix[R, J] alpha_raw;",
      "  matrix[N, J] z_raw;",
      "  vector<lower=0>[J] sigma_alpha;",
      "  vector<lower=0>[J] sigma_z;",
      "  vector[J] beta0;"
    )
  } else {
    add(
      "  // Baseline-logit model for independent observations.",
      "  // Only J logits are estimated; the unsampled source is the reference category.",
      "  matrix[N, J] z_raw;",
      "  vector<lower=0>[J] sigma_z;",
      "  vector[J] beta0;"
    )
  }
  add("")
  if (sample_tracer) {
    add(
      "  // Non-centred sampled-source mean signatures.",
      "  matrix[J, K] z_source;",
      ""
    )
  }
  if (!fix_unsampled) {
    add(
      "  // Estimated unsampled-source mean signature.",
      "  vector[K] z_unsampled;",
      "  real<lower=0> tau_unsampled;",
      ""
    )
  }
  add(
    "  // Residual stream-1 tracer SD on standardised scale.",
    "  vector[K] log_sigma_y_raw;",
    "",
    "  // Student-t degrees of freedom for stream-2 log-proportion likelihood.",
    "  real<lower=2> nu_rho;",
    "}",
    "",
    "transformed parameters {"
  )
  if (hierarchical) {
    add("  matrix[R, J] alpha;")
  }
  add(
    "  matrix[N, J] zeta;",
    "  matrix[N, J + 1] p;",
    ""
  )
  if (sample_tracer) {
    add("  matrix[J, K] mu_source;")
  }
  if (!fix_unsampled) {
    add("  vector[K] mu_unsampled;")
  }
  add(
    "  matrix[J + 1, K] all_mus;",
    "  vector<lower=0>[K] sigma_y;",
    "",
    "  // Positive-floor residual scale. This avoids zero-boundary geometry.",
    "  for (k in 1:K) {",
    "    sigma_y[k] = 0.01 + 0.05 * exp(log_sigma_y_raw[k]);",
    "  }",
    ""
  )
  if (hierarchical) {
    add(
      "  // Group-level logits.",
      "  for (g in 1:R) {",
      "    for (j in 1:J) {",
      "      alpha[g, j] = beta0[j] + sigma_alpha[j] * alpha_raw[g, j];",
      "    }",
      "  }",
      "",
      "  // Observation-level logits and simplex proportions.",
      "  for (n in 1:N) {",
      "    vector[J + 1] eta;",
      "    for (j in 1:J) {",
      "      zeta[n, j] = alpha[YR[n], j] + sigma_z[j] * z_raw[n, j];",
      "      eta[j] = zeta[n, j];",
      "    }",
      "    eta[J + 1] = 0;",
      "    p[n] = to_row_vector(softmax(eta));",
      "  }",
      ""
    )
  } else {
    add(
      "  // Observation-level logits and simplex proportions.",
      "  for (n in 1:N) {",
      "    vector[J + 1] eta;",
      "    for (j in 1:J) {",
      "      zeta[n, j] = beta0[j] + sigma_z[j] * z_raw[n, j];",
      "      eta[j] = zeta[n, j];",
      "    }",
      "    eta[J + 1] = 0;",
      "    p[n] = to_row_vector(softmax(eta));",
      "  }",
      ""
    )
  }
  if (sample_tracer) {
    add(
      "  // Sampled-source signatures: non-centred SE-scale prior.",
      "  for (j in 1:J) {",
      "    for (k in 1:K) {",
      "      mu_source[j, k] = x[j, k] + (s[j, k] / sqrt(m[j, k])) * z_source[j, k];",
      "      all_mus[j, k] = mu_source[j, k];",
      "    }",
      "  }"
    )
  } else {
    add(
      "  // Fixed sampled-source signatures.",
      "  for (j in 1:J) {",
      "    for (k in 1:K) {",
      "      all_mus[j, k] = x[j, k];",
      "    }",
      "  }"
    )
  }
  if (fix_unsampled) {
    add(
      "",
      "  // Fixed unsampled-source signature.",
      "  for (k in 1:K) {",
      "    all_mus[J + 1, k] = x_bar[k];",
      "  }"
    )
  } else {
    add(
      "",
      "  // Estimated unsampled-source signature with covariance-informed prior.",
      "  mu_unsampled = x_bar_unsampled_prior + tau_unsampled * L_x_unsampled_prior * z_unsampled;",
      "  for (k in 1:K) {",
      "    all_mus[J + 1, k] = mu_unsampled[k];",
      "  }"
    )
  }
  add("}", "", "model {")
  if (hierarchical) {
    add(
      "  // Hierarchical composition priors.",
      "  to_vector(alpha_raw) ~ std_normal();",
      "  to_vector(z_raw) ~ std_normal();",
      "  sigma_alpha ~ normal(0, 1);",
      "  sigma_z ~ normal(0, 1);",
      "  beta0 ~ normal(0, 1.5);"
    )
  } else {
    add(
      "  // Independent-observation composition priors.",
      "  to_vector(z_raw) ~ std_normal();",
      "  sigma_z ~ normal(0, 1);",
      "  beta0 ~ normal(0, 1.5);"
    )
  }
  add("")
  if (sample_tracer) {
    add("  to_vector(z_source) ~ std_normal();")
  }
  if (!fix_unsampled) {
    add(
      "  z_unsampled ~ std_normal();",
      "  tau_unsampled ~ normal(0, 1);"
    )
  }
  add(
    "  log_sigma_y_raw ~ normal(0, 0.5);",
    "  nu_rho ~ gamma(2, 0.1);",
    "",
    "  // Stream 2: log-composition likelihood on log probabilities.",
    "  for (n in 1:N) {",
    "    vector[J + 1] eta;",
    "    vector[J + 1] log_p_n;",
    "    for (j in 1:J) {",
    "      eta[j] = zeta[n, j];",
    "    }",
    "    eta[J + 1] = 0;",
    "    log_p_n = log_softmax(eta);",
    "    for (j in 1:(J + 1)) {",
    "      ln_rho[n, j] ~ student_t(nu_rho, log_p_n[j], sigma_ln_rho[n, j]);",
    "    }",
    "  }",
    "",
    "  // Stream 1: heteroscedastic tracer likelihood.",
    "  // The sampled-source process component is propagated via s[j, k].",
    "  // Unknown unsampled-source process variance is absorbed by sigma_y[k].",
    "  for (n in 1:N) {",
    "    vector[K] mu_n;",
    "    vector[K] sigma_n;",
    "    for (k in 1:K) {",
    "      real process_var = square(sigma_y[k]);",
    "      for (j in 1:J) {",
    "        process_var += square(p[n, j]) * square(s[j, k]);",
    "      }",
    "      sigma_n[k] = sqrt(process_var);",
    "      mu_n[k] = dot_product(to_vector(p[n]), col(all_mus, k));",
    "    }",
    "    Y[n] ~ normal(mu_n, sigma_n);",
    "  }",
    "}",
    "",
    "generated quantities {",
    "  matrix[N, K] mix_mean;",
    "  vector[N] log_lik;",
    "",
    "  for (n in 1:N) {",
    "    vector[K] mu_n;",
    "    vector[K] sigma_n;",
    "    for (k in 1:K) {",
    "      real process_var = square(sigma_y[k]);",
    "      for (j in 1:J) {",
    "        process_var += square(p[n, j]) * square(s[j, k]);",
    "      }",
    "      sigma_n[k] = sqrt(process_var);",
    "      mix_mean[n, k] = dot_product(to_vector(p[n]), col(all_mus, k));",
    "      mu_n[k] = mix_mean[n, k];",
    "    }",
    "    log_lik[n] = normal_lpdf(Y[n] | mu_n, sigma_n);",
    "  }",
    "}"
  )
  stan_code <- paste(out, collapse = "\n")
  if (isTRUE(write_file)) {
    if (is.null(code_path) || !nzchar(code_path)) {
      stop("`code_path` must be supplied when `write_file = TRUE`.",
           call. = FALSE)
    }
    dir.create(dirname(code_path), recursive = TRUE, showWarnings = FALSE)
    writeLines(stan_code, con = code_path, useBytes = TRUE)
  }
  stan_code
}

#' Standardise tracer data before Stan fitting
#'
#' `standardise_tracer_sdata()` standardises stream-1 tracer observations,
#' sampled-source signatures, and sampled-source standard deviations in a Stan
#' data list.
#'
#' @param sdata A named Stan data list containing at least `Y` and `x`. `Y` must
#'   be an `N x K` matrix-like object of stream-1 tracer observations and `x`
#'   must be a `J x K` matrix-like object of sampled-source mean signatures. If
#'   present, `s` is treated as a `J x K` matrix of sampled-source standard
#'   deviations.
#'
#' @return A modified `sdata` list. Elements `Y` and `x` are centred and scaled
#'   by tracer using only sampled-source signatures. If `s` is present, it is
#'   divided by the tracer scale but not centred. The returned list has
#'   `tracer_center` and `tracer_scale` attributes for downstream
#'   back-transformation.
#'
#' @details
#' Standardisation is based only on sampled-source signatures. This is a
#' deliberate production choice: users do not know the unsampled-source
#' signature, and in some applications may not know whether an unsampled source
#' exists. Each tracer is centred by the across-source mean in `sdata$x` and
#' scaled by the across-source standard deviation. Non-finite or zero scales are
#' replaced by one.
#'
#' Standard deviations are scale parameters, so `sdata$s` is divided by the
#' tracer scale without centring. The centring and scaling constants are stored
#' as attributes rather than Stan data variables.
#'
#' @importFrom stats sd
#' @keywords internal
#' @noRd
standardise_tracer_sdata <- function(sdata) {
  # Assumes:
  #   sdata$Y is N x K
  #   sdata$x is J x K
  #   sdata$s is J x K, if sampled-source process uncertainty is used
  #
  # Important production behaviour:
  #   Standardisation is based only on sampled-source signatures `sdata$x`.
  #   No unsampled-source signature is ever used here.
  if (!"Y" %in% names(sdata)) {
    stop("`sdata` must contain `Y` before standardisation.", call. = FALSE)
  }
  if (!"x" %in% names(sdata)) {
    stop("`sdata` must contain `x` before standardisation.", call. = FALSE)
  }
  Y <- as.matrix(sdata$Y)
  X <- as.matrix(sdata$x)
  tracer_center <- colMeans(X, na.rm = TRUE)
  tracer_scale <- apply(X, 2, sd, na.rm = TRUE)
  tracer_scale[!is.finite(tracer_scale) | tracer_scale == 0] <- 1
  sdata$Y <- sweep(sweep(Y, 2, tracer_center, "-"), 2, tracer_scale, "/")
  sdata$x <- sweep(sweep(X, 2, tracer_center, "-"), 2, tracer_scale, "/")
  # SDs are scale parameters: divide by scale only; do not centre.
  if ("s" %in% names(sdata)) {
    sdata$s <- sweep(as.matrix(sdata$s), 2, tracer_scale, "/")
  }
  # Keep back-transformation information for downstream summaries / plots.
  # These attributes are not Stan data variables.
  attr(sdata, "tracer_center") <- tracer_center
  attr(sdata, "tracer_scale") <- tracer_scale
  sdata
}

#' Build an unsampled-source covariance prior from sampled signatures
#'
#' `make_unsampled_prior_cov()` constructs the prior mean and covariance factor
#' used by Stan model variants that estimate the unsampled-source tracer
#' signature.
#'
#' @param sdata A standardised Stan data list containing sampled-source
#'   signatures in `sdata$x`.
#' @param ridge Numeric scalar. Ridge term added to the diagonal of the empirical
#'   sampled-source covariance matrix before Cholesky decomposition. Defaults to
#'   `0.1`.
#'
#' @return A list with two elements: `x_bar`, the mean sampled-source signature
#'   on the standardised tracer scale, and `L_x`, the lower Cholesky factor of
#'   the ridge-regularised sampled-source covariance matrix.
#'
#' @details
#' This helper is used only for model variants with `fix_unsampled = FALSE`. The
#' resulting objects are passed to Stan as `x_bar_unsampled_prior` and
#' `L_x_unsampled_prior`. The ridge term stabilises the covariance matrix when
#' the number of sampled sources is small relative to the number of tracers.
#'
#' @importFrom stats cov
#' @keywords internal
#' @noRd
make_unsampled_prior_cov <- function(sdata, ridge = 0.1) {
  if (!"x" %in% names(sdata)) {
    stop("`sdata` must contain standardised `x` before constructing the unsampled prior.", call. = FALSE)
  }
  X <- as.matrix(sdata$x)
  # Covariance among sampled-source signatures on the standardised tracer scale.
  V <- cov(X)
  # Add a ridge for stability because J can be small relative to K.
  V_reg <- V + diag(ridge^2, ncol(V))
  L <- t(chol(V_reg))
  list(x_bar = colMeans(X), L_x = L)
}

#' Standardise Stan data and append unsampled-prior data
#'
#' `standardisation_wrapper()` applies tracer standardisation to a Stan data
#' list and, when required, appends the covariance-informed prior objects for an
#' estimated unsampled-source signature.
#'
#' @param sdata A Stan data list assembled by
#'   \code{\link{mixmustr_wrangle_input}} before
#'   standardisation.
#' @param fix_unsampled Logical. If `FALSE`, `x_bar_unsampled_prior` and
#'   `L_x_unsampled_prior` are added to the returned list. If `TRUE`, those
#'   elements are omitted because fixed-unsampled Stan variants do not declare
#'   them.
#' @param unsampled_prior_ridge Numeric scalar. Ridge term passed to
#'   `make_unsampled_prior_cov()`. Defaults to `0.1`.
#'
#' @return A standardised Stan data list. The returned list also carries
#'   `tracer_center` and `tracer_scale` attributes for back-transformation.
#'
#' @details
#' This is an internal helper that centralises the preprocessing required by the
#' updated Stan code generated by \code{\link{build_stancode}}. It first
#' standardises `Y`, `x`, and, where present, `s`. It then adds unsampled-prior
#' data only for variants that estimate the unsampled-source signature.
#'
#' @keywords internal
#' @noRd
standardisation_wrapper <- function(sdata, fix_unsampled,
                                    unsampled_prior_ridge = 0.1) {
  sdata <- standardise_tracer_sdata(sdata = sdata)
  # Estimated-unsampled models need the covariance-informed prior data.
  # Fixed-unsampled models do not declare these Stan data variables.
  if (!isTRUE(fix_unsampled)) {
    prior_u <- make_unsampled_prior_cov(
      sdata = sdata, ridge = unsampled_prior_ridge
    )
    sdata$x_bar_unsampled_prior <- prior_u$x_bar
    sdata$L_x_unsampled_prior <- prior_u$L_x
  }
  sdata
}

#' @noRd
check_sd_tabs <- function(to_eval, mus, param = "SDs") {
  fct_eval <- if (param == "SDs") is.numeric else if (param == "ns") is.integer
  if (!is.data.frame(to_eval) || !is.data.frame(mus)) {
    stop("You need valid data.frames of tracers signature mean and ", param,
         ". See `tracer_parameters` for examples of each.")
  }
  if (!all(dim(to_eval) == dim(mus)) ||
        !all(names(to_eval) == names(mus)) ||
        !all(rownames(to_eval) == rownames(mus))
      ) {
    stop("Data frames of tracers signature mean and ", param, " do not match",
         " in structure.")
  }
  if (names(to_eval)[1] != "source") {
    stop("First column name in tracers signature mean and ", param, " should",
         " be `source`.")
  }
  if (!all(apply(to_eval[, -1], 2, fct_eval)) ||
        !all(apply(to_eval[, -1], 2, function(x) sum(is.na(x)) == 0))) {
    type <- if (param == "SDs") "numeric" else if (param == "ns") "integer"
    stop("Tracers signature columns should all be ", type, " and cannot",
         " contain NAs.")
  }
}

#' Create the Stan data list for a MixMustR model
#'
#' `mixmustr_wrangle_input()` validates MixMustR input data and converts the two
#' data streams and sampled-source tracer summaries into the Stan data list used
#' by the updated MixMustR model variants.
#'
#' @param data_streams_list A named list containing `df_stream_1` and
#'   `df_stream_2`. `df_stream_1` must contain stream-1 tracer observations with
#'   one column per tracer and may contain a `group` column. `df_stream_2` must
#'   contain one column per sampled source with values interpreted as stream-2
#'   mixing proportions and may contain a `group` column. An `Unsampled` column
#'   is optional; if absent, the unsampled component is reconstructed as
#'   `1 - rowSums(sampled proportions)`.
#' @param tracer_list A named list containing `mus`, `sigmas`, and `ns`. Each
#'   element must contain a `source` column followed by matching tracer columns.
#'   `mus` supplies sampled-source mean signatures, `sigmas` supplies
#'   sampled-source standard deviations, and `ns` supplies sampled-source sample
#'   sizes.
#' @param model_path Optional character scalar retained for compatibility with
#'   higher-level wrappers. The current wrangling logic does not inspect this
#'   value.
#' @param sigma_ln_rho Numeric scalar, numeric vector, or numeric matrix giving
#'   uncertainty in stream-2 log proportions. A scalar is recycled to an
#'   `N x (J + 1)` matrix. A vector must have length `J + 1` and is recycled
#'   across observations. A matrix must have dimensions `N x (J + 1)`.
#' @param sample_tracer Logical. Included for compatibility with model-choice
#'   wrappers. The returned data list always includes `x`, `s`, and `m` because
#'   the updated likelihood uses sampled-source process variance and sampled
#'   signature model variants use standard errors.
#' @param fix_unsampled Logical. Should the unsampled-source signature be fixed
#'   in the Stan model? If `FALSE`, covariance-informed unsampled-prior objects
#'   are appended to the Stan data list. If `TRUE`, those objects are omitted.
#' @param hierarchical Logical. Should hierarchical group data be included? If
#'   `TRUE`, the returned list includes `R` and `YR`. If no `group` column is
#'   available, all observations are assigned to a single group.
#' @param unsampled_prior_ridge Numeric scalar. Ridge term used when 
#'   constructing the covariance-informed unsampled-source prior. Defaults to
#'   `0.1`.
#' @param eps Numeric scalar. Small positive floor applied to stream-2
#'   proportions before renormalisation and log transformation. Defaults to
#'   `1e-12`.
#'
#' @return A named list suitable for passing to Stan. The list contains `N`,
#'   `J`, `K`, `Y`, `x`, `s`, `m`, `ln_rho`, and `sigma_ln_rho`. If
#'   `hierarchical = TRUE`, it also contains `R` and `YR`. If
#'   `fix_unsampled = FALSE`, it additionally contains `x_bar_unsampled_prior`
#'   and `L_x_unsampled_prior`. The returned list carries `tracer_center` and
#'   `tracer_scale` attributes for downstream back-transformation.
#'
#' @details
#' The function enforces consistent tracer and source ordering before building
#' the Stan data list. Tracer columns are taken from `tracer_list$mus`, excluding
#' the `source` column, and must be present in `df_stream_1`. Sampled-source
#' columns are taken from `tracer_list$mus$source` and must be present in
#' `df_stream_2`.
#'
#' Stream-2 proportions are converted to an `N x (J + 1)` composition by
#' appending or reconstructing the unsampled-source component. If `df_stream_2`
#' contains an `Unsampled` column, that column is used. Otherwise the unsampled
#' component is reconstructed as `1 - rowSums(P_sampled)`. Proportions are
#' floored at `eps`, renormalised to sum to one, and then transformed with
#' `log()`. The resulting `ln_rho` is therefore compatible with the
#' `log_softmax()` likelihood generated by \code{\link{build_stancode}}.
#'
#' Tracer quantities are standardised after the raw Stan data list is assembled.
#' Standardisation uses only sampled-source signatures. This avoids using any
#' information about the unsampled-source signature that would not be available
#' in real applications. For estimated-unsampled models, the unsampled-source
#' prior mean and covariance factor are derived from the standardised
#' sampled-source signature matrix.
#'
#' When `hierarchical = TRUE`, group indices are taken from `df_stream_1$group`
#' when available, otherwise from `df_stream_2$group` when available. If neither
#' stream contains a `group` column, all observations are assigned to a single
#' group. When `hierarchical = FALSE`, no group data are added to the Stan list.
#'
#' @examples
#' \dontrun{
#' # Hierarchical model with fixed sampled-source means and an estimated
#' # unsampled-source signature.
#' library(MixMustR)
#' sdata <- mixmustr_wrangle_input(
#'   data_streams_list = synthetic_df_convergent,
#'   tracer_list = tracer_parameters,
#'   sigma_ln_rho = 1,
#'   sample_tracer = FALSE,
#'   fix_unsampled = FALSE,
#'   hierarchical = TRUE
#' )
#'
#' names(sdata)
#' attr(sdata, "tracer_center")
#' attr(sdata, "tracer_scale")
#'
#' # Non-hierarchical model with sampled-source signatures estimated and the
#' # unsampled-source signature fixed to the sampled-source centroid.
#' sdata_no_hier <- mixmustr_wrangle_input(
#'   data_streams_list = synthetic_df_convergent,
#'   tracer_list = tracer_parameters,
#'   sigma_ln_rho = 0.5,
#'   sample_tracer = TRUE,
#'   fix_unsampled = TRUE,
#'   hierarchical = FALSE
#' )
#'
#' # Source-specific stream-2 uncertainty can be supplied as a vector of length
#' # J + 1, including the unsampled component.
#' sigma_vec <- rep(0.5, length(tracer_parameters$mus$source) + 1)
#' sdata_vec <- mixmustr_wrangle_input(
#'   data_streams_list = synthetic_df_convergent,
#'   tracer_list = tracer_parameters,
#'   sigma_ln_rho = sigma_vec,
#'   sample_tracer = FALSE,
#'   fix_unsampled = FALSE,
#'   hierarchical = TRUE
#' )
#' }
#'
#' @seealso
#'   \code{\link{build_stancode}},
#'   \code{\link{run_mixmod}},
#'   \code{\link{run_mixmustr_models}},
#'   \code{\link{mixmustr_models}}
#' @export
mixmustr_wrangle_input <- function(data_streams_list, tracer_list,
                                   model_path = NULL, sigma_ln_rho,
                                   sample_tracer, fix_unsampled, hierarchical,
                                   unsampled_prior_ridge = 0.1, eps = 1e-12) {
  # Expected structures:
  #   data_streams_list$df_stream_1: N rows, optional `group`, K tracer columns
  #   data_streams_list$df_stream_2: N rows, optional `group`, J sampled-source
  #     proportion columns; optional `Unsampled` column is allowed but not
  #     required.
  #   tracer_list$mus:    source column + K tracer means
  #   tracer_list$sigmas: source column + K tracer SDs
  #   tracer_list$ns:     source column + K tracer sample sizes
  if (!"df_stream_1" %in% names(data_streams_list)) {
    stop("`data_streams_list` must contain `df_stream_1`.", call. = FALSE)
  }
  if (!"df_stream_2" %in% names(data_streams_list)) {
    stop("`data_streams_list` must contain `df_stream_2`.", call. = FALSE)
  }
  if (!"mus" %in% names(tracer_list)) {
    stop("`tracer_list` must contain `mus`.", call. = FALSE)
  }
  if (!"sigmas" %in% names(tracer_list)) {
    stop("`tracer_list` must contain `sigmas` for the model likelihood.",
         call. = FALSE)
  }
  if (!"ns" %in% names(tracer_list)) {
    stop("`tracer_list` must contain `ns`.", call. = FALSE)
  }
  df_stream_1 <- data_streams_list$df_stream_1
  df_stream_2 <- data_streams_list$df_stream_2
  source_order <- as.character(tracer_list$mus$source)
  tracer_names <- setdiff(names(tracer_list$mus), "source")
  if (!all(tracer_names %in% names(df_stream_1))) {
    missing_tracers <- setdiff(tracer_names, names(df_stream_1))
    stop(
      "`df_stream_1` is missing tracer column(s): ",
      paste(missing_tracers, collapse = ", "), call. = FALSE
    )
  }
  if (!all(source_order %in% names(df_stream_2))) {
    missing_sources <- setdiff(source_order, names(df_stream_2))
    stop(
      "`df_stream_2` is missing sampled-source column(s): ",
      paste(missing_sources, collapse = ", "),
      call. = FALSE
    )
  }
  Y <- as.matrix(df_stream_1[, tracer_names, drop = FALSE])
  x <- as.matrix(tracer_list$mus[, tracer_names, drop = FALSE])
  s <- as.matrix(tracer_list$sigmas[, tracer_names, drop = FALSE])
  m <- as.matrix(tracer_list$ns[, tracer_names, drop = FALSE])
  storage.mode(m) <- "integer"
  P_sampled <- as.matrix(df_stream_2[, source_order, drop = FALSE])
  if ("Unsampled" %in% names(df_stream_2)) {
    P_unsampled <- as.numeric(df_stream_2[["Unsampled"]])
  } else {
    P_unsampled <- 1 - rowSums(P_sampled)
  }
  # Numerical guard: proportions are renormalised after flooring.
  P_all <- cbind(P_sampled, Unsampled = P_unsampled)
  P_all <- pmax(P_all, eps)
  P_all <- P_all / rowSums(P_all)
  ln_rho <- log(P_all)
  N <- nrow(Y)
  J <- nrow(x)
  K <- ncol(Y)
  if (length(sigma_ln_rho) == 1L) {
    sigma_ln_rho_mat <- matrix(as.numeric(sigma_ln_rho), nrow = N, ncol = J + 1)
  } else if (is.vector(sigma_ln_rho) && length(sigma_ln_rho) == (J + 1L)) {
    sigma_ln_rho_mat <- matrix(
      rep(as.numeric(sigma_ln_rho), each = N), nrow = N, ncol = J + 1,
      byrow = FALSE
    )
  } else {
    sigma_ln_rho_mat <- as.matrix(sigma_ln_rho)
    if (!identical(dim(sigma_ln_rho_mat), c(N, J + 1L))) {
      stop(
        "`sigma_ln_rho` must be a scalar, a vector of length J + 1, ",
        "or an N x (J + 1) matrix.",
        call. = FALSE
      )
    }
  }
  sdata <- list(
    N = N, J = J, K = K, Y = Y, x = x, s = s, m = m, ln_rho = ln_rho,
    sigma_ln_rho = sigma_ln_rho_mat
  )
  if (isTRUE(hierarchical)) {
    group_values <- NULL
    if ("group" %in% names(df_stream_1)) {
      group_values <- df_stream_1[["group"]]
    } else if ("group" %in% names(df_stream_2)) {
      group_values <- df_stream_2[["group"]]
    }
    if (is.null(group_values)) {
      group_values <- rep("group_1", N)
    }
    group_fac <- factor(group_values)
    sdata$R <- nlevels(group_fac)
    sdata$YR <- as.integer(group_fac)
  }
  standardisation_wrapper(
    sdata = sdata,
    fix_unsampled = fix_unsampled,
    unsampled_prior_ridge = unsampled_prior_ridge
  )
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
#' containing the run time (`timing`), the fitted Stan model (`model`), and
#' the standardised Stan data list used to fit it (`data`).
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
    models[[i]] <- list(timing = timing, model = model, data = sdata)
    names(models)[i] <- mod_name_suffix_i
  }
  models
}

#' Wrapper to run mixture model in Stan
#' 
#' @param model_path Character scalar. Path to the `.stan` file to compile and
#'   fit, typically generated by \code{\link{build_stancode}}.
#' @param mod_name_suffix A character string used to create a unique model name.
#' @param ... Additional arguments passed to \code{\link[rstan]{stan}}, such as
#'   `data`, `iter`, `chains`, or `control`.
#'
#' @return A fitted `stanfit` object, as returned by \code{\link[rstan]{stan}}.
#'
#' @importFrom rstan stan
#' 
#' @seealso
#'   \code{\link{run_mixmustr_models}}
#' 
#' @export
run_mixmod <- function(model_path, mod_name_suffix, ...) {
  stan(file = model_path, model_name = mod_name_suffix, ...)
}
