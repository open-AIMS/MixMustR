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
      "\tint<lower=1> m[J]; // source-specific sample size\n"
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
      "\tint<lower=1> YR[N]; // dummy vector of random levels\n"
    )
  }
  script <- c(script, "}\n", "parameters {\n",
    "\tvector[J + 1] zeta[N]; // zeta from softmax\n"
  )
  if (hierarchical) {
    script <- c(script,
      "\tvector[J + 1] randvec[R]; // random effects\n",
      "\treal<lower=0> randsd[J + 1]; // random effects hyper parameter\n",
      "\treal<lower=0> ext_sigma[J + 1]; // residual error, likelihood 2\n"
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
      "\treal<lower=0> sigma[K]; // residual SDs for each tracer\n"
    )
  }
  if (!sample_tracer & !fix_unsampled) {
    script <- c(script,
      "\treal sampled_mean_x[K]; // mean unsampled estimate for each tracer\n"
    )
  }
  script <- c(script, "}\n", "transformed parameters {\n",
    "\tsimplex[J + 1] p[N]; // probability vector for each site\n"
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
      "\t\t\tomega[j, k] = sqrt((s[j, k]^2 * (m[j] - 1)) / nu[j, k]);\n",
      "\t\t}\n", "\t}\n"
    )
  }
  script <- c(script, "}\n", "model {\n")
  if (sample_tracer) {
    script <- c(script,
      "\t// Sampling sources\n",
      "\tfor (j in 1:J) {\n",
      "\t\tfor (k in 1:K) {\n",
      "\t\t\ttarget += chi_square_lpdf(m[j] | nu[j, k]);\n",
      "\t\t\ttarget += normal_lpdf(x[j, k] | mu[j, k], s[j, k] / sqrt(m[j]));\n",
      "\t\t}\n",
      "\t}\n"
    )
  }
  sd_term <- ifelse(sample_tracer, "sqrt(p[n][j]^2 * omega[j, k]^2)", "sigma[k]")
  script <- c(script,
    "\t// Mixture model\n",
    "\tfor (n in 1:N) {\n",
    "\t\tfor (k in 1:K) {\n",
    "\t\t\treal phi[J + 1];\n",
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

#' @export
mixmustr_wrangle_input <- function(synth_df, mu_tab, sig_tab, model_path,
                                   sigma_ln_rho, m, sample_tracer,
                                   fix_unsampled, hierarchical, ...) {
  yobs <- synth_df$sim_data |>
    dplyr::select(dplyr::select(sig_tab, -Group) |> names())
  reff_df <- synth_df$sim_data |>
    dplyr::mutate(YR = as.numeric(as.factor(group)))
  out <- list(
    N = nrow(yobs),
    J = nrow(mu_tab), # sources
    K = ncol(yobs), # tracers
    Y = yobs,
    x = dplyr::select(mu_tab, -Group),
    ln_rho = reshape_ref_data(
      synth_df, target = "stream_2_df", order_ref = mu_tab$Group
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
  if (sample_tracer) {
    out[["s"]] <- dplyr::select(sig_tab, -Group)
    out[["m"]] <- rep(m, nrow(mu_tab))
  }
  if (!fix_unsampled & !sample_tracer) {
    out[["s"]] <- dplyr::select(sig_tab, -Group)
  }
  if (hierarchical) {
    out[["R"]] <- max(reff_df$YR)
    out[["YR"]] <- reff_df$YR
  }
  out
}

#' @export
run_all_mixmustr_models <- function(model_choices, synth_df, mu_tab, sig_tab,
                                    sigma_ln_rho, ...) {
  models <- vector(mode = "list", length = nrow(model_choices))
  for (i in seq_len(nrow(model_choices))) {
    build_stancode(model_choices$sample_tracer[i],
                   model_choices$fix_unsampled[i],
                   model_choices$hierarchical[i],
                   model_choices$code_path[i])
    model_path <- model_choices$code_path[i]
    sdata <- mixmustr_wrangle_input(
      synth_df, mu_tab, sig_tab, model_path, sigma_ln_rho = sigma_ln_rho,
      m = 100, sample_tracer = model_choices$sample_tracer[i],
      fix_unsampled = model_choices$fix_unsampled[i],
      hierarchical = model_choices$hierarchical[i]
    )
    synth_n_i <- tools::file_path_sans_ext(basename(model_path))
    timing <- system.time({
      model <- run_mixmod(model_path, synth_n_i, data = sdata, ...)
    })
    models[[i]] <- list(timing = timing, model = model)
  }
  models
}

#' @export
run_mixmod <- function(model_path, synth_n, ...) {
  model_name <- basename(model_path) |>
    tools::file_path_sans_ext() |>
    paste0("_", synth_n)
  rstan::stan(file = model_path, model_name = model_name, ...)
}
