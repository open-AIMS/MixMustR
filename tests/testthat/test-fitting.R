test_that("build_stancode returns a single character string", {
  code <- build_stancode(sample_tracer = FALSE, fix_unsampled = TRUE, hierarchical = FALSE)
  expect_type(code, "character")
  expect_length(code, 1)
})

test_that("build_stancode toggles hierarchical Stan data/parameters", {
  hier <- build_stancode(FALSE, TRUE, TRUE)
  indep <- build_stancode(FALSE, TRUE, FALSE)
  expect_match(hier, "int<lower=1> R;", fixed = TRUE)
  expect_match(hier, "array\\[N\\] int<lower=1, upper=R> YR;")
  expect_false(grepl("int<lower=1> R;", indep, fixed = TRUE))
})

test_that("build_stancode toggles the fixed vs estimated unsampled signature", {
  fixed <- build_stancode(FALSE, TRUE, FALSE)
  estimated <- build_stancode(FALSE, FALSE, FALSE)
  expect_match(fixed, "Fixed unsampled-source signature", fixed = TRUE)
  expect_match(estimated, "x_bar_unsampled_prior", fixed = TRUE)
  expect_false(grepl("x_bar_unsampled_prior", fixed, fixed = TRUE))
})

test_that("build_stancode toggles estimated sampled-source signatures", {
  estimated <- build_stancode(TRUE, TRUE, FALSE)
  fixed <- build_stancode(FALSE, TRUE, FALSE)
  expect_match(estimated, "mu_source[j, k]", fixed = TRUE)
  expect_false(grepl("mu_source[j, k]", fixed, fixed = TRUE))
})

test_that("build_stancode writes the code to code_path when requested", {
  stan_file <- tempfile(fileext = ".stan")
  code <- build_stancode(FALSE, TRUE, FALSE, code_path = stan_file)
  expect_true(file.exists(stan_file))
  expect_equal(readLines(stan_file), strsplit(code, "\n")[[1]])
})

test_that("build_stancode errors when write_file = TRUE without a code_path", {
  expect_error(
    build_stancode(FALSE, TRUE, FALSE, code_path = NULL, write_file = TRUE),
    "code_path"
  )
})

test_that("standardise_tracer_sdata centres and scales x to mean 0 / sd 1", {
  sdata <- list(
    Y = matrix(1:6, nrow = 3, ncol = 2),
    x = matrix(c(1, 2, 3, 10, 20, 30), nrow = 3, ncol = 2)
  )
  out <- MixMustR:::standardise_tracer_sdata(sdata)
  expect_equal(unname(colMeans(out$x)), c(0, 0), tolerance = 1e-8)
  expect_equal(unname(apply(out$x, 2, sd)), c(1, 1), tolerance = 1e-8)
  expect_equal(attr(out, "tracer_center"), unname(colMeans(sdata$x)), ignore_attr = TRUE)
})

test_that("standardise_tracer_sdata replaces a zero/non-finite scale with 1", {
  sdata <- list(Y = matrix(1:4, 2, 2), x = matrix(c(5, 5, 1, 2), 2, 2))
  out <- MixMustR:::standardise_tracer_sdata(sdata)
  expect_equal(attr(out, "tracer_scale")[1], 1, ignore_attr = TRUE)
})

test_that("standardise_tracer_sdata errors when Y or x is missing", {
  expect_error(MixMustR:::standardise_tracer_sdata(list(x = matrix(1, 1, 1))), "`Y`")
  expect_error(MixMustR:::standardise_tracer_sdata(list(Y = matrix(1, 1, 1))), "`x`")
})

test_that("make_unsampled_prior_cov returns a consistent mean and Cholesky factor", {
  sdata <- list(x = matrix(c(1, -1, 0, 2, -2, 0, 0.5, -0.5, 1), nrow = 3, ncol = 3))
  out <- MixMustR:::make_unsampled_prior_cov(sdata, ridge = 0.1)
  expect_equal(out$x_bar, unname(colMeans(sdata$x)), ignore_attr = TRUE)
  reconstructed <- out$L_x %*% t(out$L_x)
  expect_equal(
    reconstructed, cov(sdata$x) + diag(0.1^2, ncol(sdata$x)),
    tolerance = 1e-8, ignore_attr = TRUE
  )
})

test_that("standardisation_wrapper only adds unsampled-prior data when fix_unsampled = FALSE", {
  sdata <- list(Y = matrix(1:6, 3, 2), x = matrix(c(1, 2, 3, 10, 20, 30), 3, 2))
  fixed <- MixMustR:::standardisation_wrapper(sdata, fix_unsampled = TRUE)
  estimated <- MixMustR:::standardisation_wrapper(sdata, fix_unsampled = FALSE)
  expect_false("x_bar_unsampled_prior" %in% names(fixed))
  expect_true(all(c("x_bar_unsampled_prior", "L_x_unsampled_prior") %in% names(estimated)))
  expect_length(estimated$x_bar_unsampled_prior, ncol(sdata$x))
  expect_equal(dim(estimated$L_x_unsampled_prior), c(ncol(sdata$x), ncol(sdata$x)))
})

test_that("check_sd_tabs passes for the real tracer_parameters tables", {
  expect_silent(MixMustR:::check_sd_tabs(tracer_parameters$sigmas, tracer_parameters$mus, "SDs"))
  expect_silent(MixMustR:::check_sd_tabs(tracer_parameters$ns, tracer_parameters$mus, "ns"))
})

test_that("check_sd_tabs errors when inputs are not data frames", {
  expect_error(MixMustR:::check_sd_tabs(1, tracer_parameters$mus))
})

test_that("check_sd_tabs errors on mismatched structure", {
  bad <- tracer_parameters$sigmas[, -2]
  expect_error(MixMustR:::check_sd_tabs(bad, tracer_parameters$mus), "do not match")
})

test_that("check_sd_tabs errors when the first column isn't named 'source'", {
  a <- data.frame(src = "A", t1 = 1.0)
  b <- data.frame(src = "A", t1 = 1.0)
  expect_error(MixMustR:::check_sd_tabs(a, b), "First column name")
})

test_that("check_sd_tabs errors on NA or wrong-type values", {
  bad <- tracer_parameters$sigmas
  bad[1, 2] <- NA
  expect_error(MixMustR:::check_sd_tabs(bad, tracer_parameters$mus), "NAs")
})

test_that("mixmustr_wrangle_input produces a correctly-dimensioned Stan data list", {
  sdata <- mixmustr_wrangle_input(
    synthetic_df_convergent, tracer_parameters, sigma_ln_rho = 1,
    sample_tracer = FALSE, fix_unsampled = TRUE, hierarchical = TRUE
  )
  expect_equal(sdata$N, nrow(synthetic_df_convergent$df_stream_1))
  expect_equal(sdata$J, length(tracer_parameters$mus$source))
  expect_equal(sdata$K, ncol(tracer_parameters$mus) - 1L)
  expect_equal(sdata$R, length(unique(synthetic_df_convergent$df_stream_1$group)))
  expect_length(sdata$YR, sdata$N)
  expect_equal(dim(sdata$sigma_ln_rho), c(sdata$N, sdata$J + 1L))
  expect_true(all(sdata$sigma_ln_rho == 1))
})

test_that("mixmustr_wrangle_input omits R/YR when hierarchical = FALSE", {
  sdata <- mixmustr_wrangle_input(
    synthetic_df_convergent, tracer_parameters, sigma_ln_rho = 1,
    sample_tracer = FALSE, fix_unsampled = TRUE, hierarchical = FALSE
  )
  expect_false(any(c("R", "YR") %in% names(sdata)))
})

test_that("mixmustr_wrangle_input adds the unsampled-prior data only when fix_unsampled = FALSE", {
  fixed <- mixmustr_wrangle_input(
    synthetic_df_convergent, tracer_parameters, sigma_ln_rho = 1,
    sample_tracer = FALSE, fix_unsampled = TRUE, hierarchical = FALSE
  )
  estimated <- mixmustr_wrangle_input(
    synthetic_df_convergent, tracer_parameters, sigma_ln_rho = 1,
    sample_tracer = FALSE, fix_unsampled = FALSE, hierarchical = FALSE
  )
  expect_false("x_bar_unsampled_prior" %in% names(fixed))
  expect_true(all(c("x_bar_unsampled_prior", "L_x_unsampled_prior") %in% names(estimated)))
})

test_that("mixmustr_wrangle_input accepts a length J + 1 sigma_ln_rho vector", {
  j <- length(tracer_parameters$mus$source)
  sigma_vec <- seq_len(j + 1)
  sdata <- mixmustr_wrangle_input(
    synthetic_df_convergent, tracer_parameters, sigma_ln_rho = sigma_vec,
    sample_tracer = FALSE, fix_unsampled = TRUE, hierarchical = FALSE
  )
  expect_equal(unname(sdata$sigma_ln_rho[1, ]), as.numeric(sigma_vec))
  expect_equal(unname(sdata$sigma_ln_rho[, 1]), rep(sigma_vec[1], sdata$N))
})

test_that("mixmustr_wrangle_input errors on a sigma_ln_rho of the wrong length", {
  expect_error(
    mixmustr_wrangle_input(
      synthetic_df_convergent, tracer_parameters,
      sigma_ln_rho = seq_len(length(tracer_parameters$mus$source)),
      sample_tracer = FALSE, fix_unsampled = TRUE, hierarchical = FALSE
    ),
    "sigma_ln_rho"
  )
})

test_that("mixmustr_wrangle_input reconstructs Unsampled identically to an explicit column", {
  implicit <- mixmustr_wrangle_input(
    synthetic_df_convergent, tracer_parameters, sigma_ln_rho = 1,
    sample_tracer = FALSE, fix_unsampled = TRUE, hierarchical = FALSE
  )
  explicit_streams <- synthetic_df_convergent
  sampled <- explicit_streams$df_stream_2[, tracer_parameters$mus$source]
  explicit_streams$df_stream_2$Unsampled <- 1 - rowSums(sampled)
  explicit <- mixmustr_wrangle_input(
    explicit_streams, tracer_parameters, sigma_ln_rho = 1,
    sample_tracer = FALSE, fix_unsampled = TRUE, hierarchical = FALSE
  )
  expect_equal(implicit$ln_rho, explicit$ln_rho)
})

test_that("mixmustr_wrangle_input errors on missing required elements", {
  expect_error(
    mixmustr_wrangle_input(
      list(df_stream_2 = synthetic_df_convergent$df_stream_2), tracer_parameters,
      sigma_ln_rho = 1, sample_tracer = FALSE, fix_unsampled = TRUE, hierarchical = FALSE
    ),
    "df_stream_1"
  )
  expect_error(
    mixmustr_wrangle_input(
      synthetic_df_convergent, list(mus = tracer_parameters$mus),
      sigma_ln_rho = 1, sample_tracer = FALSE, fix_unsampled = TRUE, hierarchical = FALSE
    ),
    "sigmas"
  )
})

test_that("mixmustr_wrangle_input errors when tracer or source columns are missing", {
  bad_streams <- synthetic_df_convergent
  bad_streams$df_stream_1 <- bad_streams$df_stream_1[, "group", drop = FALSE]
  expect_error(
    mixmustr_wrangle_input(
      bad_streams, tracer_parameters, sigma_ln_rho = 1,
      sample_tracer = FALSE, fix_unsampled = TRUE, hierarchical = FALSE
    ),
    "tracer column"
  )
  bad_streams2 <- synthetic_df_convergent
  bad_streams2$df_stream_2 <- bad_streams2$df_stream_2[, "group", drop = FALSE]
  expect_error(
    mixmustr_wrangle_input(
      bad_streams2, tracer_parameters, sigma_ln_rho = 1,
      sample_tracer = FALSE, fix_unsampled = TRUE, hierarchical = FALSE
    ),
    "sampled-source column"
  )
})

test_that("run_mixmustr_models and run_mixmod return a usable stanfit", {
  skip_on_cran()
  skip_if_no_model_fit()
  expect_named(minimal_model_fit, c("timing", "model", "data"))
  expect_s4_class(minimal_model_fit$model, "stanfit")

  model_path <- mixmustr_models[6, ]$code_path
  fit2 <- run_mixmod(
    model_path, "test_run_mixmod", data = minimal_model_fit$data,
    iter = 40, warmup = 20, chains = 1, cores = 1, refresh = 0
  )
  expect_s4_class(fit2, "stanfit")
})
