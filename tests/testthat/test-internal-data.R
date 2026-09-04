test_that("simulate_mvn_mixture returns valid simplex rows", {
  out <- MixMustR:::simulate_mvn_mixture(G = 3, M = 4, J = 3, noise = 0.5, seed = 1)
  expect_named(out, c("dataset", "group_means"))
  expect_equal(nrow(out$dataset), 12)
  expect_equal(ncol(out$dataset), 4)
  expect_equal(names(out$dataset), c("group", "Source_1", "Source_2", "Source_3"))
  props <- as.matrix(out$dataset[, -1])
  expect_equal(unname(rowSums(props)), rep(1, nrow(props)), tolerance = 1e-8)
  expect_true(all(props >= 0 & props <= 1))
  expect_equal(nrow(out$group_means), 3)
})

test_that("simulate_mvn_mixture is reproducible given the same seed", {
  a <- MixMustR:::simulate_mvn_mixture(G = 2, M = 2, J = 2, seed = 42)
  b <- MixMustR:::simulate_mvn_mixture(G = 2, M = 2, J = 2, seed = 42)
  expect_identical(a$dataset, b$dataset)
})

test_that("simulate_mvn_mixture differs across seeds", {
  a <- MixMustR:::simulate_mvn_mixture(G = 2, M = 2, J = 2, seed = 1)
  b <- MixMustR:::simulate_mvn_mixture(G = 2, M = 2, J = 2, seed = 2)
  expect_false(isTRUE(all.equal(a$dataset, b$dataset)))
})

test_that("rename_sources maps Source_N columns by numeric suffix", {
  df <- data.frame(group = 1, Source_2 = 0.4, Source_1 = 0.6)
  out <- MixMustR:::rename_sources(df, c("A", "B"))
  expect_equal(out$A, 0.6)
  expect_equal(out$B, 0.4)
})

test_that("rename_sources errors on a length mismatch", {
  df <- data.frame(group = 1, Source_1 = 0.6, Source_2 = 0.4)
  expect_error(MixMustR:::rename_sources(df, "A"))
})

test_that("rm_unsampled strips the Unsampled column only from the two relevant elements", {
  x <- list(
    df_stream_1 = data.frame(z = 1),
    df_stream_2 = data.frame(a = 1, Unsampled = 2),
    stream_1_props = data.frame(b = 3, Unsampled = 4)
  )
  out <- MixMustR:::rm_unsampled(x)
  expect_false("Unsampled" %in% names(out$df_stream_2))
  expect_false("Unsampled" %in% names(out$stream_1_props))
  expect_equal(out$df_stream_1, x$df_stream_1)
})

test_that("make_unsampled_signature centres exactly on the sampled-source mean when delta = 0", {
  mus_tab <- tracer_parameters$mus
  sig <- MixMustR:::make_unsampled_signature(mus_tab, delta = 0)
  expect_equal(sig$source, "Unsampled")
  expect_equal(names(sig), names(mus_tab))
  tracer_cols <- setdiff(names(mus_tab), "source")
  expect_equal(
    as.numeric(sig[, tracer_cols]),
    unname(colMeans(mus_tab[, tracer_cols])),
    tolerance = 1e-8
  )
})

test_that("make_unsampled_signature displaces by exactly delta scaled units", {
  mus_tab <- tracer_parameters$mus
  tracer_cols <- setdiff(names(mus_tab), "source")
  s_scaled <- scale(as.matrix(mus_tab[, tracer_cols]))
  center <- attr(s_scaled, "scaled:center")
  scal <- attr(s_scaled, "scaled:scale")
  sig <- MixMustR:::make_unsampled_signature(mus_tab, delta = 3)
  displacement_scaled <- (as.numeric(sig[, tracer_cols]) - center) / scal
  expect_equal(sqrt(sum(displacement_scaled^2)), 3, tolerance = 1e-6)
})

test_that("augment_bcs_with_unsampled appends a homogenised Unsampled row", {
  mus_tab <- tracer_parameters$mus
  out <- MixMustR:::augment_bcs_with_unsampled(bcs_si, bcs_fa, mus_tab, delta = 3)
  expect_named(out, c("bcs_si", "bcs_fa"))
  expect_equal(nrow(out$bcs_si), nrow(bcs_si) + 1L)
  expect_equal(nrow(out$bcs_fa), nrow(bcs_fa) + 1L)
  expect_equal(out$bcs_si$source[nrow(out$bcs_si)], "Unsampled")
  expect_true("Taxa" %in% names(out$bcs_fa))
  expect_equal(names(out$bcs_fa)[2], "Taxa")
  sd_cols_si <- grep("\\(SD\\)$", names(out$bcs_si), value = TRUE)
  n_cols_si <- grep("\\(n\\)$", names(out$bcs_si), value = TRUE)
  expect_true(all(unlist(out$bcs_si[, sd_cols_si]) == 1))
  expect_true(all(unlist(out$bcs_si[, n_cols_si]) == 10L))
  sd_cols_fa <- grep("\\(SD\\)$", names(out$bcs_fa), value = TRUE)
  n_cols_fa <- grep("\\(n\\)$", names(out$bcs_fa), value = TRUE)
  expect_true(all(unlist(out$bcs_fa[, sd_cols_fa]) == 1))
  expect_true(all(unlist(out$bcs_fa[, n_cols_fa]) == 10L))
})

test_that("make_mixture_data combines signatures and proportions correctly", {
  mus_tab <- tracer_parameters$mus
  source_names <- c(mus_tab$source, "Unsampled")
  signatures <- MixMustR:::augment_bcs_with_unsampled(bcs_si, bcs_fa, mus_tab, delta = 3)
  props_1 <- MixMustR:::simulate_mvn_mixture(
    G = 2, M = 2, J = length(source_names), noise = 0.5, seed = 11
  )$dataset |>
    MixMustR:::rename_sources(source_names)
  props_2 <- MixMustR:::simulate_mvn_mixture(
    G = 2, M = 2, J = length(source_names), noise = 1, seed = 11
  )$dataset |>
    MixMustR:::rename_sources(source_names)
  template_df <- data.frame(group = paste0("G", rep(1:2, each = 2)))
  out <- MixMustR:::make_mixture_data(
    signatures$bcs_si, signatures$bcs_fa, props_1[, -1], props_2[, -1],
    truth_stream = 1, template_df, rand_gen = TRUE, sd_ = 1, seed = 1
  )
  expect_named(out, c("df_stream_1", "df_stream_2", "stream_1_props"))
  expect_equal(nrow(out$df_stream_1), 4)
  expect_equal(ncol(out$df_stream_1), 1 + length(setdiff(names(mus_tab), "source")))
  expect_equal(names(out$df_stream_2), c("group", source_names))
})

test_that("make_mixture_data errors on an invalid truth_stream", {
  mus_tab <- tracer_parameters$mus
  source_names <- c(mus_tab$source, "Unsampled")
  signatures <- MixMustR:::augment_bcs_with_unsampled(bcs_si, bcs_fa, mus_tab, delta = 3)
  props <- MixMustR:::simulate_mvn_mixture(
    G = 2, M = 2, J = length(source_names), seed = 11
  )$dataset |>
    MixMustR:::rename_sources(source_names)
  template_df <- data.frame(group = paste0("G", rep(1:2, each = 2)))
  expect_error(
    MixMustR:::make_mixture_data(
      signatures$bcs_si, signatures$bcs_fa, props[, -1], props[, -1],
      truth_stream = 3, template_df
    ),
    "truth_stream"
  )
})

test_that("calc_tracer_estimate returns the raw mean when rand_gen = FALSE", {
  x <- data.frame(
    source = c("A", "B"), tracer_family = "si", tracer = "t1",
    mean = c(1, 2), sd = c(0.1, 0.1), a = -Inf, b = Inf
  )
  out <- MixMustR:::calc_tracer_estimate(x, rand_gen = FALSE)
  expect_equal(out$estimate, c(1, 2))
})

test_that("calc_tracer_estimate is reproducible when rand_gen = TRUE", {
  x <- data.frame(
    source = c("A", "B"), tracer_family = "si", tracer = "t1",
    mean = c(1, 2), sd = c(0.1, 0.1), a = -Inf, b = Inf
  )
  out1 <- MixMustR:::calc_tracer_estimate(x, rand_gen = TRUE, sd_ = 0.1, seed = 5)
  out2 <- MixMustR:::calc_tracer_estimate(x, rand_gen = TRUE, sd_ = 0.1, seed = 5)
  expect_equal(out1, out2)
})

test_that("wrangle_tracer_pars reproduces the built-in tracer_parameters dataset", {
  rebuilt <- MixMustR:::wrangle_tracer_pars(bcs_si, bcs_fa)
  expect_identical(rebuilt, tracer_parameters)
})

test_that("reshape_isotope_df pivots to the expected long format", {
  out <- MixMustR:::reshape_isotope_df(bcs_si)
  expect_true(all(c("source", "tracer", "mean", "sd", "n", "tracer_family") %in% names(out)))
  mean_cols <- grep("^d\\(", names(bcs_si), value = TRUE)
  mean_cols <- mean_cols[!grepl("\\(SD\\)$|\\(n\\)$", mean_cols)]
  expect_equal(nrow(out), nrow(bcs_si) * length(mean_cols))
  expect_true(all(out$tracer_family == "si"))
})

test_that("reshape_fattyacids_df pivots to the expected long format", {
  out <- MixMustR:::reshape_fattyacids_df(bcs_fa)
  expect_true(all(c("source", "tracer", "mean", "sd", "n", "tracer_family") %in% names(out)))
  expect_true(all(out$tracer_family == "fa"))
  expect_false("Taxa" %in% names(out))
})

test_that("compare_mixing_proportions returns a ggplot", {
  p <- compare_mixing_proportions(
    synthetic_df_divergent, synthetic_df_convergent, tracer_parameters$mus
  )
  expect_s3_class(p, "ggplot")
})
