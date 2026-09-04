library(usethis)
library(MixMustR)
# Matches the 6-sampled-source, 8-tracer, divergent (sim_id = 8) row of the
# factorial simulation grid described in the accompanying manuscript (see
# MixMustR_paper_original/copilot/jobs_dev.csv).
mus <- tracer_parameters$mus
source_names <- c(mus$source, "Unsampled")
j <- length(source_names)
synthetic_data <- data.frame(group = paste0("G", rep(1:10, each = 10)))
signatures <- MixMustR:::augment_bcs_with_unsampled(bcs_si, bcs_fa, mus, delta = 3)
# Stream 1 (tracers) and stream 2 (e.g. eDNA) share group centres (same seed)
# but differ in within-group noise: low noise (0.5) vs high noise (0.001).
df_props_1 <- MixMustR:::simulate_mvn_mixture(
  G = 10, M = 10, J = j, noise = 0.5, seed = 10008
)$dataset |>
  dplyr::mutate(group = paste0("G", group)) |>
  MixMustR:::rename_sources(source_names)
df_props_2 <- MixMustR:::simulate_mvn_mixture(
  G = 10, M = 10, J = j, noise = 0.001, seed = 10008
)$dataset |>
  dplyr::mutate(group = paste0("G", group)) |>
  MixMustR:::rename_sources(source_names)
synthetic_df_divergent <- MixMustR:::make_mixture_data(
  signatures$bcs_si, signatures$bcs_fa, df_props_1[, -1], df_props_2[, -1],
  truth_stream = 1, synthetic_data, rand_gen = TRUE, sd_ = 1, seed = 20008
) |>
  MixMustR:::rm_unsampled()
usethis::use_data(synthetic_df_divergent, overwrite = TRUE)
