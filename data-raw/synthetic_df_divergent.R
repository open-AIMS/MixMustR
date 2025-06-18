library(usethis)
library(mixmustr)
synthetic_data <- data.frame(
  group = paste0("G", c(rep(1, 11), rep(2, 4), rep(3, 18), rep(4, 4), rep(5, 8), rep(6, 8), rep(7, 26)))
)
mus <- tracer_parameters$mus
sources <- mus$source
tracers <- mixmustr:::produce_mix_props(
  synthetic_data, sources, seed_1 = 1, seed_2 = 1, delta = -0.05
)
edna <- mixmustr:::produce_mix_props(
  synthetic_data, sources, seed_1 = 1, seed_2 = 3, delta = -0.05
)
synthetic_df_divergent <- mixmustr:::make_mixture_data(
  bcs_si, bcs_fa, tracers[, -1], edna[, -1], truth_stream = 1,
  synthetic_data, rand_gen = TRUE, sd_ = 1, seed = 3
)
usethis::use_data(synthetic_df_divergent, overwrite = TRUE)
