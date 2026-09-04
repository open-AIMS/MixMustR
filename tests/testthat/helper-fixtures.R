# Shared fixtures for tests. Kept intentionally small/low-iteration: these
# tests check structure and correctness of shapes, not inference quality.

small_streams <- synthetic_df_convergent
small_streams$df_stream_1 <- small_streams$df_stream_1[1:20, ]
small_streams$df_stream_2 <- small_streams$df_stream_2[1:20, ]
small_streams$stream_1_props <- small_streams$stream_1_props[1:20, ]

# Single minimally-iterated fit (fastest variant: sample_tracer = FALSE,
# fix_unsampled = TRUE, hierarchical = FALSE), reused across test files that
# need a real `stanfit` object. Graceful degradation if Stan can't compile.
minimal_model_fit <- tryCatch(
  suppressWarnings(
    run_mixmustr_models(
      mixmustr_models[6, ], small_streams, tracer_parameters,
      sigma_ln_rho = 1, iter = 60, warmup = 30, chains = 1, cores = 1,
      refresh = 0
    )[[1]]
  ),
  error = function(e) NULL
)

skip_if_no_model_fit <- function() {
  testthat::skip_if(
    is.null(minimal_model_fit),
    "Could not fit a minimal Stan model in this environment."
  )
}
