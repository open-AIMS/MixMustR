library(usethis)
library(dplyr)
mixmustr_models <- expand.grid(
  sample_tracer = c(TRUE, FALSE), fix_unsampled = c (TRUE, FALSE),
  hierarchical = c(TRUE, FALSE)
) |>
  dplyr::mutate(
    code_path = file.path("stan", paste0("sample_tracer_", sample_tracer,
      "_fix_unsampled_", fix_unsampled,
      "_hierarchical_", hierarchical, ".stan"))
  )
usethis::use_data(mixmustr_models, overwrite = TRUE)
