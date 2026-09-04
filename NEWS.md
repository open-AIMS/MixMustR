# MixMustR 0.1.0

First versioned release of `MixMustR`.

## New features

* New Stan model-generation engine (`build_stancode()`): baseline-logit simplex
  parametrisation, `log_softmax()` stream-2 likelihood, heteroscedastic
  stream-1 likelihood, and a covariance-informed prior for the estimated
  unsampled-source tracer signature.
* `mixmustr_wrangle_input()` rewritten to match: source/tracer alignment,
  optional `Unsampled` column handling, and `sigma_ln_rho` now accepts a
  scalar, a length-`J + 1` vector, or an `N x (J + 1)` matrix (previously
  `N x J`).
* `run_mixmustr_models()`/`run_mixmod()` output now also returns the
  standardised Stan data list used to fit each model (`$data`).
* `synthetic_df_convergent` and `synthetic_df_divergent` regenerated using a
  hierarchical logistic-normal simulation design (10 groups x 10 observations)
  with an explicit, out-of-hull pooled-unsampled tracer signature, matching
  the factorial simulation study described in the accompanying manuscript.

## Testing

* Added a full `testthat` suite covering every exported and internal
  function.

## Documentation

* Overhauled README, DESCRIPTION, package-level docs, dataset docs, and the
  introductory vignette for accuracy and consistency with the current
  implementation and the accompanying manuscript.

## Bug fixes

* Replaced deprecated `.data$col` usage inside `select()`/`rename()`/
  `pivot_longer()`/`pivot_wider()`/`relocate()` (tidyselect >= 1.2.0) with
  plain strings.
* Replaced deprecated `dplyr::case_match()` with `dplyr::recode_values()`.
* Replaced deprecated `ggplot2::geom_errorbarh()` with
  `ggplot2::geom_errorbar(orientation = "y")`.
* Removed dead/unused internal functions (`check_sigma_ln_rho()`, and the
  Dirichlet-based synthetic-data helpers superseded by the new simulation
  engine).
