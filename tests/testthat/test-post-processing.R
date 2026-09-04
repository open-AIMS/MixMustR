test_that("make_post_prop_long returns the expected long-format columns", {
  skip_on_cran()
  skip_if_no_model_fit()
  out <- make_post_prop_long(
    minimal_model_fit$model, tracer_parameters$mus, small_streams,
    target = "df_stream_2", n = 1
  )
  expect_true(all(
    c("N", "source", "Predicted", "ymin", "ymax", "Observed", "Variant:") %in%
      names(out)
  ))
  expect_true(all(out$`Variant:` == "1"))
  expect_setequal(unique(out$source), c(tracer_parameters$mus$source, "Unsampled"))
})

test_that("mixmustr_bayes_R2 returns a summary table with one row per source", {
  skip_on_cran()
  skip_if_no_model_fit()
  out <- mixmustr_bayes_R2(
    minimal_model_fit$model, summary = TRUE,
    data_streams_list = small_streams, target = "df_stream_2",
    order_ref = tracer_parameters$mus$source
  )
  expect_true(all(c("Source", "mean", "2.5%HDI", "97.5%HDI") %in% names(out)))
  expect_equal(nrow(out), length(tracer_parameters$mus$source) + 1L)
})

test_that("mixmustr_bayes_R2 returns raw draws when summary = FALSE", {
  skip_on_cran()
  skip_if_no_model_fit()
  out <- mixmustr_bayes_R2(
    minimal_model_fit$model, summary = FALSE,
    data_streams_list = small_streams, target = "df_stream_2",
    order_ref = tracer_parameters$mus$source
  )
  expect_true(is.matrix(out))
  expect_equal(ncol(out), length(tracer_parameters$mus$source) + 1L)
})
