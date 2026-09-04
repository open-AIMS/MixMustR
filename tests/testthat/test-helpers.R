test_that("reshape_ref_data returns an N x (J + 1) matrix summing to 1", {
  m <- reshape_ref_data(
    synthetic_df_convergent, target = "df_stream_2",
    order_ref = tracer_parameters$mus$source
  )
  expect_true(is.matrix(m))
  expect_equal(nrow(m), nrow(synthetic_df_convergent$df_stream_2))
  expect_equal(ncol(m), length(tracer_parameters$mus$source) + 1L)
  expect_equal(colnames(m)[ncol(m)], "Unsampled")
  expect_equal(unname(rowSums(m)), rep(1, nrow(m)), tolerance = 1e-8)
})

test_that("reshape_ref_data works with the stream_1_props target", {
  m <- reshape_ref_data(
    synthetic_df_convergent, target = "stream_1_props",
    order_ref = tracer_parameters$mus$source
  )
  expect_equal(ncol(m), length(tracer_parameters$mus$source) + 1L)
  expect_equal(unname(rowSums(m)), rep(1, nrow(m)), tolerance = 1e-8)
})

test_that("reshape_ref_data errors when order_ref has an unknown source", {
  expect_error(
    reshape_ref_data(
      synthetic_df_convergent, target = "df_stream_2",
      order_ref = c(tracer_parameters$mus$source, "NotASource")
    )
  )
})

test_that("trun_na_zr draws within bounds and is reproducible", {
  set.seed(1)
  x <- MixMustR:::trun_na_zr(a = 0, b = 1, mean = 0.5, sd = 0.1)
  expect_length(x, 1)
  expect_true(x >= 0 && x <= 1)
})

test_that("trun_na_zr replaces an invalid (NaN) draw with 0", {
  expect_equal(MixMustR:::trun_na_zr(a = 5, b = 1, mean = 0, sd = 1), 0)
})

test_that("assign_new_names renames columns in place", {
  df <- data.frame(a = 1, b = 2)
  out <- MixMustR:::assign_new_names(df, c("x", "y"))
  expect_equal(names(out), c("x", "y"))
  expect_equal(out$x, 1)
})

test_that("assign_new_names errors on a length mismatch", {
  df <- data.frame(a = 1, b = 2)
  expect_error(MixMustR:::assign_new_names(df, c("x", "y", "z")))
})
