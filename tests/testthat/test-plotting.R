test_that("plot_multiple_faceted_scatter_avg returns a ggplot", {
  data <- data.frame(
    Observed = c(0.2, 0.4, 0.6, 0.8),
    Predicted = c(0.25, 0.35, 0.65, 0.75),
    ymin = c(0.2, 0.3, 0.6, 0.7),
    ymax = c(0.3, 0.5, 0.7, 0.9),
    `Variant:` = c("A", "B", "A", "B"),
    source = c("Source1", "Source1", "Source2", "Source2"),
    check.names = FALSE
  )
  p <- plot_multiple_faceted_scatter_avg(data)
  expect_s3_class(p, "ggplot")
})

test_that("evaluate_uncertainty errors on mismatched or non-numeric input", {
  expect_error(evaluate_uncertainty(c(0.5, 0.5), 1))
  expect_error(evaluate_uncertainty(c("a", "b"), c(1, 1)))
})

test_that("evaluate_uncertainty returns a list(data, plot) with the expected shape", {
  pi <- c(a = 0.6, b = 0.4)
  out <- evaluate_uncertainty(pi, c(1, 1), iter = 50, seed = 1)
  expect_named(out, c("data", "plot"))
  expect_equal(nrow(out$data), 50 * length(pi))
  expect_setequal(unique(out$data$Sources), names(pi))
  expect_s3_class(out$plot, "ggplot")
})

test_that("evaluate_uncertainty uses generic names when pi is unnamed", {
  out <- evaluate_uncertainty(c(0.6, 0.4), c(1, 1), iter = 10, seed = 1)
  expect_setequal(unique(out$data$Sources), c("Source 1", "Source 2"))
})

test_that("evaluate_uncertainty is reproducible given the same seed", {
  a <- evaluate_uncertainty(c(0.6, 0.4), c(1, 1), iter = 20, seed = 7)
  b <- evaluate_uncertainty(c(0.6, 0.4), c(1, 1), iter = 20, seed = 7)
  expect_identical(a$data, b$data)
})
