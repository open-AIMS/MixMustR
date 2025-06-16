#' @importFrom dplyr %>% mutate
#' @noRd
reshape_ref_data <- function(x, target = "stream_1_df", order_ref) {
  (x[[target]][, order_ref]) %>%
    mutate(Unsampled = 1 - rowSums(.)) |>
    as.matrix() |>
    abs()
}

#' @noRd
fix_sum_to_one <- function(x) {
  new_x <- x
  to_change <- round(1 - rowSums(new_x), 2)
  for (i in seq_len(nrow(new_x))) {
    possible <- which(new_x[i, ] >= abs(to_change[i]))
    j <- sample(possible, 1)
    tested <- new_x[i, j] + to_change[i]
    new_sum <- sum(c(tested, new_x[i, -j]))
    n_ <- 0
    while (tested < 0 | tested > 1 & new_sum != 1 & n_ <= 50) {
      n_ <- n_ + 1
      possible <- which(new_x[i, ] >= abs(to_change[i]))
      j <- sample(possible, 1)
      tested <- new_x[i, j] + to_change[i]
      new_sum <- sum(c(tested, new_x[i, -j]))
    }
    new_x[i, j] <- tested
  }
  if (!all(round(rowSums(new_x), 2) == 1)) {
    stop("Adjusted row sums do not add up to 1.")
  }
  new_x
}

#' @noRd
linear_rescale <- function(x, r_out) {
  p <- (x - min(x)) / (max(x) - min(x))
  r_out[[1]] + p * (r_out[[2]] - r_out[[1]])
}

#' @noRd
trun_na_zr <- function(...) {
  truncnorm::rtruncnorm(1, ...) |>
    tidyr::replace_na(0)
}

#' @noRd
assign_new_names <- function(df_, new_names) {
  names(df_) <- new_names
  df_
}

#' @noRd
abs_log <- function(x, adjust = 0.001, constant = 0) {
  x <- abs(x)
  x[x <= 0] <- 0 + adjust
  x[x >= 1] <- 1 - adjust
  log(x) + constant
}
