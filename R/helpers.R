#' @importFrom dplyr %>% mutate rowwise c_across
#' @importFrom tidyselect everything
#' @noRd
reshape_ref_data <- function(x, target = "df_stream_2", order_ref) {
  (x[[target]][, order_ref]) |>
    rowwise() |>
    mutate(Unsampled = 1 - sum(c_across(everything()))) |>
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

#' @importFrom truncnorm rtruncnorm
#' @importFrom tidyr replace_na
#' @noRd
trun_na_zr <- function(...) {
  rtruncnorm(1, ...) |>
    replace_na(0)
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
