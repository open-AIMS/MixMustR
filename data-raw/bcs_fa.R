library(usethis)
library(dplyr)
bcs_fa <- read.csv("data-raw/bcs_fa.csv", header = TRUE, check.names = FALSE) |>
  dplyr::arrange(source)
usethis::use_data(bcs_fa, overwrite = TRUE)
