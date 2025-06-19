library(usethis)
library(dplyr)
bcs_si <- read.csv("data-raw/bcs_si.csv", header = TRUE, check.names = FALSE) |>
  dplyr::arrange(source)
usethis::use_data(bcs_si, overwrite = TRUE)
