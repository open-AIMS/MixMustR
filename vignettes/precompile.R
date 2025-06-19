# Adapted from
# https://github.com/bcgov/bcdata/blob/master/vignettes/precompile.R

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0

library(knitr)
library(tools)
library(purrr)

rm(list = ls())

# Convert *.orig to *.Rmd -------------------------------------------------
orig_files <- dir(path = "vignettes/", pattern = "*\\.Rmd\\.orig",
                  full.names = TRUE)
# need to set system variable locally first -------------------------------
Sys.setenv("NOT_CRAN" = "true")
purrr::walk(orig_files, ~knitr::knit(.x, tools::file_path_sans_ext(.x)))
# fix image html tags -------------------
new_files <- tools::file_path_sans_ext(orig_files)
for (i in seq_along(new_files)) {
  tmp <- readLines(new_files[i])
  to_fix <- which(grepl("<embed src=", tmp, fixed = TRUE))
  for (j in seq_along(to_fix)) {
    target <- (strsplit(tmp[to_fix[j]], "\"")[[1]]) |>
      grep("\\.pdf$", x = _, value = TRUE)
    new_tags <- paste0(
      "<p align=\"center\">\n", "  <img src=\"", target, "\" width = 800/>\n",
      "</p>\n"
    )
    tmp[to_fix[j]] <- new_tags
  }
  writeLines(tmp, new_files[i])
}
# Move .pdf files into correct directory so they render -------------------
images <- dir(".", pattern = "vignette-fig.*\\.pdf$")
success <- file.copy(from = images, to = file.path("vignettes", images),
                     overwrite = TRUE)
# Clean up if successful --------------------------------------------------
if (!all(success)) {
  stop("Image files were not successfully transferred to vignettes directory")
} else {
  unlink(images)
}
