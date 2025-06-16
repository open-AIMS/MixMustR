library(usethis)
library(mixmustr)
tracer_parameters <- wrangle_tracer_pars(bcs_si, bcs_fa)
usethis::use_data(tracer_parameters, overwrite = TRUE)
