library(usethis)
library(dplyr)
bcs_si <- "
source,d(13C/12C),d13C sd,d(15N/14N),d15N sd,Study
Mangrove,-27.7,1.4,4.4,2.6,Nyunja et al. (2009)
Seagrass,-15.1,3,1.2,1.3,Nyunja et al. (2009)
Saltmarsh,-24.6,4.32,2.31,1.67,Bulmer et al. (2020)
Macroalgae,-20.4,3.1,2.4,1.2,Nyunja et al. (2009)
Terrestrial grass,-13.3,0.5,1.69,0.5,Goni & Thomas (2000) Spartina
Plankton,-19.3,2.8,6.2,3,Nyunja et al. (2009)
" |>
  read.table(text = _, header = TRUE, check.names = FALSE, sep = ",") |>
  dplyr::arrange(source)
usethis::use_data(bcs_si, overwrite = TRUE)
