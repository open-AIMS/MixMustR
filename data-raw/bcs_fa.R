library(usethis)
library(dplyr)
bcs_fa <- "
source,Taxa,24:0,24:0 (SD),18:1w9,18:1w9 (SD),18:2w6,18:2w6 (SD),18:3w3,18:3w3 (SD),20:4w6,20:4w6 (SD),20:5w3,20:5w3 (SD),Study
Macroalgae,Sargassum,0.00,0.00,9.46,1.38,6.04,1.87,7.17,0.96,14.27,1.97,5.51,3.95,Khotimchenko 1991
Mangrove,A. marina,2.73,1.40,14.13,1.87,15.48,2.00,28.30,2.52,0.10,0.21,0.00,0.00,Meziane et al 2006
Saltmarsh,Salicornia bigelovii,3.30,3.30,9.10,7.64,18.35,4.45,24.85,32.60,0.00,0.00,0.00,0.00,Weete et al 1970
Seagrass,\"C. augustata, A antartica, H. uniervis\",0.00,0.00,4.00,1.00,27.00,5.00,33.00,12.00,1.00,0.50,1.00,1.00,Belicka et al 2012
Plankton,sediment trap 200-1430m,0.00,0.00,18.50,10.16,0.56,0.60,1.42,0.26,0.38,0.85,4.68,3.97,Burns et al 2003
Terrestrial grass,Phragmites australis,0.29,0.09,1.99,0.14,12.53,0.91,40.31,3.52,1.08,1.02,0.16,0.06,Wang et al 2014
" |>
  read.table(text = _, header = TRUE, check.names = FALSE, sep = ",") |>
  dplyr::arrange(source)
usethis::use_data(bcs_fa, overwrite = TRUE)
