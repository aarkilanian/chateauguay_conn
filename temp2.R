require(sf)
require(sfnetworks)
require(dci)

riv <- read_sf("../chateauguay_conn/data/inputs_new.gpkg", layer = "riv_bigcomplex3") %>%
  st_zm()

riv_in <- import_rivers(riv)

a <- enforce_dendritic(riv_in)

b <- enforce_dendritic(riv_in, correct = TRUE)
