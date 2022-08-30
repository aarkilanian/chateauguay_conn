##### Import and clean river data #####

# Read in GRHQ data

# CRHQ Data has been joined by mean with occ data (urban = ant, ag = cau, cfo, cgi, cin, cpi)

# Filter out imagined and looping rivers
###
# Select based on TYPECE
# Select based on ISOLE
# Select based on ENABLED
# Select based on PRIORITE
# Select based on FONCTION (kept 1 and 4, maybe 2?)

# Add single sink point and connect all single watershed sinks

##### Enforce dendritic geometry #####

# This section is done iteratively until the dendritic correction returns no errors

# Read in corrected GRHQ rivers
chat_riv <- read_sf(dsn = "../../../../../media/alex/ADATA SD700/data/Chateauguay/inputs.gpkg",
                    layer = "rivers_final")

# Correct out of bounds values
chat_riv$ag_mean[chat_riv$ag_mean > 100] <- 100

# Rescale land use % values to range of penalties
# For urban: 0.2 to 0.8
# For agric: 0.1 to 0.7
chat_riv$ag_pen <- 0.1 + chat_riv$ag_mean * 0.6 / 100
chat_riv$urb_pen <- 0.2 + chat_riv$ant_mean * 0.6 / 100
# Calculate total penalty
chat_riv$tot_pen <- chat_riv$ag_pen + chat_riv$urb_pen

# Convert penalties to weights
chat_riv$tot_w <- 1 - chat_riv$tot_pen

# Calculate alternative weights for invasives (0.6-1)
chat_riv$tot_w_inv <- 0.55 + chat_riv$tot_w * 0.4 / 0.7

# Repair geometry column
names(chat_riv)[names(chat_riv) == "geom"] <- "geometry"
st_geometry(chat_riv) <- "geometry"

##### Calculate restoration weights #####
set.seed(4523)

# 25%
chat_riv$tot_w_25 <- chat_riv$tot_w
rest_pot <- which(!is.na(chat_riv$tot_w_25))
restor <- sample(rest_pot, size = 0.25 * length(rest_pot))
chat_riv$tot_w_25[restor] <- 1
chat_riv$tot_w_inv_25 <- if_else(chat_riv$tot_w_25 == 1, 1,chat_riv$tot_w_inv)
# 50%
chat_riv$tot_w_50 <- chat_riv$tot_w
rest_pot <- which(!is.na(chat_riv$tot_w_50))
restor <- sample(rest_pot, size = 0.5 * length(rest_pot))
chat_riv$tot_w_50[restor] <- 1
chat_riv$tot_w_inv_50 <- if_else(chat_riv$tot_w_50 == 1, 1,chat_riv$tot_w_inv)
# 75%
chat_riv$tot_w_75 <- chat_riv$tot_w
rest_pot <- which(!is.na(chat_riv$tot_w_75))
restor <- sample(rest_pot, size = 0.75 * length(rest_pot))
chat_riv$tot_w_75[restor] <- 1
chat_riv$tot_w_inv_75 <- if_else(chat_riv$tot_w_75 == 1, 1,chat_riv$tot_w_inv)

##### Import points #####

# Import dams
chat_dams <- read_sf(dsn = "../../../../../media/alex/ADATA SD700/data/Chateauguay/inputs.gpkg",
                     layer = "dams")
chat_dams$bartype <- "dam"
# Import culverts
chat_culv <- read_sf(dsn = "../../../../../media/alex/ADATA SD700/data/Chateauguay/inputs.gpkg",
                     layer = "culverts")
chat_culv$bartype <- "culvert"
# Join barriers
chat_bars <- chat_dams %>%
  select(bartype, geom) %>%
  bind_rows(chat_culv %>% select(bartype, geom)) %>%
  st_zm()
# Remove any overlapping geometries
chat_bars <- distinct(chat_bars)
# Repair geometry column
names(chat_bars)[names(chat_bars) == "geom"] <- "geometry"
st_geometry(chat_bars) <- "geometry"

# Import outlet
chat_out <- read_sf(dsn = "../../../../../media/alex/ADATA SD700/data/Chateauguay/inputs.gpkg",
                    layer = "outlet")
# Repair geometry column
names(chat_out)[names(chat_out) == "geom"] <- "geometry"
st_geometry(chat_out) <- "geometry"

##### Associate passability values for barriers #####

chat_bars <- chat_bars %>%
  mutate(pass_simple = case_when(
    bartype == "culvert" ~ 0.8,
    bartype == "dam" ~ 0.2
  ))

##### Prepare DCI inputs & river network #####

# Import barriers
chat_bars <- import_points(chat_bars, type = "barriers")

# Import outlet
chat_outlet <- import_points(chat_out, type = "outlet")

# Import rivers
chat_riv <- import_rivers(chat_riv)

# Prepare network
chat_net <- river_net(rivers = chat_riv, barriers = chat_bars,
                      outlet = chat_outlet, check = FALSE)



##### Calculate some connectivity #####

# Potamodromous DCI
#
# 0%
full_pot_0 <- calculate_dci(net = chat_net,
                          form = "potamodromous",
                          pass = "pass_simple",
                          weight = 'tot_w',
                          n.cores = 12)
full_pot_25 <- calculate_dci(net = chat_net,
                            form = "potamodromous",
                            pass = "pass_simple",
                            weight = 'tot_w_25',
                            n.cores = 12)
full_pot_50 <- calculate_dci(net = chat_net,
                            form = "potamodromous",
                            pass = "pass_simple",
                            weight = 'tot_w_50',
                            n.cores = 12)
full_pot_75 <- calculate_dci(net = chat_net,
                            form = "potamodromous",
                            pass = "pass_simple",
                            weight = 'tot_w_75',
                            n.cores = 12)

# Diadromous DCI
#
# 2 seconds
full_dia_0 <- calculate_dci(net = chat_net,
                          form = "diadromous",
                          pass = "pass_simple",
                          weight = "tot_w_inv",
                          n.cores = 12)
full_dia_25 <- calculate_dci(net = chat_net,
                          form = "diadromous",
                          pass = "pass_simple",
                          weight = "tot_w_inv_25",
                          n.cores = 12)
full_dia_50 <- calculate_dci(net = chat_net,
                          form = "diadromous",
                          pass = "pass_simple",
                          weight = "tot_w_inv_50",
                          n.cores = 12)
full_dia_75 <- calculate_dci(net = chat_net,
                          form = "diadromous",
                          pass = "pass_simple",
                          weight = "tot_w_inv_75",
                          n.cores = 12)

# Retrieved mean dispersal across a large meta-study
# https://royalsocietypublishing.org/doi/10.1098/rspb.2017.2214
# Assume average 500m/day dispersal
# Over a 4 week yearly spawning period that is 14,000m
#
# 4 minutes
thresh_pot_0 <- calculate_dci(net = chat_net,
                          form = "potamodromous",
                          pass = "pass_simple",
                          weight = 'tot_w',
                          threshold = 14000,
                          n.cores = 12)
thresh_pot_25 <- calculate_dci(net = chat_net,
                            form = "potamodromous",
                            pass = "pass_simple",
                            weight = 'tot_w_25',
                            threshold = 14000,
                            n.cores = 12)
thresh_pot_50 <- calculate_dci(net = chat_net,
                            form = "potamodromous",
                            pass = "pass_simple",
                            weight = 'tot_w_50',
                            threshold = 14000,
                            n.cores = 12)
thresh_pot_75 <- calculate_dci(net = chat_net,
                            form = "potamodromous",
                            pass = "pass_simple",
                            weight = 'tot_w_75',
                            threshold = 14000,
                            n.cores = 12)

# Calculate mean rankings
mean_rank <- (pot_0_rank+ pot_25_rank+ pot_50_rank+ pot_75_rank+
             dia_0_rank+ dia_25_rank+ dia_50_rank+ dia_75_rank+
             thresh_0_rank+ thresh_25_rank+ thresh_50_rank+ thresh_75_rank)/12

##### Export results #####

rivers <- export_dci(chat_net, full_pot_0, "rivers") %>%
  select(-DCI, -DCI_rel) %>%
  dplyr::left_join(full_pot_0, by = c("member.label" = "segment")) %>%
  dplyr::left_join(full_pot_25, by = c("member.label" = "segment"), suffix = c(".fullp00", ".fullp25")) %>%
  dplyr::left_join(full_pot_50, by = c("member.label" = "segment")) %>%
  dplyr::left_join(full_pot_75, by = c("member.label" = "segment"), suffix = c(".fullp50", ".fullp75")) %>%
  dplyr::left_join(full_dia_0, by = c("member.label" = "segment")) %>%
  dplyr::left_join(full_dia_25, by = c("member.label" = "segment"), suffix = c(".fulld00", ".fulld25")) %>%
  dplyr::left_join(full_dia_50, by = c("member.label" = "segment")) %>%
  dplyr::left_join(full_dia_75, by = c("member.label" = "segment"), suffix = c(".fulld50", ".fulld75")) %>%
  dplyr::left_join(thresh_pot_0, by = c("member.label" = "segment")) %>%
  dplyr::left_join(thresh_pot_25, by = c("member.label" = "segment"), suffix = c(".thresh00", ".thresh25")) %>%
  dplyr::left_join(thresh_pot_50, by = c("member.label" = "segment")) %>%
  dplyr::left_join(thresh_pot_75, by = c("member.label" = "segment"), suffix = c(".thresh50", ".thresh75")) %>%
  select(O_STRAHLER, riv_length, ant_mean:DCI_rel.thresh75) %>%
  select(-ag_pen, -urb_pen, -tot_pen, -rivID, -geometry.y, -fid.y)

barriers <- sf::st_as_sf(activate(chat_net, nodes)) %>%
  dplyr::filter(.data$type %in% c("barrier", "outlet")) %>%
  dplyr::left_join(full_pot_0, by = c("member.label" = "segment")) %>%
  dplyr::left_join(full_pot_25, by = c("member.label" = "segment"), suffix = c(".fullp00", ".fullp25")) %>%
  dplyr::left_join(full_pot_50, by = c("member.label" = "segment")) %>%
  dplyr::left_join(full_pot_75, by = c("member.label" = "segment"), suffix = c(".fullp50", ".fullp75")) %>%
  dplyr::left_join(full_dia_0, by = c("member.label" = "segment")) %>%
  dplyr::left_join(full_dia_25, by = c("member.label" = "segment"), suffix = c(".fulld00", ".fulld25")) %>%
  dplyr::left_join(full_dia_50, by = c("member.label" = "segment")) %>%
  dplyr::left_join(full_dia_75, by = c("member.label" = "segment"), suffix = c(".fulld50", ".fulld75")) %>%
  dplyr::left_join(thresh_pot_0, by = c("member.label" = "segment")) %>%
  dplyr::left_join(thresh_pot_25, by = c("member.label" = "segment"), suffix = c(".thresh00", ".thresh25")) %>%
  dplyr::left_join(thresh_pot_50, by = c("member.label" = "segment")) %>%
  dplyr::left_join(thresh_pot_75, by = c("member.label" = "segment"), suffix = c(".thresh50", ".thresh75")) %>%
  select(-type, -fid, -node.label)

# Compute ranks
barriers <- barriers %>%
  # Compute ranks
  mutate(rankp0 = rank(DCI.fullp00)) %>%
  mutate(rankp25 = rank(DCI.fullp25)) %>%
  mutate(rankp50 = rank(DCI.fullp50)) %>%
  mutate(rankp75 = rank(DCI.fullp75)) %>%
  mutate(rankd0 = rank(DCI.fulld00)) %>%
  mutate(rankd25 = rank(DCI.fulld25)) %>%
  mutate(rankd50 = rank(DCI.fulld50)) %>%
  mutate(rankd75 = rank(DCI.fulld75)) %>%
  mutate(rankth0 = rank(DCI.thresh00)) %>%
  mutate(rankth25 = rank(DCI.thresh25)) %>%
  mutate(rankth50 = rank(DCI.thresh50)) %>%
  mutate(rankth75 = rank(DCI.thresh75)) %>%
  # Compute percentage rank
  mutate(prankp0 = 100 - ((rankp0 - min(rankp0)) / (max(rankp0) - min(rankp0)) * 100)) %>%
  mutate(prankp25 = 100 - ((rankp25 - min(rankp25)) / (max(rankp25) - min(rankp25)) * 100)) %>%
  mutate(prankp50 = 100 - ((rankp50 - min(DCI.fullp50, na.rm = T)) / (max(DCI.fullp50, na.rm = T) - min(DCI.fullp50, na.rm = T)) * 100)) %>%
  mutate(prankp75 = 100 - ((rankp75 - min(rankp50)) / (max(rankp50) - min(rankp50)) * 100)) %>%
  mutate(prankd0 = (rankd0 - min(rankd0)) / (max(rankd0) - min(rankd0)) * 100) %>%
  mutate(prankd25 = (rankd25 - min(rankd25)) / (max(rankd25) - min(rankd25)) * 100) %>%
  mutate(prankd50 = (rankd50 - min(rankd50)) / (max(rankd50) - min(rankd50)) * 100) %>%
  mutate(prankd75 = (rankd75 - min(rankd75)) / (max(rankd75) - min(rankd75)) * 100) %>%
  mutate(prankth0 = 100 - ((rankth0 - min(rankth0)) / (max(rankth0) - min(rankth0)) * 100)) %>%
  mutate(prankth25 = 100 - ((rankth25 - min(rankth25)) / (max(rankth25) - min(rankth25)) * 100)) %>%
  mutate(prankth50 = 100 - ((rankth50 - min(rankth50)) / (max(rankth50) - min(rankth50)) * 100)) %>%
  mutate(prankth75 = 100 - ((rankth75 - min(rankth75)) / (max(rankth75) - min(rankth75)) * 100)) %>%
  # Calculate combined score
  mutate(comb_score = (0.0625*prankp0 + 0.0625*prankp25 + 0.0625*prankp50 + 0.0625*prankp75 +
                              0.125*prankd0 + 0.125*prankd25 + 0.125*prankd50 + 0.125*prankd75 +
                              0.0625*prankth0 + 0.0625*prankth25 + 0.0625*prankth50 + 0.0625*prankth75)) %>%
  # Compute final rank
  mutate(final_rank = rank(comb_score))

# Combine ranks to rivers
rivers <- rivers %>%
  left_join(barriers %>% select(member.label, final_rank) %>% as.data.frame())
