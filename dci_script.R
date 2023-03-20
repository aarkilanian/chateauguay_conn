##### Import packages #####

# Installation instructions (if required)

# install.packages("dplyr")
# install.packages("sf")
# install.packages("sfnetworks")
# install.packages("devtools")
# devtools::install_github("aarkilanian/dci")

# Import packages
require(dplyr)
require(sf)
require(dci)

##### Import data #####

# Import rivers from MELCC's GRHQ (RH_L layer)
# https://www.donneesquebec.ca/recherche/dataset/grhq
rivers <- read_sf("")

# Import dams from MELCC
# https://www.cehq.gouv.qc.ca/barrages/default.asp#version-telechargeable
dams <- read_sf("")

# Import culverts from MTQ
# https://www.donneesquebec.ca/recherche/dataset/structure
culverts <- read_sf("")

# Import waterfalls from MELCC' GRHQ (C_hyd_p layer)
falls <- read_sf("")

##### Prepare rivers #####

# Clean river data
rivers <- rivers %>%
  # Discard Z dimension (sf cannot work with Z or M data)
  st_zm() %>%
  # Filter rivers with ENABLED set to 1
  filter(ENABLED == 1) %>%

# Clip rivers to Quebec boundary and add rivers to connect outlets
# Recommended to be done in GIS software
# Export rivers
st_write(rivers, "")
# Import clipped rivers with all outlets flowing to a single final outlet
rivers <- read_sf("")

# Read in clipped rivers
rivers <- read_sf("")

# Initial import to dci to remove small fragments
rivers <- import_rivers(rivers)

# Detect topological errors
rivers_err <- enforce_dendritic(rivers, correct = FALSE)
# Export river errors
st_write(rivers_err, "")
# Import final corrected rivers
rivers <- read_sf("")
# Should be corrected then rerun until no errors remain

##### Prepare barriers & outlet #####

# Import outlet
outlet <- read_sf("")

# Prepare dams
dams <- dams %>%
# Match river data projection
  st_transform(st_crs(rivers)) %>%
  # Rename ID
  rename("original_id" = "ide_strct") %>%
  mutate(original_id = as.character(original_id)) %>%
  # Assign culvert type
  mutate(type = "culvert") %>%
  # Assign passability
  mutate(pass = 0.8) %>%
  # Select relevant columns
  select(original_id, type, pass, geom)

# Prepare culverts
culverts <- culverts %>%
  # Match river data projection
  st_transform(st_crs(rivers)) %>%
  # Rename ID
  rename("original_id" = "Numéro.ba") %>%
  mutate(original_id = as.character(original_id)) %>%
  # Assign culvert type
  mutate(type = "dam") %>%
  # Assign passability
  mutate(pass = 0.2) %>%
  # Select relevant columns
  select(original_id, type, pass, geom)

# Prepare waterfalls
falls <- falls %>%
  # Match river data projection
  st_transform(st_crs(rivers)) %>%
  # Rename ID
  rename("original_id" = "ID_RHP") %>%
  mutate(original_id = as.character(original_id)) %>%
  # Assign culvert type
  mutate(type = "falls") %>%
  # Assign passability
  mutate(pass = 0) %>%
  # Select relevant columns
  select(original_id, type, pass, geom)

# Combine barriers
all_bar <- bind_rows(culverts, dams, falls)

# Snap barriers to rivers
# Should be done in GIS software
# Export barriers
st_write(all_bar, "")
# Load snapped barriers
all_bar <- ""

##### Compute weights #####

# In GIS software summarize environmental data over rivers
# For example using the CRHQ
# https://www.donneesquebec.ca/recherche/dataset/crhq
# Export rivers
st_write(rivers, "")
# Import joined rivers
rivers <- read_sf("")

# Convert summarized variables to penalties
# An example of this using anthropogenic and agricultural land use:

# Convert land use to weighting
rivers <- rivers %>%
  # Compute agricultural penalty
  mutate(
    penalty_agr = case_when(
      agr_mean < 30 ~ 0,
      agr_mean >= 30 & agr_mean < 50 ~ 0.1,
      agr_mean >= 50 & agr_mean < 70 ~ 0.2,
      agr_mean >= 70 & agr_mean < 90 ~ 0.3,
      agr_mean >= 90 ~ 0.4
    )
  ) %>%
  # Compute urban penalty
  mutate(
    penalty_urb = case_when(
      ant_mean < 30 ~ 0,
      ant_mean >= 30 & ant_mean < 50 ~ 0.2,
      ant_mean >= 50 & ant_mean < 70 ~ 0.3,
      ant_mean >= 70 & ant_mean < 90 ~ 0.4,
      ant_mean >= 90 ~ 0.5
    )
  ) %>%
  # Combine penalties to get weighting
  mutate(weight = 1 - (penalty_urb + penalty_agr)) %>%
  # Set minimum weight to 0.05
  mutate(weight = pmax(weight, 0.05))

# Set imagined outlet line weights to 0
# These should be identified with a unique variable, in this case it's "included"
rivers <- rivers %>%
  mutate(weight = if_else(included == "0", 0, weight))

##### Compute restoration scenarios #####

# Name the recently calculated weight as the 0% restoration weight
rivers <- rivers %>%
  rename(weight_0 = weight)

# Determine restoration targets for each scenario
rest_25 <- sample(1:nrow(rivers), size = 0.25*nrow(rivers))
rest_50 <- sample(1:nrow(rivers), size = 0.50*nrow(rivers))
rest_75 <- sample(1:nrow(rivers), size = 0.75*nrow(rivers))

# Compute restoration scenarios
rivers <- rivers %>%
  # 25% restoration
  mutate(weight_25 = if_else(row_number() %in% rest_25, 1, weight_0)) %>%
  # 50% restoration
  mutate(weight_50 = if_else(row_number() %in% rest_50, 1, weight_0)) %>%
  # 75% restoration
  mutate(weight_75 = if_else(row_number() %in% rest_75, 1, weight_0))

# Compute invasive habitat quality based on original quality measure
rivers <- rivers %>%
  # 0% restoration
  mutate(weight_inv_0 = pmin(weight_0 + 0.1, 1)) %>%
  # 25% restoration
  mutate(weight_inv_25 = pmin(weight_25 + 0.1, 1)) %>%
  # 50% restoration
  mutate(weight_inv_50 = pmin(weight_50 + 0.1, 1)) %>%
  # 75% restoration
  mutate(weight_inv_75 = pmin(weight_75 + 0.1, 1))

##### Calculate connectivity #####

# Import rivers
rivers_fin <- import_rivers(rivers)

# Import barriers
bars_fin <- import_points(all_bar, "barriers")

# Import outlet
outlet_fin <- import_points(outlet, "outlet")

# Generate river_net object
net <- river_net(rivers = chat_riv,
                 barriers = chat_bars,
                 outlet = chat_out,
                 check = FALSE)

# Calculating long-term potamodromous DCI (no distance threshold)
# 0% restoration
full_pot_0 <- calculate_dci(net = net,
                            form = "potamodromous",
                            pass = "pass",
                            weight = 'weight_0',
                            n.cores = 12)
# 25% restoration
full_pot_25 <- calculate_dci(net = net,
                             form = "potamodromous",
                             pass = "pass",
                             weight = 'weight_25',
                             n.cores = 12)
# 50% restoration
full_pot_50 <- calculate_dci(net = net,
                             form = "potamodromous",
                             pass = "pass",
                             weight = 'weight_50',
                             n.cores = 12)
# 75% restoration
full_pot_75 <- calculate_dci(net = net,
                             form = "potamodromous",
                             pass = "pass",
                             weight = 'weight_75',
                             n.cores = 12)

# Calculating short-term potamodromous DCI
# 0% restoration
thresh_pot_0 <- calculate_dci(net = net,
                              form = "potamodromous",
                              pass = "pass",
                              weight = 'weight_0',
                              threshold = 14000,
                              n.cores = 12)
# 25% restoration
thresh_pot_25 <- calculate_dci(net = net,
                               form = "potamodromous",
                               pass = "pass",
                               weight = 'weight_25',
                               threshold = 14000,
                               n.cores = 12)
# 50% restoration
thresh_pot_50 <- calculate_dci(net = net,
                               form = "potamodromous",
                               pass = "pass",
                               weight = 'weight_50',
                               threshold = 14000,
                               n.cores = 12)
# 75% restoration
thresh_pot_75 <- calculate_dci(net = net,
                               form = "potamodromous",
                               pass = "pass",
                               weight = 'weight_75',
                               threshold = 14000,
                               n.cores = 12)

# Calculating diadromous invasive DCI
# 0% restoration
full_dia_0 <- calculate_dci(net = net,
                            form = "diadromous",
                            pass = "pass",
                            weight = "weight_inv_0",
                            n.cores = 12)
# 25% restoration
full_dia_25 <- calculate_dci(net = net,
                             form = "diadromous",
                             pass = "pass",
                             weight = "weight_inv_25",
                             n.cores = 12)
# 50% restoration
full_dia_50 <- calculate_dci(net = net,
                             form = "diadromous",
                             pass = "pass",
                             weight = "weight_inv_50",
                             n.cores = 12)
# 75% restoration
full_dia_75 <- calculate_dci(net = net,
                             form = "diadromous",
                             pass = "pass",
                             weight = "weight_inv_75",
                             n.cores = 12)

# Export results onto rivers
rivers <- export_dci(net, full_pot_0, "rivers") %>%
  select(-DCI, -DCI_rel) %>%
  dplyr::left_join(full_pot_0, by = c("member.label" = "segment")) %>%
  dplyr::left_join(full_pot_25, by = c("member.label" = "segment"), 
                   suffix = c("_pot_0", "_pot_25")) %>%
  dplyr::left_join(full_pot_50, by = c("member.label" = "segment")) %>%
  dplyr::left_join(full_pot_75, by = c("member.label" = "segment"), 
                   suffix = c("_pot_50", "_pot_75")) %>%
  dplyr::left_join(full_dia_0, by = c("member.label" = "segment")) %>%
  dplyr::left_join(full_dia_25, by = c("member.label" = "segment"), 
                   suffix = c("_dia_0", "_dia_25")) %>%
  dplyr::left_join(full_dia_50, by = c("member.label" = "segment")) %>%
  dplyr::left_join(full_dia_75, by = c("member.label" = "segment"), 
                   suffix = c("_dia_50", "_dia_75")) %>%
  dplyr::left_join(thresh_pot_0, by = c("member.label" = "segment")) %>%
  dplyr::left_join(thresh_pot_25, by = c("member.label" = "segment"), 
                   suffix = c("_thr_0", "_thr_25")) %>%
  dplyr::left_join(thresh_pot_50, by = c("member.label" = "segment")) %>%
  dplyr::left_join(thresh_pot_75, by = c("member.label" = "segment"), 
                   suffix = c("_thr_50", "_thr_75")) %>%
  select(O_STRAHLER, riv_length, ant_mean:DCI_rel.thresh75) %>%
  select(-ag_pen, -urb_pen, -tot_pen, -rivID, -geometry.y, -fid.y)

# Export results onto barriers
barriers <- sf::st_as_sf(activate(net, nodes)) %>%
  dplyr::filter(.data$type %in% c("barrier", "outlet")) %>%
  dplyr::left_join(full_pot_0, by = c("member.label" = "segment")) %>%
  dplyr::left_join(full_pot_25, by = c("member.label" = "segment"), 
                   suffix = c("_pot_0", "_pot_25")) %>%
  dplyr::left_join(full_pot_50, by = c("member.label" = "segment")) %>%
  dplyr::left_join(full_pot_75, by = c("member.label" = "segment"), 
                   suffix = c("_pot_50", "_pot_75")) %>%
  dplyr::left_join(full_dia_0, by = c("member.label" = "segment")) %>%
  dplyr::left_join(full_dia_25, by = c("member.label" = "segment"), 
                   suffix = c("_dia_0", "_dia_25")) %>%
  dplyr::left_join(full_dia_50, by = c("member.label" = "segment")) %>%
  dplyr::left_join(full_dia_75, by = c("member.label" = "segment"), 
                   suffix = c("_dia_50", "_dia_75")) %>%
  dplyr::left_join(thresh_pot_0, by = c("member.label" = "segment")) %>%
  dplyr::left_join(thresh_pot_25, by = c("member.label" = "segment"), 
                   suffix = c("_thr_0", "_thr_25")) %>%
  dplyr::left_join(thresh_pot_50, by = c("member.label" = "segment")) %>%
  dplyr::left_join(thresh_pot_75, by = c("member.label" = "segment"),
                   suffix = c("_thr_50", "_thr_75")) %>%
  select(-type, -fid, -node.label)

##### Compute ranks from connectivity results #####

# Compute ranks on barriers
barriers <- barriers %>%
  mutate(rankp0 = rank(-DCI_pot_0)) %>%
  mutate(rankp25 = rank(-DCI_pot_25)) %>%
  mutate(rankp50 = rank(-DCI_pot_50)) %>%
  mutate(rankp75 = rank(-DCI_pot_75)) %>%
  mutate(rankd0 = rank(DCI_dia_0)) %>%
  mutate(rankd25 = rank(DCI_dia_25)) %>%
  mutate(rankd50 = rank(DCI_dia_50)) %>%
  mutate(rankd75 = rank(DCI_dia_75)) %>%
  mutate(rankth0 = rank(-DCI_thr_0)) %>%
  mutate(rankth25 = rank(-DCI_thr_25)) %>%
  mutate(rankth50 = rank(-DCI_thr_50)) %>%
  mutate(rankth75 = rank(-DCI_thr_75)) %>%
  # Compute rank product
  mutate(rank_product = rankp0 * rankp25 * rankp50 * rankp75 *
           rankd0 * rankd25 * rankd50 * rankd75 *
           rankth0 * rankth25 * rankth50 * rankth75 / 12) %>%
  # Compute final rank
  mutate(final_rank = rank(rank_product))
