library(tidyverse)
library(raster)
library(sf)
library(stars)
library(giscoR)
library(data.table)
library(rgbif)

select <- dplyr::select

# European raster CGRS----------------------------------------------------------.

# CGRS grid available from https://www.eea.europa.eu/data-and-maps/data/common-european-chorological-grid-reference-system-cgrs
sf_cgrs <- st_read(dsn = "Data/cgrs_grid/cgrs_grid.shp")
sf_cgrs <- sf_cgrs[-11267, ] # something wrong with this one square

# European countries (NUTS) ----------------------------------------------------.

sf_europe <- gisco_get_countries(region = "Europe")

# to cut off overseas areas
cut_pol <-  st_sfc(st_polygon(list(cbind(c(-5, 40, 45, -40, -5),
                                         c(25, 30, 75, 65, 25)))),
                   crs = 4326) %>% 
  st_transform(crs = st_crs(sf_europe)) # harmonize coordinate system

# cut off overseas areas and other countries
sf_europe <- st_intersection(sf_europe, cut_pol)
sf_europe <- sf_europe %>% 
  # Exclusions following Schweiger et al. 2014
  filter(!NAME_ENGL %in% c("Russian Federation", "Belarus", "Ukraine", "Iceland",
                           "Svalbard and Jan Mayen", "Faroes")) %>% 
  st_transform(crs = st_crs(sf_cgrs))

# read and process Wordclim data ###############################################
################################################################################.

# data available from https://www.worldclim.org/data/worldclim21.html

# Europen subset of CGRS grid
sf_cgrs_europe <- sf_cgrs %>% 
  st_filter(sf_europe)

# mean number of days in the included months
d_days <- data.frame()
for (mn_i in 1:12){
  out <- c()
  for (yr_i in 1970:2020){
    out <- c(out, lubridate::days_in_month(as.Date(paste(yr_i, mn_i, 1, sep = "-"))))
  }
  d_days <- data.frame(month = mn_i, days = mean(out)) %>% 
    bind_rows(d_days, .)
}


# read through the twelve monthly files
dir <- "Data/Wordclim/wc2.1_2.5m_tavg"
cut <- st_buffer(st_union(sf_cgrs_europe), 100000)

r_wordclim <- raster(list.files(dir, full.names = T)[1]) # first month
r_wordclim <- crop(r_wordclim, extent(st_bbox(cut)))
r_wordclim <- r_wordclim * d_days$days[1]

# add other months
for (i in 2:12){
  file <- list.files(dir, full.names = T)[i]
  
  r_target <- raster(file)
  r_target <- crop(r_target, extent(st_bbox(cut)))
  
  r_wordclim <- r_wordclim + r_target *  d_days$days[i]
}

r_wordclim <- r_wordclim / sum(d_days$days) # mean temperature

# convert to sf
sf_wordclim <- st_as_sf(st_as_stars(r_wordclim))
sf_wordclim <- st_transform(sf_wordclim, crs = st_crs(sf_cgrs_europe))

names(sf_wordclim)[1] <- "Tniche"

# convert to sf
sf_wordclim <- st_as_sf(st_as_stars(r_wordclim))
sf_wordclim <- st_transform(sf_wordclim, crs = st_crs(sf_cgrs_europe))

names(sf_wordclim)[1] <- "Tniche"

# Average T over CGRS grid #####################################################
################################################################################.

sf_wordclim_cgrs <- st_join(sf_wordclim, sf_cgrs_europe, 
                            join = st_within) # choose only those completely within CGRS grid

sf_wordclim_cgrs <- sf_wordclim_cgrs %>% 
  filter(!is.na(CGRSNAME)) # points outside grid / between two grid cells

# mean across grid
sf_wordclim_cgrs <- sf_wordclim_cgrs %>% 
  group_by(CGRSNAME) %>% 
  summarise(Tniche = mean(Tniche),
            n = n())

# add to CGRS data
d_cgrs_europe <- sf_cgrs_europe %>% 
  as.data.frame() %>% 
  left_join(sf_wordclim_cgrs, by = "CGRSNAME")


sf_cgrs_europe$Tniche <- d_cgrs_europe$Tniche
sf_cgrs_europe <- sf_cgrs_europe %>% 
  filter(!is.na(Tniche))

################################################################################.
# BUTTERFLIES ---------- #######################################################
################################################################################.

# Read data ####################################################################
################################################################################.

# GBIF records -----------------------------------------------------------------.
# https://doi.org/10.15468/dl.t6ha3h

d_gbif_full <- fread("Data/GBIF/0071956-210914110416597/occurrence.csv")

# Import GBIF species data and determine category ##############################
################################################################################.

# list of species to be analysed
specieslist <- fread("Other/specieslist.txt") %>% 
  filter(group == "butterflies") %>% 
  select(species) %>% 
  deframe()

sf_europe_buff <- st_buffer(sf_europe, 20000)

sp_i <- specieslist[1]
d_temperature_niche <- data.frame()
for (sp_i in specieslist[!specieslist %in% d_temperature_niche$species]){
  print(paste0("Start ", sp_i, " at: ", Sys.time()))
  
  d_gbif <- d_gbif_full %>% 
    filter(species == sp_i)
  
  # transform to sf object
  sf_gbif <- d_gbif %>% 
    st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
             crs = "epsg:4121") %>% 
    st_transform(crs = st_crs(sf_cgrs_europe))
  
  # exclude non-european observations
  sel <- st_intersects(sf_gbif, sf_europe_buff, sparse = T)
  sel <- which(unlist(lapply(sel, length)) > 0)
  sf_gbif <- sf_gbif[sel, ]
  
  sf_gbif <- st_join(sf_gbif, sf_cgrs_europe)
  
  out <- sf_gbif %>% 
    as.data.frame() %>%
    select(CGRSNAME, Tniche) %>% 
    distinct() %>% 
    summarise(Tniche_mean = mean(Tniche, na.rm = T),
              Tniche_sd = sd(Tniche, na.rm = T),
              Tniche_q05 = quantile(Tniche, .05, na.rm = T),
              Tniche_q95 = quantile(Tniche, .95, na.rm = T))
  
  out$species <- sp_i
  
  if (class(out) == "data.frame"){
    d_temperature_niche <- d_temperature_niche %>% 
      bind_rows(out)
  }
}


saveRDS(d_temperature_niche, file = "Data/d_temperature_niche_butterflies.rds")

################################################################################.
# GRASSHOPPERS -------------- ##################################################
################################################################################.

# Read data ####################################################################
################################################################################.

# GBIF records -----------------------------------------------------------------.
# https://doi.org/10.15468/dl.reemkv

d_gbif_full <- fread("Data/GBIF/0083203-210914110416597/occurrence.csv")

# Import GBIF species data and determine category ##############################
################################################################################.

# list of species to be analysed
specieslist <- fread("Other/specieslist.txt") %>% 
  filter(group == "grasshoppers") %>% 
  select(species) %>% 
  deframe()

sf_europe_buff <- st_buffer(sf_europe, 20000)

sp_i <- specieslist[1]
d_temperature_niche <- data.frame()
for (sp_i in specieslist[!specieslist %in% d_temperature_niche$species]){
  print(paste0("Start ", sp_i, " at: ", Sys.time()))
  
  
  d_gbif <- d_gbif_full %>% 
    filter(species == sp_i)
  
  # transform to sf object
  sf_gbif <- d_gbif %>% 
    st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
             crs = "epsg:4121") %>% 
    st_transform(crs = st_crs(sf_cgrs_europe))
  
  # exclude non-european observations
  sel <- st_intersects(sf_gbif, sf_europe_buff, sparse = T)
  sel <- which(unlist(lapply(sel, length)) > 0)
  sf_gbif <- sf_gbif[sel, ]
  
  sf_gbif <- st_join(sf_gbif, sf_cgrs_europe)
  
  out <- sf_gbif %>% 
    as.data.frame() %>%
    select(CGRSNAME, Tniche) %>% 
    distinct() %>% 
    summarise(Tniche_mean = mean(Tniche, na.rm = T),
              Tniche_sd = sd(Tniche, na.rm = T),
              Tniche_q05 = quantile(Tniche, .05, na.rm = T),
              Tniche_q95 = quantile(Tniche, .95, na.rm = T))
  
  out$species <- sp_i
  
  if (class(out) == "data.frame"){
    d_temperature_niche <- d_temperature_niche %>% 
      bind_rows(out)
  }
  
  
  print(paste0(sp_i, " done at: ", Sys.time()))
}

saveRDS(d_temperature_niche, file = "Data/d_temperature_niche_grasshoppers.rds")

################################################################################.
# DRAGONFLIES -------------- ###################################################
################################################################################.

# Read data ####################################################################
################################################################################.

# GBIF records -----------------------------------------------------------------.
# https://doi.org/10.15468/dl.czbrmq

d_gbif_full <- fread("Data/GBIF/0079117-210914110416597/occurrence.csv")

# Import GBIF species data and determine category ##############################
################################################################################.

# list of species to be analysed
specieslist <- fread("Other/specieslist.txt") %>% 
  filter(group == "dragonflies") %>% 
  select(species) %>% 
  deframe()

sf_europe_buff <- st_buffer(sf_europe, 20000)

sp_i <- specieslist[1]
d_temperature_niche <- data.frame()
for (sp_i in specieslist[!specieslist %in% d_temperature_niche$species]){
  print(paste0("Start ", sp_i, " at: ", Sys.time()))
  
  
  d_gbif <- d_gbif_full %>% 
    filter(species == sp_i)
  
  # transform to sf object
  sf_gbif <- d_gbif %>% 
    st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
             crs = "epsg:4121") %>% 
    st_transform(crs = st_crs(sf_cgrs_europe))
  
  # exclude non-european observations
  sel <- st_intersects(sf_gbif, sf_europe_buff, sparse = T)
  sel <- which(unlist(lapply(sel, length)) > 0)
  sf_gbif <- sf_gbif[sel, ]
  
  sf_gbif <- st_join(sf_gbif, sf_cgrs_europe)
  
  out <- sf_gbif %>% 
    as.data.frame() %>%
    select(CGRSNAME, Tniche) %>% 
    distinct() %>% 
    summarise(Tniche_mean = mean(Tniche, na.rm = T),
              Tniche_sd = sd(Tniche, na.rm = T),
              Tniche_q05 = quantile(Tniche, .05, na.rm = T),
              Tniche_q95 = quantile(Tniche, .95, na.rm = T))
  
  out$species <- sp_i
  
  if (class(out) == "data.frame"){
    d_temperature_niche <- d_temperature_niche %>% 
      bind_rows(out)
  }
}

saveRDS(d_temperature_niche, file = "Data/d_temperature_niche_dragonflies.rds")