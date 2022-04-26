# calculate the proportion of GBIF observations that are in Switzerland.
# To correct sampling bias, use data at the ISEA10 grid
# (grid with area ~ 865km2, also used by IUCN red list)

library(tidyverse)
library(sf)
library(giscoR)
library(data.table)
library(rgbif)

select <- dplyr::select

# Read data ####################################################################
################################################################################.

# ISEA10 raster ----------------------------------------------------------------.
# accessible at https://www.iucnredlist.org/resources/spatialtoolsanddata
sf_isea10 <- st_read(dsn = "Data/IUCN_ISEA10/ISEA10_Global.shp")

# determine proportion of ISEA grid cells in Switzerland ########################
################################################################################.

sf_ch_border <- raster::getData(name = 'GADM', country = 'CHE', level = 0) %>%
  st_as_sf() %>%
  st_transform(crs = st_crs(sf_isea10))

# first cut, to improve computation time (not really necessary)
sf_ch_border_buff <- st_buffer(sf_ch_border, .5)

sf_use_s2(FALSE)
sf_isea10_sub <- sf_isea10 %>% 
  st_filter(sf_ch_border_buff)

sf_isea10_sub_ch <- st_intersection(sf_isea10_sub, sf_ch_border)
sf_isea10_sub_ch$area <- as.numeric(st_area(sf_isea10_sub_ch))

d_isea10_sub_ch <- sf_isea10_sub_ch %>% 
  as.data.frame() %>% 
  select(HexID, area) %>% 
  rename(ch_area = area)

sf_isea10$area <- as.numeric(st_area(sf_isea10))

# Butterflies ##################################################################
################################################################################.

# https://doi.org/10.15468/dl.t6ha3h
d_gbif_full <- fread("Data/GBIF/0071956-210914110416597/occurrence.csv")

# list of species to be analysed
specieslist <- fread("Other/specieslist.txt") %>% 
  filter(group == "butterflies") %>% 
  select(species) %>% 
  deframe()

d_prop_ch <- data.frame()

for (sp_i in specieslist[!specieslist %in% d_prop_ch$species]){
  print(paste0("Start ", sp_i, " at: ", Sys.time()))
  
  d_gbif <- d_gbif_full %>% 
    filter(species == sp_i)
    
  # transform to sf object
    sf_gbif <- d_gbif %>% 
      st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
               crs = "epsg:4121") %>% 
      st_transform(crs = st_crs(sf_isea10))
    
    # intersect with ISEA10 grid
    sf_use_s2(FALSE)
    sf_gbif <- st_intersection(sf_gbif, sf_isea10) 

    # calculate area proportion of recorded ISEA grid that are in CH
    out <- sf_gbif %>% 
      as.data.frame() %>% 
      select(HexID, area) %>% 
      distinct() %>% 
      left_join(d_isea10_sub_ch, by = "HexID") %>% 
      summarise(area_tot = sum(area),
                ch_area_tot = sum(ch_area, na.rm = T)) %>% 
      mutate(prop_ch = ch_area_tot / area_tot,
             species = sp_i)
    
    if (class(out) == "data.frame"){
      d_prop_ch <- d_prop_ch %>% 
        bind_rows(out)
    }
}

saveRDS(d_prop_ch, file = "Data/d_prop_ch_butterflies.rds")


# Grasshoppers #################################################################
################################################################################.

# https://doi.org/10.15468/dl.reemkv 
d_gbif_full <- fread("Data/GBIF/0083203-210914110416597/occurrence.csv")

# list of species to be analysed
specieslist <- fread("Other/specieslist.txt") %>% 
  filter(group == "grasshoppers") %>% 
  select(species) %>% 
  deframe()

d_prop_ch <- data.frame()

for (sp_i in specieslist[!specieslist %in% d_prop_ch$species]){
  print(paste0("Start ", sp_i, " at: ", Sys.time()))
  
  d_gbif <- d_gbif_full %>% 
    filter(species == sp_i)
  
  # transform to sf object
  sf_gbif <- d_gbif %>% 
    st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
             crs = "epsg:4121") %>% 
    st_transform(crs = st_crs(sf_isea10))
  
  # intersect with ISEA10 grid
  sf_use_s2(FALSE)
  sf_gbif <- st_intersection(sf_gbif, sf_isea10) 
  
  # calculate area proportion of recorded ISEA grid that are in CH
  out <- sf_gbif %>% 
    as.data.frame() %>% 
    select(HexID, area) %>% 
    distinct() %>% 
    left_join(d_isea10_sub_ch, by = "HexID") %>% 
    summarise(area_tot = sum(area),
              ch_area_tot = sum(ch_area, na.rm = T)) %>% 
    mutate(prop_ch = ch_area_tot / area_tot,
           species = sp_i)
  
  if (class(out) == "data.frame"){
    d_prop_ch <- d_prop_ch %>% 
      bind_rows(out)
  }
}

saveRDS(d_prop_ch, file = "Data/d_prop_ch_grasshoppers.rds")

# Dragonflies ##################################################################
################################################################################.

# https://doi.org/10.15468/dl.czbrmq
d_gbif_full <- fread("Data/GBIF/0079117-210914110416597/occurrence.csv")

# list of species to be analysed
specieslist <- fread("Other/specieslist.txt") %>% 
  filter(group == "dragonflies") %>% 
  select(species) %>% 
  deframe()

d_prop_ch <- data.frame()

for (sp_i in specieslist[!specieslist %in% d_prop_ch$species]){
  print(paste0("Start ", sp_i, " at: ", Sys.time()))
  
  d_gbif <- d_gbif_full %>% 
    filter(species == sp_i)
  
  # transform to sf object
  sf_gbif <- d_gbif %>% 
    st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
             crs = "epsg:4121") %>% 
    st_transform(crs = st_crs(sf_isea10))
  
  # intersect with ISEA10 grid
  sf_use_s2(FALSE)
  sf_gbif <- st_intersection(sf_gbif, sf_isea10) 
  
  
  # calculate area proportion of recorded ISEA grid that are in CH
  out <- sf_gbif %>% 
    as.data.frame() %>% 
    select(HexID, area) %>% 
    distinct() %>% 
    left_join(d_isea10_sub_ch, by = "HexID") %>% 
    summarise(area_tot = sum(area),
              ch_area_tot = sum(ch_area, na.rm = T)) %>% 
    mutate(prop_ch = ch_area_tot / area_tot,
           species = sp_i)
  
  if (class(out) == "data.frame"){
    d_prop_ch <- d_prop_ch %>% 
      bind_rows(out)
  }
}

saveRDS(d_prop_ch, file = "Data/d_prop_ch_dragonflies.rds")
