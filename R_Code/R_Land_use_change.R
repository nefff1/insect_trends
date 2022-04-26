library(tidyverse); theme_set(theme_classic())
library(data.table)
library(sf)
library(sfheaders)
library(rstanarm)

rename <- dplyr::rename
options(mc.cores = 4)

################################################################################.
# READ AND EDIT DATA ----------- ###############################################
################################################################################.

# Read national agricultural statistics data -----------------------------------.

# Total agricultural area per year and municipality
d_AgrArea <- fread("Data/d_AgrArea.csv")

# Grassland area per year and municipality
d_grasslands <- fread("Data/d_grasslands.csv")

# Number of livestock units per year and municipality
d_LSU <- fread("Data/d_LSU.csv")

# Area of different crops per year and municipality
d_cultures <- fread("Data/d_cultures.csv")

# Average insecticide application rates ----------------------------------------.
# from Caradima et al. (2019)
# used to calculate the insecticide application rate (IAR)

d_index_factors <- c(grains = .03, 
                     corn = .01,
                     beets = .07,
                     rapeseed = 1.82,
                     potatoes = .44,
                     vineyards = .37,
                     vegetables = 2.66,
                     legumes = .38,
                     orchards = 3.1) %>% 
  data.frame(index_factor = .) %>% 
  rownames_to_column("culture")

# records data -----------------------------------------------------------------.
# (to determine which squares were included in analyses)

d_records_butter <- fread("Data/d_records_butterflies.csv")
d_records_odo <- fread("Data/d_records_dragonflies.csv")
d_records_ortho<- fread("Data/d_records_grasshoppers.csv")

d_records <- d_records_butter %>% 
  mutate(group = "butterflies") %>% 
  bind_rows(d_records_odo %>% 
              mutate(group = "dragon- and damselflies"),
            d_records_ortho %>% 
              mutate(group = "grasshoppers"))

d_km2_ID <- d_records %>% 
  mutate(CX_km2 = as.integer(substr(km2_ID, 1, 7)),
         CY_km2 = as.integer(substr(km2_ID, 11, 18))) %>% 
  select(km2_ID, CX_km2, CY_km2) %>% 
  distinct()

# zone data for squares km2 ----------------------------------------------------.
d_zone <- fread("Data/d_zone.csv")

d_zone <- d_zone %>% 
  mutate(CX_km2 = as.numeric(substr(km2_ID, 1, 7)),
         CY_km2 = as.numeric(substr(km2_ID, 11, 18)),
         CX_km2_center = CX_km2 + 500,
         CY_km2_center = CY_km2 + 500)

# Arealstatistik data ----------------------------------------------------------.
# Switzerland-wide land-use data
# https://www.bfs.admin.ch/bfs/de/home/statistiken/kataloge-datenbanken/publikationen.assetdetail.256981.html

d_nolu <- fread("Data/ag-b-00.03-37-nolu04/AREA_NOLU04_46_161114.csv")
d_nolu <- d_nolu %>% 
  mutate(X = X + 2000000,
         Y = Y + 1000000) %>% 
  mutate(CX_km2 = floor(X / 1000) * 1000,
         CY_km2 = floor(Y / 1000) * 1000,
         km2_ID = paste(CX_km2, CY_km2, sep = " / "))

# municipalities of Switzerland as of 1.1.2021 ---------------------------------.

# available from https://www.bfs.admin.ch/bfs/de/home/dienstleistungen/geostat/geodaten-bundesstatistik/administrative-grenzen/generalisierte-gemeindegrenzen.assetdetail.15724821.html
sf_municipalities <- st_read(dsn = "Data/ag-b-00.03-875-gg21/ggg_2021-LV95/shp/g1g21_01012021.shp") %>% 
  st_transform(crs = 2056)

# edit bioclimate zone data ####################################################
################################################################################.

# create sf object
radius <- 500

sf_zone <- d_zone %>% 
  mutate(xlow = CX_km2_center - radius,
         ylow = CY_km2_center - radius,
         xhigh = CX_km2_center + radius,
         yhigh = CY_km2_center + radius,
         X1 = xlow, Y1 = ylow,
         X2 = xlow, Y2 = yhigh,
         X3 = xhigh, Y3 = yhigh,
         X4 = xhigh, Y4 = ylow,
         X5 = xlow, Y5 = ylow) %>% 
  pivot_longer(c(X1, Y1, X2, Y2, X3, Y3, X4, Y4, X5, Y5)) %>% 
  mutate(dim = substr(name, 1, 1),
         nr = substr(name, 2, 2)) %>% 
  select(-name) %>% 
  spread(dim, value) %>% 
  sf_polygon(polygon_id = "km2_ID", x = "X", y = "Y") 

sf_zone$zone <- d_zone$zone
st_crs(sf_zone) <- 2056 # define coordinate reference system

# Determine which area proportions of each commune belong to which zone ########
################################################################################.

# use data from Arealstastik to do so

# create sf of Arealstatistik data
sf_nolu <- d_nolu %>% 
  mutate(ha_ID = paste(X, Y, sep = " / ")) %>% 
  mutate(xlow = X,
         ylow = Y,
         xhigh = X + 100,
         yhigh = Y + 100,
         X1 = xlow, Y1 = ylow,
         X2 = xlow, Y2 = yhigh,
         X3 = xhigh, Y3 = yhigh,
         X4 = xhigh, Y4 = ylow,
         X5 = xlow, Y5 = ylow) %>% 
  pivot_longer(c(X1, Y1, X2, Y2, X3, Y3, X4, Y4, X5, Y5)) %>% 
  mutate(dim = substr(name, 1, 1),
         nr = substr(name, 2, 2)) %>% 
  select(-name) %>% 
  spread(dim, value) %>% 
  sf_polygon(polygon_id = "ha_ID", x = "X", y = "Y") 

st_crs(sf_nolu) <- 2056 

# intersect with municipalities
sf_municipalities_areal <- st_intersection(sf_nolu, sf_municipalities)  
sf_municipalities_areal$area <- st_area(sf_municipalities_areal) # determine area

# convert to data frame
d_municipalities_areal <- as.data.frame(sf_municipalities_areal) %>% 
  select(-geometry) %>% 
  mutate(area = as.numeric(area))

# determine share of agricultural area (AgrArea) for different zone categories per commune
d_municipalities_AgrArea <-
  d_municipalities_areal %>% 
  left_join(d_nolu %>% 
              left_join(d_zone, by = "km2_ID") %>% 
              mutate(ha_ID = paste(X, Y, sep = " / ")),
            by = "ha_ID") %>% 
  left_join(d_km2_ID %>% 
              mutate(in_records = T),
            by = "km2_ID") %>% 
  # AgrArea categories (including buildings, orchards, crops, grasslands, ...)
  mutate(AgrArea85 = LU85_46 %in% c(107, 201, 202, 203, 221, 222, 223), 
         AgrArea97 = LU97_46 %in% c(107, 201, 202, 203, 221, 222, 223),
         AgrArea09 = LU09_46 %in% c(107, 201, 202, 203, 221, 222, 223)) %>%
  group_by(GMDNR, GMDNAME, zone, in_records) %>% 
  summarise(area_tot = sum(area),
            AgrArea85 = sum(as.numeric(AgrArea85) * area),
            FJ85_mean = mean(FJ85),
            AgrArea97 = sum(as.numeric(AgrArea97) * area),
            FJ97_mean = mean(FJ97),
            AgrArea09 = sum(as.numeric(AgrArea09) * area),
            FJ09_mean = mean(FJ09),
            .groups = "drop") %>% 
  group_by(GMDNR, GMDNAME) %>% 
  mutate(prop_AgrArea85 =  AgrArea85 / sum(AgrArea85),
         prop_AgrArea97 =  AgrArea97 / sum(AgrArea97),
         prop_AgrArea09 =  AgrArea09 / sum(AgrArea09)) %>% 
  ungroup()

# municipalities in which no records were available (used later)
d_municipalities_no_sample <- d_municipalities_AgrArea %>% 
  group_by(GMDNR, GMDNAME) %>% 
  summarise(sampled = !all(is.na(in_records)),
            .groups = "drop") %>% 
  filter(!sampled)


# exclude observations from squares for which no records were available
d_municipalities_AgrArea <- d_municipalities_AgrArea %>% 
  filter(in_records) %>% 
  select(-in_records) %>% 
  # very few cases with no AgrArea at al:
  mutate(prop_AgrArea85 = ifelse(is.na(prop_AgrArea85), 0, prop_AgrArea85),
         prop_AgrArea97 = ifelse(is.na(prop_AgrArea97), 0, prop_AgrArea97),
         prop_AgrArea09 = ifelse(is.na(prop_AgrArea09), 0, prop_AgrArea09))


# temporal resolution of AgrArea prop ------------------------------------------.

# assume no changes pre 85 and post 09
# this should have a very minor impact, the proportion between regions are
# not changing drastically on such time scales

d_municipalities_AgrArea_res <- data.frame()
for (i in 1:nrow(d_municipalities_AgrArea)){
  d_target <- d_municipalities_AgrArea[i, ]
  
  # assume no change pre 85 ----------------------------------------------------.
  
  d_out <- data.frame(year = 1950:floor(d_target$FJ85_mean),
                      prop_AgrArea = d_target$prop_AgrArea85)
  
  # interpolate 85 - 97 --------------------------------------------------------.
  
  mod1 <- data.frame(prop_AgrArea = c(d_target$prop_AgrArea85, d_target$prop_AgrArea97),
                     year = c(d_target$FJ85_mean, d_target$FJ97_mean)) %>% 
    lm(prop_AgrArea ~ year, data = .)
  
  pred1 <- predict(mod1, newdata = data.frame(year = seq(floor(d_target$FJ85_mean) + 1,
                                                         floor(d_target$FJ97_mean), 1)))
  
  d_out <- data.frame(prop_AgrArea = pred1,
                      year = seq(floor(d_target$FJ85_mean) + 1,
                                 floor(d_target$FJ97_mean), 1)) %>% 
    bind_rows(d_out, .)
  
  # interpolate 97 - 09 --------------------------------------------------------.
  
  mod2 <- data.frame(prop_AgrArea = c(d_target$prop_AgrArea97, d_target$prop_AgrArea09),
                     year = c(d_target$FJ97_mean, d_target$FJ09_mean)) %>% 
    lm(prop_AgrArea ~ year, data = .)
  
  pred2 <- predict(mod2, newdata = data.frame(year = seq(floor(d_target$FJ97_mean) + 1,
                                                         floor(d_target$FJ09_mean), 1)))
  
  d_out <- data.frame(prop_AgrArea = pred2,
                      year = seq(floor(d_target$FJ97_mean) + 1,
                                 floor(d_target$FJ09_mean), 1)) %>% 
    bind_rows(d_out, .)
  
  # assume no change past 09 ---------------------------------------------------.
  
  d_out <- data.frame(year = (floor(d_target$FJ09_mean) + 1):2020,
                      prop_AgrArea = d_target$prop_AgrArea09) %>% 
    bind_rows(d_out, .)
  
  # output ---------------------------------------------------------------------.
  
  d_municipalities_AgrArea_res <- data.frame(d_target[, c("GMDNR", "GMDNAME", "zone")], 
                                   d_out) %>% 
    bind_rows(d_municipalities_AgrArea_res, .)
}



# assemble data per bioclimatic zone ###########################################
################################################################################.

d_AgrArea_reg <- d_AgrArea %>% 
  filter(!BFS_Nr %in% d_municipalities_no_sample$GMDNR) %>% 
  left_join(d_municipalities_AgrArea_res, by = c(Jahr = "year",
                                       BFS_Nr = "GMDNR")) %>% 
  mutate(value = value * prop_AgrArea) %>% 
  group_by(zone, Jahr) %>% 
  summarise(AgrArea = sum(value),
            .groups = "drop")


d_LSU_reg <- d_LSU %>% 
  filter(!BFS_Nr %in% d_municipalities_no_sample$GMDNR) %>% 
  left_join(d_municipalities_AgrArea_res, by = c(Jahr = "year",
                                       BFS_Nr = "GMDNR")) %>% 
  mutate(LSU = LSU * prop_AgrArea) %>% 
  group_by(zone, Jahr) %>% 
  summarise(LSU = sum(LSU),
            .groups = "drop")


d_grasslands_reg <- d_grasslands %>% 
  filter(!BFS_Nr %in% d_municipalities_no_sample$GMDNR) %>% 
  left_join(d_municipalities_AgrArea_res, by = c(Jahr = "year",
                                       BFS_Nr = "GMDNR")) %>% 
  mutate(value = value * prop_AgrArea) %>% 
  group_by(zone, Jahr) %>% 
  summarise(AgrArea = sum(value),
            .groups = "drop")

# add pre 1975 data with an approximation
# assume no change in the proportion of grasslands in the total agricultural area pre 1975
d_grasslands_reg <- d_AgrArea_reg %>% 
  left_join(d_grasslands_reg, by = c("zone", "Jahr"),
            suffix = c("", "_GL")) %>% 
  mutate(prop_GL = AgrArea_GL / AgrArea) %>% 
  group_by(zone) %>% 
  mutate(prop_GL = ifelse(is.na(prop_GL), prop_GL[Jahr == 1975], prop_GL),
         AgrArea_GL = AgrArea * prop_GL) %>% 
  select(-c(AgrArea, prop_GL))


d_culturess_comb_reg <- d_culturess %>% 
  filter(!BFS_Nr %in% d_municipalities_no_sample$GMDNR) %>% 
  left_join(d_municipalities_AgrArea_res, by = c(Jahr = "year",
                                       BFS_Nr = "GMDNR")) %>% 
  # weight Area by the proportion of agricultural land in area having species records
  # in that region. This is an approximative approach
  mutate(Area = Area * prop_AgrArea) %>%
  group_by(zone, Jahr, culture) %>% 
  summarise(Area = sum(Area),
            .groups = "drop") %>% 
  left_join(d_index_factors, by = "culture") %>% 
  filter(!is.na(index_factor)) %>% 
  mutate(index_contr = Area * index_factor)

# Calculate land-use change variables ##########################################
################################################################################.

# time intervals in which the variables are estimated
timestep <- 5; set.seed(12) 
timestep <- 10; set.seed(17) 

d_as_gam_smry <- data.frame() # summary data frame (will be the main output)

for (region_i in unique(d_AgrArea_reg$zone)){
  # AgrArea prop ---------------------------------------------------------------.
  
  target_AgrArea_prop <- d_AgrArea_reg %>% 
    filter(zone == region_i) %>% 
    left_join(d_municipalities_AgrArea %>% 
                group_by(zone) %>% 
                summarise(area_tot = sum(area_tot) / 10000, .groups = "drop"),
              by = "zone") %>% 
    mutate(AgrArea_prop = AgrArea / area_tot)
  
  mod_AgrArea_prop <- stan_gamm4(AgrArea_prop ~  s(Jahr), data = target_AgrArea_prop,
                       chains = 4, iter = 2000)
  
  pred_AgrArea_prop <- posterior_epred(mod_AgrArea_prop, newdata = data.frame(Jahr = min(target_AgrArea_prop$Jahr):max(target_AgrArea_prop$Jahr)))
  
  d_pred_AgrArea_prop <- as.data.frame(pred_AgrArea_prop)
  d_pred_AgrArea_prop$run <- 1:nrow(d_pred_AgrArea_prop)
  d_pred_AgrArea_prop <- d_pred_AgrArea_prop %>%
    gather(Jahr, value, -run) %>%
    mutate(Jahr = as.integer(Jahr),
           Jahr = Jahr + min(target_AgrArea_prop$Jahr) - 1,
           zone = region_i,
           var = "AgrArea_prop")
    
  # AgrArea -------------------------------------------------------------------------.
  target_AgrArea <- d_AgrArea_reg %>% 
    filter(zone == region_i)
  
  mod_AgrArea <- stan_gamm4(AgrArea ~  s(Jahr), data = target_AgrArea,
                       chains = 4, iter = 2000)
  
  pred_AgrArea <- posterior_epred(mod_AgrArea, newdata = data.frame(Jahr = min(target_AgrArea$Jahr):max(target_AgrArea$Jahr)))
  
  d_pred_AgrArea <- as.data.frame(pred_AgrArea)
  d_pred_AgrArea$run <- 1:nrow(d_pred_AgrArea)
  d_pred_AgrArea <- d_pred_AgrArea %>%
    gather(Jahr, value, -run) %>%
    mutate(Jahr = as.integer(Jahr),
           Jahr = Jahr + min(target_AgrArea$Jahr) - 1,
           zone = region_i,
           var = "AgrArea")
  
  # grassland area -------------------------------------------------------------.
  target_GL <- d_grasslands_reg %>% 
    filter(zone == region_i)
  
  mod_GL <- stan_gamm4(AgrArea_GL ~  s(Jahr), data = target_GL,
                       chains = 4, iter = 2000)
  
  pred_GL <- posterior_epred(mod_GL, newdata = data.frame(Jahr = min(target_GL$Jahr):max(target_GL$Jahr)))
  
  d_pred_GL <- as.data.frame(pred_GL)
  d_pred_GL$run <- 1:nrow(d_pred_GL)
  d_pred_GL <- d_pred_GL %>%
    gather(Jahr, value, -run) %>%
    mutate(Jahr = as.integer(Jahr),
           Jahr = Jahr + min(target_GL$Jahr) - 1,
           zone = region_i,
           var = "GL")
  
  # LSU-------------------------------------------------------------------------.
  
  target_LSU = d_LSU_reg %>% 
    filter(zone == region_i)
  
  mod_LSU <- stan_gamm4(LSU ~  s(Jahr), data = target_LSU,
                       chains = 4, iter = 2000)
  
  pred_LSU <- posterior_epred(mod_LSU, newdata = data.frame(Jahr = min(target_LSU$Jahr):max(target_LSU$Jahr)))
  
  d_pred_LSU <- as.data.frame(pred_LSU)
  d_pred_LSU$run <- 1:nrow(d_pred_LSU)
  d_pred_LSU <- d_pred_LSU %>%
    gather(Jahr, value, -run) %>%
    mutate(Jahr = as.integer(Jahr),
           Jahr = Jahr + min(target_LSU$Jahr) - 1,
           zone = region_i,
           var = "LSU")

  # LSU / GL -------------------------------------------------------------------.
  
  # remove 1955 which is not covered in LSU
  pred_GL_sub <- pred_GL[, -1]
  
  pred_LSU_GL <- pred_LSU / pred_GL_sub
  
  d_pred_LSU_GL <- as.data.frame(pred_LSU_GL)
  d_pred_LSU_GL$run <- 1:nrow(d_pred_LSU_GL)
  d_pred_LSU_GL <- d_pred_LSU_GL %>%
    gather(Jahr, value, -run) %>%
    mutate(Jahr = as.integer(Jahr),
           Jahr = Jahr + 1955,
           zone = region_i,
           var = "LSU/GL")

  # IAR ---------------------------------------------------------------------.
  
  target_IAR <- d_culturess_comb_reg %>% 
    filter(zone == region_i) %>% 
    group_by(Jahr, zone) %>% 
    summarise(IAR = sum(index_contr),
              .groups = "drop") %>% 
    left_join(d_municipalities_AgrArea %>% 
                group_by(zone) %>% 
                summarise(area_tot = sum(area_tot) / 10000, .groups = "drop"),
              by = "zone") %>% 
    mutate(IAR = IAR / area_tot)
  
  mod_IAR <- stan_gamm4(IAR ~  s(Jahr), data = target_IAR,
                            chains = 4, iter = 2000)
  
  pred_IAR <- posterior_epred(mod_IAR, newdata = data.frame(Jahr = min(target_IAR$Jahr):max(target_IAR$Jahr)))
  
  d_pred_IAR <- as.data.frame(pred_IAR)
  d_pred_IAR$run <- 1:nrow(d_pred_IAR)
  d_pred_IAR <- d_pred_IAR %>%
    gather(Jahr, value, -run) %>%
    mutate(Jahr = as.integer(Jahr),
           Jahr = Jahr + min(target_IAR$Jahr) - 1,
           zone = region_i,
           var = "IAR")

  # ----------------------------------------------------------------------------.
  # determine state variables --------------------------------------------------.
  # ----------------------------------------------------------------------------.
  
  # cut to years of interest (only 1980 and later used in main analyses)
  pred_AgrArea_prop <- pred_AgrArea_prop[, 1955:2020 %in% c(1970:2020)]
  pred_LSU_GL <- pred_LSU_GL[, 1956:2020 %in% c(1970:2020)]
  pred_IAR <- pred_IAR[, 1969:2020 %in% c(1970:2020)]
  
  for (i in 1:(50/timestep)){
    
    # AgrArea_prop ------------------------------------------------------------------.
    
    # change
    change_mean <- mean(pred_AgrArea_prop[, i * timestep + 1] - pred_AgrArea_prop[, (i-1) * timestep + 1])
    
    # mean
    mean_mean <- mean(apply(pred_AgrArea_prop[, ((i-1) * timestep + 1):(i * timestep + 1)], 1, mean))
    
    d_as_gam_smry <- data.frame(var = "AgrArea_prop",
                                zone = region_i,
                                year_start = 1969 + (i-1) * timestep + 1,
                                year_end = 1969 + i * timestep + 1,
                                change_mean, mean_mean) %>% 
      bind_rows(d_as_gam_smry, .)
    
    
    # LSU / GL -----------------------------------------------------------------.
    
    # change
    change_mean <- mean(pred_LSU_GL[, i * timestep + 1] - pred_LSU_GL[, (i-1) * timestep + 1])
    
    # mean
    mean_mean <- mean(apply(pred_LSU_GL[, ((i-1) * timestep + 1):(i * timestep + 1)], 1, mean))
    
    d_as_gam_smry <- data.frame(var = "LSU/GL",
                                zone = region_i,
                                year_start = 1969 + (i-1) * timestep + 1,
                                year_end = 1969 + i * timestep + 1,
                                change_mean, 
                                mean_mean) %>% 
      bind_rows(d_as_gam_smry, .)
    

    # IAR -------------------------------------------------------------------.
    
    # change
    change_mean <- mean(pred_IAR[, i * timestep + 1] - pred_IAR[, (i-1) * timestep + 1])
    
    # mean
    mean_mean <- mean(apply(pred_IAR[, ((i-1) * timestep + 1):(i * timestep + 1)], 1, mean))
    
    d_as_gam_smry <- data.frame(var = "IAR",
                                zone = region_i,
                                year_start = 1969 + (i-1) * timestep + 1,
                                year_end = 1969 + i * timestep + 1,
                                change_mean, 
                                mean_mean) %>% 
      bind_rows(d_as_gam_smry, .)
  }
}

# export output
saveRDS(d_as_gam_smry, paste0("Data/d_as_gam_smry_", timestep, ".rds"))
