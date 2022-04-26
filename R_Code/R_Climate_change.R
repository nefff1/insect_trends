library(raster)
library(tidyverse); theme_set(theme_classic())
library(sf)
library(sfheaders)
library(parallel)
library(data.table)
library(lubridate)
library(rstanarm)
options(mc.cores = 4)

select <- dplyr::select

# Read raw data ################################################################
################################################################################.

# Data frame containing monthly mean temperature (TabsM) and precipitation sum (RhiresM)
# Spatially resolved at a 1.25 degree grid (CX and CY)
# covers years 1961-2020
# Data from MeteoSwiss
d_climate <- fread("Data/d_climate.csv")

# add number of days per month
d_climate <- d_climate %>% 
  mutate(date_dummy = paste(year, as.numeric(month), 1, sep = "-"),
         ndays = days_in_month(as.Date(date_dummy)))

# bioclimatic regions
d_zone <- fread('Data/d_zone.csv')

# records data
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

# Build yearly variables #######################################################
################################################################################.

# aggregate data per 3 months
d_quartermean <- data.frame()
for (month_i in 1:10){
  d_quartermean <- d_climate %>% 
    filter(as.numeric(month) %in% month_i:(month_i + 2)) %>% 
    group_by(year, CX, CY) %>% 
    summarise(Tmean_quarter = mean(TabsM),
              Psum_quarter = sum(RhiresM),
              .groups = "drop") %>% 
    mutate(month = month.abb[month_i]) %>% 
    bind_rows(d_quartermean, .)
}
d_quartermean <- d_quartermean %>% 
  mutate(month = factor(month, levels = month.abb))

# aggregate data per grid cell per year
d_climate_smry <-
  d_climate %>% 
  left_join(d_quartermean,  by = c("month", "year", "CX", "CY")) %>% 
  group_by(year, CX, CY) %>%
  summarise(T_all = sum(TabsM * ndays) / sum(ndays),
            T_seas = sd(TabsM),
            P_warmest_quarter = Psum_quarter[Tmean_quarter == max(Tmean_quarter, na.rm = T) & !is.na(Psum_quarter)],
            .groups = "drop") 

rm(list = c("d_quartermean", "d_climate")) # cleaning

# aggregate by bioclimatic zone ################################################
################################################################################.

# approach: check whether there is a sampling square (1km x 1km) intersecting with a climate square

d_km2_ID <- d_records %>% 
  mutate(CX_km2 = as.integer(substr(km2_ID, 1, 7)),
         CY_km2 = as.integer(substr(km2_ID, 11, 18))) %>% 
  select(km2_ID, CX_km2, CY_km2) %>% 
  distinct()
sf_km2_ID <- d_km2_ID %>% 
  mutate(xlow = CX_km2,
         ylow = CY_km2,
         xhigh = CX_km2 + 1000,
         yhigh = CY_km2 + 1000,
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
st_crs(sf_km2_ID) <- 2056 # define coordinate reference system

#  approximate grid cell width of climate data
radius_x <- 1600 / 2
radius_y <- 2320 / 2

# create sf polygons
sf_climate_smry_sub <-
  d_climate_smry %>% 
  filter(year == 2000) %>% # dummy year
  mutate(ID = paste(CX, CY),
         xlow = CX - radius_x,
         ylow = CY - radius_y,
         xhigh = CX + radius_x,
         yhigh = CY + radius_y,
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
  sf_polygon(polygon_id = "ID", x = "X", y = "Y") 

st_crs(sf_climate_smry_sub) <- 2056 # define coordinate reference system

# select those squares that intersect with a sampled km2_ID
sel <- st_intersects(sf_climate_smry_sub, st_union(sf_km2_ID), prepared = F)
sel <- sapply(sel, function(x) length(x) > 0)

sel_ids <- paste(round(d_climate_smry$CX[d_climate_smry$year == 2000][sel]),
                 round(d_climate_smry$CY[d_climate_smry$year == 2000][sel]),
                 sep = " / ")

# add bioclimatic zone info and summarise per bioclimatic zone

d_climate_smry_reg <- d_climate_smry %>% 
  mutate(ID = paste(round(CX), round(CY), sep = " / ")) %>% 
  filter(ID %in% sel_ids) %>% 
  select(-ID) %>% 
  mutate(CX_km2 = floor(CX / 1000) * 1000,
         CY_km2 = floor(CY / 1000) * 1000,
         km2_ID = paste(CX_km2, CY_km2, sep = " / ")) %>% 
  left_join(d_zone, by = "km2_ID") %>% 
  filter(!is.na(zone)) %>% # border region, should not have a significant influence on means
  select(-c(km2_ID, CX_km2, CY_km2, CX, CY)) %>% 
  group_by(year, zone) %>% 
  summarise_all(.funs = ~ mean(.)) %>%
  ungroup()

# Calculate climate change variables ###########################################
################################################################################.

vars <- c("T_all", "T_seas", "P_warmest_quarter")

# 5 year intervals -------------------------------------------------------------.
# always include also preceding 5 years

set.seed(21)
timestep <- 5
steps <- seq(1980, 2019, timestep)

d_climate_lm_smry <- data.frame()
l_plots <- list()

for (var in vars){
  out <- data.frame()
  d_pred_comb <- data.frame()
  for (region_i in unique(d_climate_smry_reg$zone)){
    for (step_i in steps){
      target <- d_climate_smry_reg %>% 
        filter(year >= step_i - 5, year <= step_i + timestep) %>%
        rename(target_var := !! enquo(var)) %>% 
        group_by(zone) %>% 
        ungroup() %>% 
        filter(zone == region_i)
      
      mod <- stan_glm(target_var ~ year, data = target,
                     chains = 4, iter = 2000, 
                     family = "gaussian",
                     prior_intercept = normal(mean(deframe(target[, "target_var"])), 5, autoscale = T),
                     prior = normal(0, 5, autoscale = T))
      
      pred <- posterior_epred(mod)
      
      d_pred <- as.data.frame(pred)
      d_pred$run <- 1:nrow(d_pred)
      d_pred <- d_pred %>% 
        gather(year, value, -run) %>%
        mutate(year = as.integer(year),
               year = year + step_i - 5 - 1)
      d_pred_comb <- d_pred %>% 
        mutate(zone = region_i,
               step = step_i) %>% 
        bind_rows(d_pred_comb, .)
      
      out <- d_pred %>% 
        group_by(run) %>% 
        summarise(change = value[year == step_i + 5] - value[year == step_i - 5],
                  .groups = "drop") %>% 
        summarise(change_mean = mean(change)) %>% 
        mutate(var = var,
               zone = region_i,
               year_start = step_i - 5,
               year_end = step_i + 5) %>% 
        bind_rows(out, .)

    }
  }
  
  d_climate_lm_smry <- d_climate_lm_smry %>% 
    bind_rows(out, .)
}

# export data
saveRDS(d_climate_lm_smry, paste0("Data/d_climate_lm_smry_", timestep, ".rds"))


# 10 years intervals -----------------------------------------------------------.
# only include the years of interest

set.seed(217)
timestep <- 10
steps <- seq(1980, 2019, timestep)

d_climate_lm_smry <- data.frame()
l_plots <- list()

for (var in vars){
  out <- data.frame()
  d_pred_comb <- data.frame()
  for (region_i in unique(d_climate_smry_reg$zone)){
    for (step_i in steps){
      target <- d_climate_smry_reg %>% 
        filter(year >= step_i, year <= step_i + timestep) %>%
        rename(target_var := !! enquo(var)) %>% 
        group_by(zone) %>% 
        ungroup() %>% 
        filter(zone == region_i)
      
      mod <- stan_glm(target_var ~ year, data = target,
                      chains = 4, iter = 2000, 
                      family = "gaussian",
                      prior_intercept = normal(mean(deframe(target[, "target_var"])), 5, autoscale = T),
                      prior = normal(0, 5, autoscale = T))
      
      pred <- posterior_epred(mod)
      
      d_pred <- as.data.frame(pred)
      d_pred$run <- 1:nrow(d_pred)
      d_pred <- d_pred %>% 
        gather(year, value, -run) %>%
        mutate(year = as.integer(year),
               year = year + step_i - 1)
      d_pred_comb <- d_pred %>% 
        mutate(zone = region_i,
               step = step_i) %>% 
        bind_rows(d_pred_comb, .)
      
      
      out <- d_pred %>% 
        group_by(run) %>% 
        summarise(change = value[year == step_i + timestep] - value[year == step_i],
                  .groups = "drop") %>% 
        summarise(change_mean = mean(change)) %>% 
        mutate(var = var,
               zone = region_i,
               year_start = step_i,
               year_end = step_i + timestep) %>% 
        bind_rows(out, .)
      
    }
  }

  d_climate_lm_smry <- d_climate_lm_smry %>% 
    bind_rows(out, .)
}


# export data
saveRDS(d_climate_lm_smry, paste0("Data/Output/d_climate_lm_smry_", timestep, ".rds"))
