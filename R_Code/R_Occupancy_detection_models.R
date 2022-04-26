# Initialise system ############################################################
################################################################################.

rm(list = ls()); graphics.off()
Sys.setlocale('LC_ALL','C.UTF-8')
options(max.print = 500)

library(tidyverse)
library(cmdstanr)
library(parallel)
library(data.table)
library(bayestestR)
library(posterior)

source('R_Code/f_occ_det.R')
stan_mod <- cmdstan_model('Stan_Code/Stan_occ_det_cmdstan.stan')

options(mc.cores = 4)

# information of the model run
sp_index_i <- as.numeric(commandArgs(trailingOnly = TRUE)[1]) # species index
iter <- as.numeric(commandArgs(trailingOnly = TRUE)[2]) # number of iterations to be run

# read either butterfly, dragonfly or grasshopper records data
d_records <- fread('Data/d_records.csv')

# species to be analysed
sel_species <- d_records %>% 
  group_by(species) %>% 
  summarise(n_year = length(unique(year))) %>% 
  filter(n_year > (max(d_records$year) - min(d_records$year)) * .25,
         # exclude aggregates for which better information on species level is available
         !species %in% c("Aricia agestis aggr.", "Colias hyale aggr.",
                        "Leptidea sinapis aggr.", "Pieris napi aggr.",
                        "Pieris rapae aggr.", "Boloria pales aggr.",
                        "Coenonympha gardetta aggr.",
                        "Cordulegaster bidentata/boltonii",
                        "Sympetrum striolatum/vulgatum")) %>% 
  select(species) %>% 
  deframe()

sp_i <- sel_species[sp_index_i]

# projects with restricted species scope
d_rest_projects <- fread('Data/projects_restricted.csv')

# experts
d_experts <- fread('Data/experts.csv')

# biogeographic regions data (per square km2)
d_biogeo <- fread('Data/biogeographic_regions.csv')

# mean elevation above sealevel data per square km2
d_elevation <- fread('Data/mean_elevation_km2.csv')

# data on bioclimatic zones (per square km2)
d_zone <- fread('Data/bioclimatic_zone.csv')

# sites data set (combination of square and observation year)
d_sites <- expand.grid(km2_ID = unique(d_records$km2_ID),
                       year = seq(min(d_records$year), max(d_records$year), 1)) %>% 
  arrange(km2_ID, year)

# range of year days to be included (phenology)
daylimits <- d_records %>% 
  filter(species == sp_i) %>% 
  select(yday) %>% 
  distinct() %>% # weight all days equally
  summarise(lower = quantile(yday, .05),
            upper = quantile(yday, .95))

# visits data set
d_visits <-
  d_records %>% 
  filter(yday >= daylimits$lower,
         yday <= daylimits$upper) %>% 
  # add information for restricted projects
  left_join(d_rest_projects, by = 'PROJET') %>% 
  # exclude records from restricted projects that did not aim at the focal species:
  rowwise() %>% 
  filter((!species %in% strsplit(focal_species, ' \\| ')[[1]]) |
           any(strsplit(focal_species, ' \\| ')[[1]] == sp_i)) %>% 
  ungroup() %>% 
  group_by(km2_ID, year, yday, visit_ID) %>% 
  summarise(pres = ifelse(any(species == sp_i), 1, 0),
            # define list length factor
            list_length = length(unique(species)),
            list_length_cat = cut(list_length, breaks = c(0, 1, 3, Inf),
                                  labels = c('single', 'short', 'long')),
            # determine source category
            source = ifelse(sp_i %in% unlist(strsplit(focal_species, ' \\| ')),
                            'Project_targeted', # targeted project
                            ifelse(any(PROJET == 'LRPAP'), 'Project_RL', # red list project
                                   ifelse(any(PROJET != ''), 'Project',
                                          ifelse(any(LEG %in% d_experts$LEG),
                                                 'CitSc_expert', 'CitSc')))), # expert or naturalist observation
            .groups = 'drop') %>% 
  mutate(km2_ID = as.factor(km2_ID),
         fact_year = as.factor(year),
         yday_z = as.numeric(scale(yday)),
         source = as.factor(source)) %>% 
  arrange(km2_ID, year, yday)


# add environmental data #######################################################
################################################################################.

d_sites <- d_sites %>% 
  left_join(d_biogeo %>% 
              select(km2_ID, biogeo12), by = 'km2_ID') %>% 
  left_join(d_zone, by = 'km2_ID') %>% 
  mutate_at(vars(biogeo12, zone), as.factor) %>% 
  left_join(d_elevation %>% 
              filter(km2_ID %in% d_sites$km2_ID) %>% 
              select(km2_ID, elevation) %>% 
              mutate(elevation_z = as.numeric(scale(elevation))), by = 'km2_ID')

# subset based on zones ########################################################
################################################################################.

# only consider zones, where species was reported at least once

sel_zone <- d_visits %>% 
  filter(pres == 1) %>% 
  left_join(d_sites, by = 'km2_ID') %>%
  select(zone) %>% 
  distinct()

d_sites <- d_sites %>% 
  filter(zone %in% sel_zone$zone) %>% 
  droplevels()

d_visits <- d_visits %>% 
  filter(km2_ID %in% d_sites$km2_ID)


d_sites <- d_sites %>% 
  arrange(km2_ID, year)

d_visits <- d_visits %>% 
  arrange(km2_ID, year)


# run MCMC chains ##############################################################
################################################################################.

f_occ_det(d_sites = d_sites, d_visits = d_visits, 
          formula_occ = ~ elevation_z + I(elevation_z^2) + (1 | biogeo12) + (1 | km2_ID),  # formula occurrence probability
          formula_det = ~ yday_z + I(yday_z^2) + list_length_cat + source + (1 | fact_year), # formula detection probability
          var_site = 'km2_ID', 
          var_year = 'year', 
          var_presabs = 'pres',
          var_area = 'zone',
          stan_mod = stan_mod,
          iter_warmup = floor(iter / 2),
          iter_sampling = iter - floor(iter / 2),
          chains = 4,
          output_dir = 'Output',
          output_basename = gsub(' ', '_', sp_i))

# post-sampling processing #####################################################
################################################################################.

# extract occupancy estimates --------------------------------------------------.

dir <- 'Output'


files <- list.files(paste0(dir, gsub(' ', '_', sp_i)),
                    full.names = T)

fit <- read_cmdstan_csv(files,
                        variables = c('z'), 
                        format = 'draws_list')

d_mcmc_z <- fit$post_warmup_draws[, grepl(paste(paste0('^z'), 
                                           collapse = '|'), 
                                     names(x)), drop = F] %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("var")


# aggregate per zone and year --------------------------------------------------.

d_z_yearmean_reg <- data.frame()
for (i in unique(d_sites$year)){
  for (j in unique(d_sites$zone)) {
    # site*year combinations of year i
    sel <- which(d_sites$year == i & d_sites$zone == j) 
    
    d_z_yearmean_reg <- data.frame(year = i, zone = j,
                                   n_squares = length(sel),
                                   iter = 1:ncol(d_mcmc_z[, -1]),
                                   z_mean = apply(d_mcmc_z[sel, -1], 2, mean)) %>%
      bind_rows(d_z_yearmean_reg, .)
  }
}

# create output ----------------------------------------------------------------.

d_z_yearmean_reg %>% 
  mutate(iter = paste0("occ_m_i", formatC(iter, width = 4, flag = "0")),
         species = sp_i) %>% 
  pivot_wider(names_from = iter, values_from = z_mean) %>% 
  select(species, year, zone, n_squares, everything()) %>% 
  fwrite(file = "Data/occ_means.csv",
         append = T)
