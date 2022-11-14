# load packages and data #######################################################
################################################################################.

library(data.table)
library(tidyverse); theme_set(theme_classic())
library(bayestestR)
library(parallel)
library(cowplot)
library(sf)
library(sfheaders)
library(lmerTest)
library(rstan)
library(shinystan)
library(bayesplot)
library(patchwork)
library(grid)
library(rstanarm)
library(ggnewscale)
library(itsadug)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

extract <- rstan::extract

# temperature niche ------------------------------------------------------------.
# Output of R_Temperature_niche.R

d_temperature_niche_butter <- readRDS("Data/d_temperature_niche_butterflies.rds")
d_temperature_niche_odo <- readRDS("Data/d_temperature_niche_dragonflies.rds")
d_temperature_niche_ortho <- readRDS("Data/d_temperature_niche_grasshoppers.rds")

d_temperature_niche <- d_temperature_niche_butter %>% 
  mutate(group = "butterflies") %>% 
  bind_rows(d_temperature_niche_odo %>% 
              mutate(group = "dragonflies"), 
            d_temperature_niche_ortho %>% 
              mutate(group = "grasshoppers"))

# specialisation ---------------------------------------------------------------.
# Raw data available from http://www.cscf.ch/cscf/de/home/projekte/fauna-indicativa.html
# Processed data available from supplementary material

d_spec <- readRDS("Data/habitat_spec_faunind.rds")

# exclude species that were not analysed
d_spec <- d_spec %>% 
  filter(species %in% d_temperature_niche$species)

d_spec <- d_spec %>% 
  group_by(group) %>% 
  mutate(spec_index = habitat_spec_index - min(habitat_spec_index),
         spec_index = spec_index / max(spec_index)) %>% 
  ungroup()

# Records data -----------------------------------------------------------------.
# data are protected and are not available due to data privacy laws 

d_records_butter <- fread("Data/d_records_butterflies.csv")
d_records_odo <- fread("Data/d_records_dragonflies.csv")
d_records_ortho<- fread("Data/d_records_grasshoppers.csv")

d_records <- d_records_butter %>% 
  mutate(group = "butterflies") %>% 
  bind_rows(d_records_odo %>% 
              mutate(group = "dragonflies"),
            d_records_ortho %>% 
              mutate(group = "grasshoppers"))

# species lists ----------------------------------------------------------------.
v_splist <- fread("Other/specieslist.txt") %>% 
  select(species, group) %>% 
  deframe()

sel_agricultural_species <- readLines("Other/Agricultural_species.txt")

# proportion of distribution in CH ---------------------------------------------.
# Output of R_Proportion_Switzerland.R

d_prop_ch_butter <- readRDS("Data/d_prop_ch_butterflies.rds")
d_prop_ch_ortho <- readRDS("Data/d_prop_ch_grasshoppers.rds")
d_prop_ch_odo <- readRDS("Data/d_prop_ch_dragonflies.rds")

d_prop_ch <- d_prop_ch_butter %>% 
  mutate(group = "butterflies") %>% 
  bind_rows(d_prop_ch_odo %>% 
              mutate(group = "dragonflies"),
            d_prop_ch_ortho %>% 
              mutate(group = "grasshoppers")) %>% 
  mutate(group = factor(group, levels = names(v_col_groups)))

# bigeographic regions and elevation classes -----------------------------------.
# Raw data available from https://data.geo.admin.ch (ch.bafu.biogeographische_regionen, ch.swisstopo.swissalti3d)

d_zone <- fread('Data/d_zone.csv')
d_zone <- d_zone %>% 
  mutate(biogeo5 = substr(zone, 1, 
                          unlist(regexpr("_", zone)) - 1),
         elevation_cat = substr(zone, 
                             unlist(regexpr("_", zone)) + 1, 
                             nchar(zone))) %>% 
  mutate(elevation_cat = factor(elevation_cat, levels = c("low", "high")),
         biogeo5 = factor(biogeo5, levels = c("Jura", "Plateau", 
                                              "NorthernAlps", "CentralAlps",
                                              "SouthernAlps")))

# mean occupancy values --------------------------------------------------------.
# available from https://doi.org/10.16904/envidat.355

d_occ_means <- fread("Data/occ_means.csv")

# climate change variables -----------------------------------------------------.
# Output of R_Climate_change.R

d_climate_lm_smry <- readRDS("Data/d_climate_lm_smry_5.rds")

# rearrange climate data
d_climate_lm_smry <- d_climate_lm_smry %>% 
  gather(mean_sd, value, -c(var, zone, year_start, year_end)) %>% 
  mutate(var = paste(var, mean_sd, sep = "_")) %>% 
  select(-mean_sd) %>% 
  spread(var, value) %>% 
  rename(year_start_model = year_start) %>% 
  mutate(year_start = year_start_model + 5)

# land-use change variables ----------------------------------------------------.
# Output of R_Land_use_change.R

d_as_gam_smry <- readRDS("Data/d_as_gam_smry_5.rds")

# rearrange land use data
d_as_gam_smry <- d_as_gam_smry %>%
  select(-c(mean_mean, mean_sd)) %>% 
  gather(mean_sd, value, -c(var, zone, year_start, year_end)) %>% 
  mutate(var = paste(var, mean_sd, sep = "_"),
         var = gsub("/", "", var)) %>% 
  select(-mean_sd) %>% 
  spread(var, value)

d_drivers <- d_climate_lm_smry %>% 
  full_join(d_as_gam_smry, by = c("zone", "year_end", "year_start")) %>% 
  filter(year_start > 1979)

# set global parameters ########################################################
################################################################################.

v_col_groups <- c(butterflies = "#CF5C36", 
                  grasshoppers = "#EFC88B", 
                  dragonflies = "#06D6A0")

v_labels_groups <- c(butterflies = "Butterflies", 
                     grasshoppers = "Grasshoppers", 
                     dragonflies = "Dragonflies")

d_n_squares <- d_occ_means %>% 
  mutate(group = v_splist[species]) %>% 
  select(group, zone, n_squares) %>% 
  distinct()

# load and summarise occupancy estimates #######################################
################################################################################.

l_trends <- list()
for (sp_i in unique(d_occ_means$species)){
  d_occ_means_target <- d_occ_means %>% 
    filter(species == sp_i)
  
  d_occ_means_target_agg <- d_occ_means_target %>%
    filter(year %in% 1980:2020) %>%
    group_by(year) %>%
    summarise_at(vars(contains("occ_m")),
                 ~sum(n_squares * .)) %>%
    mutate_at(vars(contains("occ_m")),
              ~. / d_n_squares$n_squares[d_n_squares$group == v_splist[sp_i]]) %>% 
    pivot_longer(cols = -year, names_to = "iter", values_to = "occ_m") %>%
    mutate(iter = substr(iter, 8, 11))
  
  l_trends$lm_trends_40 <-
    d_occ_means_target_agg %>%
    group_by(iter) %>%
    summarise(trend = lm(occ_m ~ year)$coefficients[2],
              .groups = "drop") %>%
    mutate(species = sp_i) %>%
    bind_rows(l_trends$lm_trends_40, .)
  
  l_trends$mean_pred_40 <-
    d_occ_means_target_agg %>%
    group_by(iter) %>%
    summarise(z_pred = predict(lm(occ_m ~ year)),
              year = year,
              .groups = "drop") %>%
    mutate(species = sp_i) %>%
    bind_rows(l_trends$mean_pred_40, .)
  
  l_trends$occ_CH40_mean <- d_occ_means_target_agg %>%
    group_by(iter) %>%
    summarise(occ_m = mean(occ_m),
              .groups = "drop") %>%
    mutate(species = sp_i) %>%
    bind_rows(l_trends$occ_CH40_mean, .)
  
  
  l_trends$z_mean_reg <- d_occ_means_target %>%
    filter(year %in% 1980:2020) %>%
    mutate(occ_m = rowMeans(select(., starts_with("occ_m")))) %>%
    select(species, year, zone, n_squares, occ_m) %>%
    bind_rows(l_trends$z_mean_reg, .)
  
  l_trends$lm_trends_5yr_reg <-
    d_occ_means_target %>%
    filter(year %in% 1980:2019) %>%
    select(-n_squares) %>%
    pivot_longer(cols = -c(year, species, zone), names_to = "iter", values_to = "occ_m") %>%
    mutate(iter = substr(iter, 8, 11),
           year_start = floor(year/5) * 5) %>%
    bind_rows(d_occ_means_target %>%
                filter(year %in% seq(1980, 2020, 5)[-1]) %>%
                select(-n_squares) %>%
                pivot_longer(cols = -c(year, species, zone), 
                             names_to = "iter", values_to = "occ_m") %>% 
                mutate(iter = substr(iter, 8, 11),
                       year_start = year - 5)) %>% 
    group_by(species, biogeo_elevation, year_start, iter) %>%
    summarise(trend = lm(occ_m ~ year)$coefficients[2],
              .groups = "drop") %>%
    bind_rows(l_trends$lm_trends_5yr_reg, .)
  
  l_trends$lm_trends_10yr_reg <-
    d_occ_means_target %>% 
    filter(year %in% 1980:2019) %>% 
    select(-n_squares) %>% 
    pivot_longer(cols = -c(year, species, zone), names_to = "iter", values_to = "occ_m") %>% 
    mutate(iter = substr(iter, 8, 11),
           year_start = floor(year/10) * 10) %>%
    bind_rows(d_occ_means_target %>%
                filter(year %in% seq(1980, 2020, 10)[-1]) %>%
                select(-n_squares) %>%
                pivot_longer(cols = -c(year, species, zone), 
                             names_to = "iter", values_to = "occ_m") %>% 
                mutate(iter = substr(iter, 8, 11),
                       year_start = year - 10)) %>% 
    group_by(species, zone, year_start, iter) %>%
    summarise(trend = lm(occ_m ~ year)$coefficients[2],
              .groups = "drop") %>% 
    bind_rows(l_trends$lm_trends_10yr_reg, .)
}

l_trends <- lapply(l_trends, function(x) mutate(x, group = v_splist[species]))

# Overall trend plot ###########################################################
################################################################################.

# function for bootstrapping
f_boot <- function(iter, dat) {
  sel <- sample(seq_len(nrow(dat)), nrow(dat), replace = T)
  which(sort(dat$mean[sel]) > 0)[1] - .5
}

# global mean occupancy (also considering not occupied regions)
d_z_mean_glob <-
  l_trends$z_mean_reg %>% 
  group_by(group, species, A) %>% 
  summarise(z_mean_glob = sum(occ_m * n_squares) / sum(d_n_squares$n_squares[d_n_squares$group %in% group]),
            .groups = "drop") %>% 
  group_by(group, species) %>% 
  summarise(z_mean = mean(z_mean_glob),
            .groups = "drop")

# main panel
set.seed(22)
p1 <-
  l_trends$lm_trends_40 %>% 
  mutate(trend = trend * 40,
         group = factor(group, levels = names(v_col_groups))) %>% 
  group_by(group, species) %>% 
  summarise(mean = mean(trend),
            lower = hdi(trend)$CI_low,
            upper = hdi(trend)$CI_high,
            .groups = "drop") %>% 
  group_by(group) %>% 
  arrange(mean) %>% 
  mutate(index = seq_along(mean),
         median = median(index),
         posneg = which(mean > 0)[1] - .5) %>%
  ungroup() %>% 
  {ggplot(data = .) +
      geom_rect(data = data.frame(group = factor(names(v_col_groups),
                                                 levels = names(v_col_groups))),
                xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
                aes(fill = group),
                alpha = .25) +
      geom_rect(data = . %>% group_by(group) %>%
                  group_map(~ data.frame(t(quantile(sapply(1:1000, f_boot, .),
                                                    probs = c(0.05, 0.95), na.rm = T)))) %>%
                  bind_rows() %>%
                  mutate(group = factor(names(v_col_groups),
                                        levels = names(v_col_groups))),
                aes(xmin = X5., xmax = X95.), ymin = -Inf, ymax = Inf, fill = "#210203", alpha = .2) +
      geom_hline(yintercept = 0, lty = 2) +
      geom_vline(aes(xintercept = median), col = "#210203", lty = 2) +
      geom_vline(aes(xintercept = posneg), col = "#210203") +
      geom_segment(aes(x = index, y = lower, yend = upper, xend = index), color = "grey50") +
      geom_point(aes(x = index, y = mean), size = .75) +
      scale_fill_manual(values = v_col_groups, guide = "none") +
      theme(axis.text.x = element_blank(),
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            plot.margin = margin(0, 0, 0, 0, "cm"),
            panel.spacing.x = unit(.5, "cm"),
            strip.background = element_blank(),
            strip.text = element_blank()) +
      xlab("Species") +
      ylab("Trend estimate") +
      scale_x_reverse() +
      facet_grid(~ group, scales = "free_x", space = "free_x")}

# mean occupancy panel
p2 <-
  l_trends$lm_trends_40 %>% 
  group_by(species) %>% 
  summarise(mean = mean(trend),
            .groups = "drop") %>% 
  left_join(d_z_mean_glob, 
            by = "species") %>% 
  mutate(group = factor(group, levels = names(v_col_groups))) %>% 
  group_by(group) %>% 
  arrange(mean) %>% 
  mutate(index = seq_along(mean)) %>%  
  ggplot() +
  geom_rect(data = data.frame(group = factor(names(v_col_groups),
                                             levels = names(v_col_groups))),
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
            aes(fill = group),
            alpha = .25) +
  geom_col(aes(x = index, y = z_mean)) +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y.right = element_text(angle = 0, vjust = .5),
        panel.border = element_rect(colour = 1, fill = NA),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        panel.spacing.x = unit(.5, "cm"),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  scale_y_continuous(position = "right", breaks = c(0, 1)) +
  scale_fill_manual(values = v_col_groups, guide = "none") +
  coord_cartesian(ylim = c(0, 1)) +
  xlab("Species") +
  ylab("Mean occupancy") +
  scale_x_reverse() +
  facet_grid(~ group, scales = "free_x", space = "free_x")

plot_grid(p1, p2, ncol = 1, rel_heights = c(1, .12),
          align = "v", axis = "lr") 

# effect size inlet ------------------------------------------------------------.

d_gradient <- data.frame(trend = c(-0.4, -0.2, 0, 0.2, 0.4),
                         z_mean = mean(d_z_mean_glob$z_mean)) %>% 
  mutate(z_end = z_mean + trend) 


d_gradient %>% 
  select(trend, z_mean, z_end) %>% 
  gather(var, value, -trend) %>% 
  mutate(x = ifelse(var == "z_mean", 1980, 2019),
         x = ifelse(value < 0, 1980 - 40 * mean(d_z_mean_glob$z_mean) / trend, x),
         value = ifelse(value < 0, 0, value)) %>% 
  ggplot(aes(x = x, y = value)) +
  geom_line(aes(group = trend), size = 1) +
  geom_text(aes(label = formatC(trend, digits = 2, format = "f")), 
            data = function(y) {y[y$var == "z_end", ]},
            hjust = -.1, size = 8 / ggplot2:::.pt) +
  coord_cartesian(ylim = c(0, NA), xlim = c(1980, 2025), clip = "off") +
  annotate(geom = "text", x = 2000.5, y = -.1, label = "Year",
           size = 9 / ggplot2:::.pt) +
  scale_x_continuous(breaks = c(1980, 2019)) +
  scale_y_continuous(breaks = c(0, .25, .5)) +
  theme(axis.title.y = element_text(size = 9),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 8),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA)) +
  xlab("Year") +
  ylab("Mean occupancy")


# Overall trend numbers ########################################################
################################################################################.

# number of species with positive and negative 40-year trends ------------------.

f_n_posneg <- function(x) {
  l_trends$lm_trends_40 %>% 
    mutate(trend = trend * 40) %>% 
    group_by(group, species) %>% 
    summarise(trend = trend[sample(seq_len(length(trend)), 1)],
              .groups = "drop") %>% 
    summarise(pos_all = sum(trend > 0),
              neg_all = sum(trend < 0),
              pos_butter = sum((trend > 0)[group == "butterflies"]),
              neg_butter = sum((trend < 0)[group == "butterflies"]),
              pos_ortho = sum((trend > 0)[group == "grasshoppers"]),
              neg_ortho = sum((trend < 0)[group == "grasshoppers"]),
              pos_odo = sum((trend > 0)[group == "dragonflies"]),
              neg_odo = sum((trend < 0)[group == "dragonflies"]),
              .groups = "drop") %>% 
    pivot_longer(everything(), names_to = c("sign", "group"), names_sep = "_", 
                 values_to = "n") %>% 
    group_by(group) %>% 
    mutate(prop = n / sum(n)) %>% 
    ungroup()
}

set.seed(65)
lapply(1:5000, f_n_posneg) %>% 
  bind_rows() %>% 
  group_by(group, sign) %>% 
  summarise_all(.funs = list(point_estimate = ~mean(.),
                             CI_low = ~hdi(.)$CI_low,
                             CI_high = ~hdi(.)$CI_high))

# mean 40-year trend per quartile ----------------------------------------------.

f_trend_mean <- function(x) {
  l_trends$lm_trends_40 %>% 
    mutate(trend = trend * 40) %>% 
    group_by(species) %>% 
    summarise(trend = trend[sample(seq_len(length(trend)), 1)],
              .groups = "drop") %>% 
    mutate(upper_quartile = quantile(trend, .75),
           lower_quartile = quantile(trend, .25),
           cat = ifelse(trend >= upper_quartile, "upper_quartile",
                        ifelse(trend <= lower_quartile, "lower_quartile",
                               "inbetween"))) %>% 
    group_by(cat) %>% 
    summarise(mean_trend = mean(trend))
}

set.seed(23)
lapply(1:5000, f_trend_mean) %>% 
  bind_rows() %>% 
  group_by(cat) %>% 
  summarise(point_estimate = mean(mean_trend),
            CI_low = hdi(mean_trend)$CI_low,
            CI_high = hdi(mean_trend)$CI_high,
            .groups = "drop")

# proportion distribution change -----------------------------------------------.

d_prop_change <-
  l_trends$mean_pred_40 %>% 
  mutate(z_pred = ifelse(z_pred < 0, 0, z_pred),
         z_pred = ifelse(z_pred > 1, 1, z_pred)) %>%
  group_by(species, iter) %>% 
  summarise(prop_change = (z_pred[A == 2020] - z_pred[A == 1980]) / mean(z_pred),
            .groups = "drop")

f_prop_change <- function(x) {
  l_trends$lm_trends_40 %>% 
    mutate(trend = trend * 40) %>% 
    group_by(species) %>% 
    summarise(iter_sel = iter[sample(seq_len(length(iter)), 1)],
              trend = trend[iter == iter_sel],
              .groups = "drop")%>% 
    mutate(upper_quartile = quantile(trend, .75),
           lower_quartile = quantile(trend, .25),
           cat = ifelse(trend >= upper_quartile, "upper_quartile",
                        ifelse(trend <= lower_quartile, "lower_quartile",
                               "inbetween")),
           iter = iter_sel) %>% 
    left_join(d_prop_change, by = c("species", "iter")) %>% 
    group_by(cat) %>% 
    summarise(mean_prop_change = mean(prop_change))
}

set.seed(7)
lapply(1:5000, f_prop_change) %>% 
  bind_rows() %>% 
  group_by(cat) %>% 
  summarise(point_estimate = mean(mean_prop_change),
            CI_low = hdi(mean_prop_change)$CI_low,
            CI_high = hdi(mean_prop_change)$CI_high,
            .groups = "drop")

# mean occupancy difference between species with positive and negative trends --.

f_sign_occ_m <- function(x) {
  l_trends$lm_trends_40 %>% 
    mutate(trend = trend * 40) %>% 
    group_by(species) %>% 
    summarise(iter_sel = iter[sample(seq_len(length(iter)), 1)],
              trend = trend[iter == iter_sel],
              .groups = "drop") %>% 
    mutate(sign = ifelse(trend > 0, "pos", "neg"),
           iter = iter_sel) %>% 
    left_join(l_trends$occ_CH40_mean, by = c("species", "iter")) %>% 
    group_by(sign) %>% 
    summarise(mean_occ_m = mean(occ_m)) %>% 
    pivot_wider(names_from = "sign", values_from = "mean_occ_m") %>% 
    mutate(diff = pos - neg) %>% 
    pivot_longer(everything(), names_to = "what", "value")
}

set.seed(62)
lapply(1:5000, f_sign_occ_m) %>% 
  bind_rows() %>% 
  group_by(what) %>% 
  summarise(point_estimate = mean(value),
            CI_low = hdi(value)$CI_low,
            CI_high = hdi(value)$CI_high,
            .groups = "drop")

# compare trends to proportion of distribution in Switzerland ##################
################################################################################.

d_chprop <- l_trends$lm_trends_40 %>% 
  mutate(trend = trend * 40,
         group = factor(group, levels = names(v_col_groups))) %>% 
  group_by(group, species) %>% 
  summarise(mean = mean(trend),
            lower = hdi(trend)$CI_low,
            upper = hdi(trend)$CI_high,
            .groups = "drop") %>% 
  left_join(d_prop_ch, by = c("group", "species"))

set.seed(12)
mod_butter <- stan_glm(mean ~ log(prop_ch), data = filter(d_chprop, group == "butterflies"),
                       prior_intercept = normal(location = mean(d_chprop$mean), scale = 5, autoscale = T),
                       prior = normal(location = 0, scale = 5, autoscale = T))
mod_ortho <- stan_glm(mean ~ log(prop_ch), data = filter(d_chprop, group == "grasshoppers"),
                      prior_intercept = normal(location = mean(d_chprop$mean), scale = 5, autoscale = T),
                      prior = normal(location = 0, scale = 5, autoscale = T))
mod_odo <- stan_glm(mean ~ log(prop_ch), data = filter(d_chprop, group == "dragonflies"),
                    prior_intercept = normal(location = mean(d_chprop$mean), scale = 5, autoscale = T),
                    prior = normal(location = 0, scale = 5, autoscale = T))
mod_all <- stan_glm(mean ~ log(prop_ch), data = d_chprop,
                    prior_intercept = normal(location = mean(d_chprop$mean), scale = 5, autoscale = T),
                    prior = normal(location = 0, scale = 5, autoscale = T))


mod_butter <- as.matrix(mod_butter)
mod_ortho <- as.matrix(mod_ortho)
mod_odo <- as.matrix(mod_odo)
mod_all <- as.matrix(mod_all)

m_intercept_butter <- mean(mod_butter[, 1])
m_slope_butter <- mean(mod_butter[, 2])
ymin_m_butter <- m_intercept_butter + m_slope_butter * 
  min(log(d_chprop$prop_ch[d_chprop$group == "butterflies"]))
ymax_m_butter <- m_intercept_butter + m_slope_butter * 
  max(log(d_chprop$prop_ch[d_chprop$group == "butterflies"]))

m_intercept_ortho <- mean(mod_ortho[, 1])
m_slope_ortho <- mean(mod_ortho[, 2])
ymin_m_ortho <- m_intercept_ortho + m_slope_ortho * 
  min(log(d_chprop$prop_ch[d_chprop$group == "grasshoppers"]))
ymax_m_ortho <- m_intercept_ortho + m_slope_ortho * 
  max(log(d_chprop$prop_ch[d_chprop$group == "grasshoppers"]))

m_intercept_odo <- mean(mod_odo[, 1])
m_slope_odo <- mean(mod_odo[, 2])
ymin_m_odo <- m_intercept_odo + m_slope_odo * 
  min(log(d_chprop$prop_ch[d_chprop$group == "dragonflies"]))
ymax_m_odo <- m_intercept_odo + m_slope_odo * 
  max(log(d_chprop$prop_ch[d_chprop$group == "dragonflies"]))

d_mod <- as.data.frame(mod_butter) %>% 
  mutate(xmin = min(log(d_chprop$prop_ch[d_chprop$group == "butterflies"])),
         xmax = max(log(d_chprop$prop_ch[d_chprop$group == "butterflies"])),
         ymin = `(Intercept)` + `log(prop_ch)` * xmin,
         ymax = `(Intercept)` + `log(prop_ch)` * xmax) %>% 
  mutate(group = "butterflies") %>% 
  bind_rows(as.data.frame(mod_ortho) %>% 
              mutate(xmin = min(log(d_chprop$prop_ch[d_chprop$group == "grasshoppers"])),
                     xmax = max(log(d_chprop$prop_ch[d_chprop$group == "grasshoppers"])),
                     ymin = `(Intercept)` + `log(prop_ch)` * xmin,
                     ymax = `(Intercept)` + `log(prop_ch)` * xmax) %>% 
              mutate(group = "grasshoppers")) %>% 
  bind_rows(as.data.frame(mod_odo) %>% 
              mutate(xmin = min(log(d_chprop$prop_ch[d_chprop$group == "dragonflies"])),
                     xmax = max(log(d_chprop$prop_ch[d_chprop$group == "dragonflies"])),
                     ymin = `(Intercept)` + `log(prop_ch)` * xmin,
                     ymax = `(Intercept)` + `log(prop_ch)` * xmax) %>% 
              mutate(group = "dragonflies")) %>% 
  mutate(group = factor(group, levels = names(v_col_groups)))

d_slope_mean <- d_chprop %>% 
  group_by(group) %>% 
  summarise(xmin = min(log(prop_ch)),
            xmax = max(log(prop_ch))) %>% 
  left_join(data.frame(group = c("butterflies", "grasshoppers",
                                 "dragonflies"),
                       ymin = c(ymin_m_butter, ymin_m_ortho, ymin_m_odo),
                       ymax = c(ymax_m_butter, ymax_m_ortho, ymax_m_odo)),
            by = "group") %>% 
  mutate(group = factor(group, levels = names(v_col_groups)))


d_chprop %>% 
  mutate(group = factor(group, levels = names(v_col_groups))) %>% 
  ggplot(aes(x = log(prop_ch), y = mean)) +
  geom_segment(data = d_mod, aes(x = xmin, xend = xmax, y = ymin, yend = ymax, col = group), 
               alpha = .05, size = .1) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(aes(fill = group), pch = 21) +
  geom_segment(data = d_slope_mean,
               aes(x = xmin, xend = xmax,
                   y = ymin, yend = ymax), col = 1) +
  scale_x_continuous(breaks = c(log(0.001), log(0.01), log(0.1), log(1)),
                     labels = c(0.001, 0.01, 0.10, 1.00)) +
  scale_colour_manual(values = v_col_groups) +
  scale_fill_manual(values = v_col_groups) +
  xlab("Proportion of distribution in Switzerland") +
  ylab("Trend estimate") +
  facet_grid(~ group, scales = "free_x",
             labeller = as_labeller(v_labels_groups)) +
  theme(legend.position = "none")

hdi(mod_butter[, 2], .95) %>% str
hdi(mod_ortho[, 2], .95) %>% str
hdi(mod_odo[, 2], .95) %>% str
hdi(mod_all[, 2], .95) %>% str

# Red list status ##############################################################
################################################################################.

d_RL <- fread("Data/RL_status.csv") %>% 
  mutate(RL_status = ifelse(RL_status %in% c("DD", "NE"), "DD/NE", RL_status),
         RL_status = factor(RL_status, levels = c("RE", "CR", "EN", "VU", "NT", "LC", "DD/NE")))

l_trends$lm_trends_40 %>% 
  mutate(trend = trend * 40,
         group = factor(group, levels = names(v_col_groups))) %>% 
  group_by(group, species) %>% 
  summarise(mean = mean(trend),
            lower = hdi(trend)$CI_low,
            upper = hdi(trend)$CI_high,
            .groups = "drop") %>% 
  left_join(d_RL, by = c("group", "species")) %>%
  ggplot(aes(x = RL_status, y = mean)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_jitter(aes(col = group), width = .2) +
  geom_segment(data = function(x) x %>% 
                 group_by(RL_status) %>% 
                 summarise(mean = mean(mean)),
               aes(x = as.numeric(RL_status) - .2 - 1,
                   xend = as.numeric(RL_status) + .2 - 1,
                   yend = mean),
               size = 1) +
  scale_color_manual(values = v_col_groups,
                     labels = v_labels_groups,
                     name = "Group") +
  ylab("Trend estimate") +
  xlab("Red list status")


# STATS ON INTERVAL-WISE TRENDS ---------------- ###############################
################################################################################.

# summarise trend data #########################################################
################################################################################.

d_trends <- l_trends$lm_trends_5yr_reg %>% 
  group_by(group, species, zone, year_start) %>% 
  summarise(trend_mean = mean(trend),
            .groups = "drop")

# full data set for regression model ###########################################
################################################################################.

d_trend_comb <- d_trends %>% 
  mutate(elevation = ifelse(grepl("low$", zone), "low", "high"),
         elevation = factor(elevation, levels = c("low", "high"))) %>% 
  filter(year_start > 1979) %>% 
  mutate(sign_trend = sign(trend_mean),
         trend_mean_sqrt = sqrt(abs(trend_mean)) * sign_trend,
         year_start_f = as.factor(year_start))

scaling_parameters <- list(
  trend_mean = list(mean = mean(d_trend_comb$trend_mean),
                    sd = sd(d_trend_comb$trend_mean)),
  trend_mean_sqrt = list(mean = mean(d_trend_comb$trend_mean_sqrt),
                         sd = sd(d_trend_comb$trend_mean_sqrt)),
  spec_index = list(mean = mean(d_spec$spec_index),
                    sd = sd(d_spec$spec_index)),
  Tniche_mean = list(mean = mean(d_temperature_niche$Tniche_mean),
                   sd = sd(d_temperature_niche$Tniche_mean)),
  T_all_change_mean = list(mean = mean(d_drivers$T_all_change_mean),
                           sd = sd(d_drivers$T_all_change_mean)),
  T_seas_change_mean = list(mean = mean(d_drivers$T_seas_change_mean),
                            sd = sd(d_drivers$T_seas_change_mean)),
  P_warmest_quarter_change_mean = list(mean = mean(d_drivers$P_warmest_quarter_change_mean),
                                       sd = sd(d_drivers$P_warmest_quarter_change_mean)),
  AgrArea_prop_change_mean = list(mean = mean(d_drivers$AgrArea_prop_change_mean),
                             sd = sd(d_drivers$AgrArea_prop_change_mean)),
  LSUGL_change_mean = list(mean = mean(d_drivers$LSUGL_change_mean),
                           sd = sd(d_drivers$LSUGL_change_mean)),
  IAR_change_mean = list(mean = mean(d_drivers$IAR_change_mean),
                            sd = sd(d_drivers$IAR_change_mean))
)



# add predictor and scale data set
d_trend_comb_z <- d_trend_comb %>% 
  mutate(trend_mean = (trend_mean - scaling_parameters$trend_mean$mean) /
           scaling_parameters$trend_mean$sd,
         trend_sd = trend_sd / scaling_parameters$trend_mean$sd,
         trend_mean_sqrt = (trend_mean_sqrt - scaling_parameters$trend_mean_sqrt$mean) /
           scaling_parameters$trend_mean_sqrt$sd) %>% 
  left_join(d_spec %>% 
              mutate(spec_index = (spec_index - scaling_parameters$spec_index$mean) /
                       scaling_parameters$spec_index$sd), by = c("group", "species")) %>% 
  left_join(d_temperature_niche %>% 
              mutate(Tniche_mean = (Tniche_mean - scaling_parameters$Tniche_mean$mean) /
                       scaling_parameters$Tniche_mean$sd), 
            by = c("group", "species")) %>% 
  left_join(d_drivers %>% 
              mutate(T_all_change_mean = (T_all_change_mean - scaling_parameters$T_all_change_mean$mean) /
                       scaling_parameters$T_all_change_mean$sd,
                     T_seas_change_mean = (T_seas_change_mean - scaling_parameters$T_seas_change_mean$mean) /
                       scaling_parameters$T_seas_change_mean$sd,
                     P_warmest_quarter_change_mean = (P_warmest_quarter_change_mean - scaling_parameters$P_warmest_quarter_change_mean$mean) /
                       scaling_parameters$P_warmest_quarter_change_mean$sd,
                     AgrArea_prop_change_mean = (AgrArea_prop_change_mean - scaling_parameters$AgrArea_prop_change_mean$mean) /
                       scaling_parameters$AgrArea_prop_change_mean$sd,
                     LSUGL_change_mean = (LSUGL_change_mean - scaling_parameters$LSUGL_change_mean$mean) /
                       scaling_parameters$LSUGL_change_mean$sd,
                     LSUGL_change_mean = (LSUGL_change_mean - scaling_parameters$LSUGL_change_mean$mean) /
                       scaling_parameters$LSUGL_change_mean$sd,
                     IAR_change_mean = (IAR_change_mean - scaling_parameters$IAR_change_mean$mean) /
                       scaling_parameters$IAR_change_mean$sd), 
            by = c("zone", "year_start")) %>% 
  mutate(species = as.factor(species),
         zone = as.factor(zone),
         group = factor(group, levels = names(v_col_groups)))

# add predictors
d_trend_comb <- d_trend_comb %>% 
  left_join(d_spec, by = c("group", "species")) %>% 
  left_join(d_temperature_niche, by = c("group", "species")) %>%
  left_join(d_drivers, by = c("zone", "year_start")) %>% 
  mutate(species = as.factor(species),
         zone = as.factor(zone),
         group = factor(group, levels = names(v_col_groups)))

# data set only including "agricultural" species for regression model ##########
################################################################################.

d_trend_comb_agrsp <- d_trend_comb %>% 
  filter(species %in% sel_agricultural_species) %>% 
  select(species, zone, year_start, trend_mean, trend_sd, group, elevation, 
         sign_trend, trend_mean_sqrt, year_start_f)

scaling_parameters_agrsp <- list(
  trend_mean = list(mean = mean(d_trend_comb_agrsp$trend_mean),
                    sd = sd(d_trend_comb_agrsp$trend_mean)),
  trend_mean_sqrt = list(mean = mean(d_trend_comb_agrsp$trend_mean_sqrt),
                         sd = sd(d_trend_comb_agrsp$trend_mean_sqrt)),
  spec_index = list(mean = mean(d_spec$spec_index[d_spec$species %in% sel_agricultural_species]),
                    sd = sd(d_spec$spec_index[d_spec$species %in% sel_agricultural_species])),
  Tniche_mean = list(mean = mean(d_temperature_niche$Tniche_mean[d_temperature_niche$species %in% sel_agricultural_species]),
                   sd = sd(d_temperature_niche$Tniche_mean[d_temperature_niche$species %in% sel_agricultural_species])),
  T_all_change_mean = list(mean = mean(d_drivers$T_all_change_mean),
                           sd = sd(d_drivers$T_all_change_mean)),
  T_seas_change_mean = list(mean = mean(d_drivers$T_seas_change_mean),
                            sd = sd(d_drivers$T_seas_change_mean)),
  P_warmest_quarter_change_mean = list(mean = mean(d_drivers$P_warmest_quarter_change_mean),
                                       sd = sd(d_drivers$P_warmest_quarter_change_mean)),
  AgrArea_prop_change_mean = list(mean = mean(d_drivers$AgrArea_prop_change_mean),
                             sd = sd(d_drivers$AgrArea_prop_change_mean)),
  LSUGL_change_mean = list(mean = mean(d_drivers$LSUGL_change_mean),
                           sd = sd(d_drivers$LSUGL_change_mean)),
  IAR_change_mean = list(mean = mean(d_drivers$IAR_change_mean),
                            sd = sd(d_drivers$IAR_change_mean))
)


# add predictor and scale data set
d_trend_comb_agrsp_z <- d_trend_comb_agrsp %>% 
  mutate(trend_mean = (trend_mean - scaling_parameters_agrsp$trend_mean$mean) /
           scaling_parameters_agrsp$trend_mean$sd,
         trend_sd = trend_sd / scaling_parameters_agrsp$trend_mean$sd,
         trend_mean_sqrt = (trend_mean_sqrt - scaling_parameters_agrsp$trend_mean_sqrt$mean) /
           scaling_parameters_agrsp$trend_mean_sqrt$sd) %>% 
  left_join(d_spec %>% 
              filter(species %in% sel_agricultural_species) %>% 
              mutate(spec_index = (spec_index - scaling_parameters_agrsp$spec_index$mean) /
                       scaling_parameters_agrsp$spec_index$sd), by = c("group", "species")) %>% 
  left_join(d_temperature_niche %>% 
              filter(species %in% sel_agricultural_species) %>% 
              mutate(Tniche_mean = (Tniche_mean - scaling_parameters_agrsp$Tniche_mean$mean) /
                       scaling_parameters_agrsp$Tniche_mean$sd), 
            by = c("group", "species")) %>% 
  left_join(d_drivers %>% 
              mutate(T_all_change_mean = (T_all_change_mean - scaling_parameters_agrsp$T_all_change_mean$mean) /
                       scaling_parameters_agrsp$T_all_change_mean$sd,
                     T_seas_change_mean = (T_seas_change_mean - scaling_parameters_agrsp$T_seas_change_mean$mean) /
                       scaling_parameters_agrsp$T_seas_change_mean$sd,
                     P_warmest_quarter_change_mean = (P_warmest_quarter_change_mean - scaling_parameters_agrsp$P_warmest_quarter_change_mean$mean) /
                       scaling_parameters_agrsp$P_warmest_quarter_change_mean$sd,
                     AgrArea_prop_change_mean = (AgrArea_prop_change_mean - scaling_parameters_agrsp$AgrArea_prop_change_mean$mean) /
                       scaling_parameters_agrsp$AgrArea_prop_change_mean$sd,
                     LSUGL_change_mean = (LSUGL_change_mean - scaling_parameters_agrsp$LSUGL_change_mean$mean) /
                       scaling_parameters_agrsp$LSUGL_change_mean$sd,
                     IAR_change_mean = (IAR_change_mean - scaling_parameters_agrsp$IAR_change_mean$mean) /
                       scaling_parameters_agrsp$IAR_change_mean$sd), 
            by = c("zone", "year_start")) %>% 
  mutate(species = as.factor(species),
         zone = as.factor(zone),
         group = factor(group, levels = names(v_col_groups)[1:2]))

# add predictors
d_trend_comb_agrsp <- d_trend_comb_agrsp %>% 
  left_join(d_spec %>% 
              filter(species %in% sel_agricultural_species),
            by = c("group", "species")) %>% 
  left_join(d_temperature_niche %>% 
              filter(species %in% sel_agricultural_species),
            by = c("group", "species")) %>%
  left_join(d_drivers, by = c("zone", "year_start")) %>% 
  mutate(species = as.factor(species),
         zone = as.factor(zone),
         group = factor(group, levels = names(v_col_groups)[1:2]))

# run regression model (version 1) #############################################
################################################################################.

# restricted version (non-"agricultural" species excluded for some parameter estimates)

d_trend_comb_z <- d_trend_comb_z %>% 
  arrange(zone, year_start_f, species) %>% 
  droplevels()

# arrange data
l_data1 <- list(
  N_obs = nrow(d_trend_comb_z),
  N_obs_agrspec = nrow(d_trend_comb_z[d_trend_comb_z$species %in% sel_agricultural_species, ]),
  N_spec = nlevels(d_trend_comb_z$species),
  N_agrspec = n_distinct(d_trend_comb_z$species[d_trend_comb_z$species %in% sel_agricultural_species]),
  N_reg = nlevels(d_trend_comb_z$zone),
  N_int= nlevels(d_trend_comb_z$year_start_f),
  N_group = nlevels(d_trend_comb_z$group),
  
  group = as.numeric(d_trend_comb_z$group),
  spec = as.numeric(d_trend_comb_z$species),
  spec_agrspec = as.numeric(droplevels(d_trend_comb_z$species[d_trend_comb_z$species %in% sel_agricultural_species])),
  reg = as.numeric(d_trend_comb_z$zone),
  
  T_all = d_trend_comb_z$T_all_change_mean,
  T_seas = d_trend_comb_z$T_seas_change_mean,
  P = d_trend_comb_z$P_warmest_quarter_change_mean,
  AgrArea_prop = d_trend_comb_z$AgrArea_prop_change_mean[d_trend_comb_z$species %in% sel_agricultural_species],
  LSUGL = d_trend_comb_z$LSUGL_change_mean[d_trend_comb_z$species %in% sel_agricultural_species],
  IAR = d_trend_comb_z$IAR_change_mean,
  
  Tniche = d_trend_comb_z$Tniche_mean,
  spec_index = d_trend_comb_z$spec_index,
  
  elevation = as.numeric(d_trend_comb_z$elevation == "high"),
  
  year_start = model.matrix(~ year_start_f, data = d_trend_comb_z)[, -1],
  
  index_agrspec = which(d_trend_comb_z$species %in% sel_agricultural_species),
  
  y =  d_trend_comb_z$trend_mean_sqrt,
  
  y_mean = mean(d_trend_comb_z$trend_mean_sqrt)
)



# run model --------------------------------------------------------------------.
set.seed(612)
fit1 <- stan(file = "Stan_Code/Stan_regression_restricted.stan.stan",
     data = l_data1,
     chains = 4, iter = 2000)

# run regression model (version 2) #############################################
################################################################################.

# only 276 "agricultural" species included

d_trend_comb_agrsp_z <- d_trend_comb_agrsp_z %>% 
  arrange(zone, year_start_f, species) %>% 
  droplevels()

# arrange data
l_data2 <- list(
  N_obs = nrow(d_trend_comb_agrsp_z),
  N_spec = nlevels(d_trend_comb_agrsp_z$species),
  N_reg = nlevels(d_trend_comb_agrsp_z$zone),
  N_int= nlevels(d_trend_comb_agrsp_z$year_start_f),
  N_group = nlevels(d_trend_comb_agrsp_z$group),
  
  group = as.numeric(d_trend_comb_agrsp_z$group),
  spec = as.numeric(d_trend_comb_agrsp_z$species),
  reg = as.numeric(d_trend_comb_agrsp_z$zone),
  
  T_all = d_trend_comb_agrsp_z$T_all_change_mean,
  T_seas = d_trend_comb_agrsp_z$T_seas_change_mean,
  P = d_trend_comb_agrsp_z$P_warmest_quarter_change_mean,
  AgrArea_prop = d_trend_comb_agrsp_z$AgrArea_prop_change_mean,
  LSUGL = d_trend_comb_agrsp_z$LSUGL_change_mean,
  IAR = d_trend_comb_agrsp_z$IAR_change_mean,
  
  Tniche = d_trend_comb_agrsp_z$Tniche_mean,
  spec_index = d_trend_comb_agrsp_z$spec_index,
  
  elevation = as.numeric(d_trend_comb_agrsp_z$elevation == "high"),
  
  year_start = model.matrix(~ year_start_f, data = d_trend_comb_agrsp_z)[, -1],
  
  y =  d_trend_comb_agrsp_z$trend_mean_sqrt,
  
  y_mean = mean(d_trend_comb_agrsp_z$trend_mean_sqrt)
)



# run model --------------------------------------------------------------------.
set.seed(877)
fit2 <- stan(file = "Stan_Code/Stan_regression_full.stan.stan",
             data = l_data2,
             chains = 4, iter = 2000)

# run regression model (version 3) #############################################
################################################################################.

# full version with all 390 species

d_trend_comb_z <- d_trend_comb_z %>% 
  arrange(zone, year_start_f, species) %>% 
  droplevels()

# arrange data
l_data3 <- list(
  N_obs = nrow(d_trend_comb_z),
  N_spec = nlevels(d_trend_comb_z$species),
  N_reg = nlevels(d_trend_comb_z$zone),
  N_int= nlevels(d_trend_comb_z$year_start_f),
  N_group = nlevels(d_trend_comb_z$group),
  
  group = as.numeric(d_trend_comb_z$group),
  spec = as.numeric(d_trend_comb_z$species),
  reg = as.numeric(d_trend_comb_z$zone),
  
  T_all = d_trend_comb_z$T_all_change_mean,
  T_seas = d_trend_comb_z$T_seas_change_mean,
  P = d_trend_comb_z$P_warmest_quarter_change_mean,
  AgrArea_prop = d_trend_comb_z$AgrArea_prop_change_mean,
  LSUGL = d_trend_comb_z$LSUGL_change_mean,
  IAR = d_trend_comb_z$IAR_change_mean,
  
  Tniche = d_trend_comb_z$Tniche_mean,
  spec_index = d_trend_comb_z$spec_index,
  
  elevation = as.numeric(d_trend_comb_z$elevation == "high"),
  
  year_start = model.matrix(~ year_start_f, data = d_trend_comb_z)[, -1],
  
  y =  d_trend_comb_z$trend_mean_sqrt,
  
  y_mean = mean(d_trend_comb_z$trend_mean_sqrt)
)



# run model --------------------------------------------------------------------.
set.seed(72)
fit3 <- stan(file = "Stan_Code/Stan_regression_full.stan.stan",
            data = l_data3,
            chains = 4, iter = 2000)


# model diagnostics ############################################################
################################################################################.

# exemplary for model version 1

# extract residuals:
d_y_hat <- extract(fit1, pars = "y_hat")$y_hat
d_resid <- sweep(d_y_hat, 2, l_data1$y)

# QQ plots ---------------------------------------------------------------------.
l_qqplots <- list()
for (i in sample(seq_len(nrow(d_resid)), 16, replace = F)) {
  l_qqplots[[as.character(i)]] <- data.frame(residuals = d_resid[i, ]) %>%
    ggplot(aes(sample = residuals)) +
    stat_qq() +
    stat_qq_line() +
    xlab("Observed") +
    ylab("Predicted") +
    labs(title = paste0("Normal QQ-plot, residuals; rep ", i))
}
plot_grid(plotlist = l_qqplots, nrow = 4)

# Residuals against predictor variables ----------------------------------------.
l_predplots <- list()
for (var_i in c("T_all_change_mean", "T_seas_change_mean", "P_warmest_quarter_change_mean",
                "AgrArea_prop_change_mean", "LSUGL_change_mean", "IAR_change_mean", "Tniche_mean", 
                "spec_index", "year_start")) {
  for (i in sample(seq_len(nrow(d_resid)), 3, replace = F)) {
    p <- d_trend_comb_z %>% 
      rename("var_focal" := sym(var_i)) %>% 
      mutate(resid = d_resid[i, ]) %>% 
      ggplot(aes(x = var_focal, y = resid)) +
      
      geom_hline(yintercept = 0, lty = 2) +
      geom_point() +     
      stat_smooth(method = "loess", method.args = list(family = "symmetric"), col = "red", size = .5,
                  formula = 'y ~ x') +
      xlab(var_i) +
      ylab("Residuals") +
      labs(title = paste0("rep ", i))
    
    l_predplots <- c(l_predplots, list(p))
  }
}
plot_grid(plotlist = l_predplots, nrow = 3, byrow = F)

# ACF plot ---------------------------------------------------------------------.

d_acf <- data.frame()
for (i in seq_len(nrow(d_resid))) {
  acf <- acf_plot(d_resid[i, ], split_by = list(d_trend_comb_z$zone, 
                                                d_trend_comb_z$species),
                  plot = F)
  d_acf <- data.frame(ACF = acf,
                      row.names = NULL) %>% 
    mutate(lag = 0:7,
           run = i) %>% 
    bind_rows(d_acf, .)
}

d_acf %>% 
  group_by(lag) %>% 
  summarise(ACF_mean = mean(ACF),
            ACF_lower = hdi(ACF, ci = .95)$CI_low,
            ACF_upper = hdi(ACF, ci = .95)$CI_high,
            ACF_sd = sd(ACF),
            .groups = "drop") %>% 
  ggplot(aes(x = lag)) +
  geom_hline(yintercept = c(-1, 0, 1), lty = 2) +
  geom_col(aes(y = ACF_mean), fill = "#537D8D") +
  geom_linerange(aes(ymin = ACF_mean - ACF_sd,
                     ymax = ACF_mean + ACF_sd)) +
  scale_x_continuous(breaks = c(0, 2, 4, 6)) +
  coord_cartesian(ylim = c(-1, 1)) +
  xlab("Lag") +
  ylab("ACF mean")


# Tabular model results ########################################################
################################################################################.

# exemplary for model version 1

l_fit <- extract(fit1, pars = c("mu_a", "b_T_all", "b_T_seas", "b_P",  
                               "b_AgrArea_prop", "b_LSUGL", "b_IAR", 
                               "b_Tniche", "b_spec_index", 
                               "b_elevation",
                               "b_T_all_Tniche", "b_T_seas_Tniche", "b_P_Tniche",
                               "b_AgrArea_prop_spec_index", "b_LSUGL_spec_index", "b_IAR_spec_index",
                               "b_T_all_AgrArea_prop", "b_T_seas_AgrArea_prop", "b_P_AgrArea_prop", 
                               "b_T_all_LSUGL", "b_T_seas_LSUGL", "b_P_LSUGL",
                               "b_T_all_IAR", "b_T_seas_IAR", "b_P_IAR",
                               "b_Tniche_elevation", "b_spec_index_elevation", 
                               "b_T_all_elevation", "b_T_seas_elevation", "b_P_elevation", 
                               "b_AgrArea_prop_elevation", "b_LSUGL_elevation",
                               "b_int"))


d_hdi <- data.frame()


# climate main -----------------------------------------------------------------.

for (cl_var_i in c("T_all", "T_seas", "P")) {
  
  hdi_target <- hdi(l_fit[[paste0("b_", cl_var_i)]], ci = c(.8, .95))
  d_hdi <- data.frame(mean = mean(l_fit[[paste0("b_", cl_var_i)]]),
                      lower95 = hdi_target$CI_low[2],
                      lower80 = hdi_target$CI_low[1],
                      upper80 = hdi_target$CI_high[1],
                      upper95 = hdi_target$CI_high[2]) %>% 
    mutate(cl_var = cl_var_i, elevation = "low", category = "climate_main") %>% 
    bind_rows(d_hdi, .)
  
  hdi_target <- hdi(l_fit[[paste0("b_", cl_var_i)]] + l_fit[[paste0("b_", cl_var_i, "_elevation")]], ci = c(.8, .95))
  d_hdi <- data.frame(mean = mean(l_fit[[paste0("b_", cl_var_i)]] + l_fit[[paste0("b_", cl_var_i, "_elevation")]]),
                      lower95 = hdi_target$CI_low[2],
                      lower80 = hdi_target$CI_low[1],
                      upper80 = hdi_target$CI_high[1],
                      upper95 = hdi_target$CI_high[2]) %>% 
    mutate(cl_var = cl_var_i, elevation = "high", category = "climate_main") %>% 
    bind_rows(d_hdi, .)
  
}


# land use main ----------------------------------------------------------------.

for (lu_var_i in c("AgrArea_prop", "LSUGL")) {
  
  hdi_target <- hdi(l_fit[[paste0("b_", lu_var_i)]], ci = c(.8, .95))
  d_hdi <- data.frame(mean = mean(l_fit[[paste0("b_", lu_var_i)]]), 
                      lower95 = hdi_target$CI_low[2],
                      lower80 = hdi_target$CI_low[1],
                      upper80 = hdi_target$CI_high[1],
                      upper95 = hdi_target$CI_high[2]) %>% 
    mutate(lu_var = lu_var_i, elevation = "low", category = "land_use_main") %>% 
    bind_rows(d_hdi, .)
  
  hdi_target <- hdi(l_fit[[paste0("b_", lu_var_i)]] + l_fit[[paste0("b_", lu_var_i, "_elevation")]], ci = c(.8, .95))
  d_hdi <- data.frame(mean = mean(l_fit[[paste0("b_", lu_var_i)]] + l_fit[[paste0("b_", lu_var_i, "_elevation")]]), 
                      lower95 = hdi_target$CI_low[2],
                      lower80 = hdi_target$CI_low[1],
                      upper80 = hdi_target$CI_high[1],
                      upper95 = hdi_target$CI_high[2]) %>% 
    mutate(lu_var = lu_var_i, elevation = "high", category = "land_use_main") %>% 
    bind_rows(d_hdi, .)
  
}

# no altitude interaction for IAR
hdi_target <- hdi(l_fit[["b_IAR"]], ci = c(.8, .95))
d_hdi <- data.frame(mean = mean(l_fit[["b_IAR"]]),
                    lower95 = hdi_target$CI_low[2],
                    lower80 = hdi_target$CI_low[1],
                    upper80 = hdi_target$CI_high[1],
                    upper95 = hdi_target$CI_high[2]) %>% 
  mutate(lu_var = "IAR", elevation = NA, category = "land_use_main") %>% 
  bind_rows(d_hdi, .)


# climate x land use -----------------------------------------------------------.

for (cl_var_i in c("T_all", "T_seas", "P")) {
  for (lu_var_i in c("AgrArea_prop", "LSUGL", "IAR")) {
    
    hdi_target <- hdi(l_fit[[paste0("b_", cl_var_i, "_", lu_var_i)]], ci = c(.8, .95))
    d_hdi <- data.frame(mean = mean(l_fit[[paste0("b_", cl_var_i, "_", lu_var_i)]]), 
                        lower95 = hdi_target$CI_low[2],
                        lower80 = hdi_target$CI_low[1],
                        upper80 = hdi_target$CI_high[1],
                        upper95 = hdi_target$CI_high[2]) %>%
      mutate(cl_var = cl_var_i, lu_var = lu_var_i, category = "climate_land_use") %>% 
      bind_rows(d_hdi, .)
    
  }
}

# T niche main -----------------------------------------------------------------.

hdi_target <- hdi(l_fit$b_Tniche, ci = c(.8, .95))
d_hdi <- data.frame(mean = mean(l_fit$b_Tniche), 
                    lower95 = hdi_target$CI_low[2],
                    lower80 = hdi_target$CI_low[1],
                    upper80 = hdi_target$CI_high[1],
                    upper95 = hdi_target$CI_high[2]) %>%
  mutate(trT_var = "Tniche",  elevation = "low", category = "Tniche_main") %>% 
  bind_rows(d_hdi, .)

hdi_target <- hdi(l_fit$b_Tniche + l_fit$b_Tniche_elevation, ci = c(.8, .95))
d_hdi <- data.frame(mean = mean(l_fit$b_Tniche + l_fit$b_Tniche_elevation), 
                    lower95 = hdi_target$CI_low[2],
                    lower80 = hdi_target$CI_low[1],
                    upper80 = hdi_target$CI_high[1],
                    upper95 = hdi_target$CI_high[2]) %>%
  mutate(trT_var = "Tniche",  elevation = "high", category = "Tniche_main") %>% 
  bind_rows(d_hdi, .)


# T niche x climate ------------------------------------------------------------.

for (cl_var_i in c("T_all", "T_seas", "P")) {
  
  hdi_target <- hdi(l_fit[[paste0("b_", cl_var_i, "_Tniche")]], ci = c(.8, .95))
  d_hdi <- data.frame(mean = mean(l_fit[[paste0("b_", cl_var_i, "_Tniche")]]), 
                      lower95 = hdi_target$CI_low[2],
                      lower80 = hdi_target$CI_low[1],
                      upper80 = hdi_target$CI_high[1],
                      upper95 = hdi_target$CI_high[2]) %>%
    mutate(cl_var = cl_var_i, trT_var = "Tniche", category = "climate_Tniche") %>% 
    bind_rows(d_hdi, .)
  
}

# specialisation ---------------------------------------------------------------.

hdi_target <- hdi(l_fit$b_spec_index, ci = c(.8, .95))
d_hdi <- data.frame(mean = mean(l_fit$b_spec_index), 
                    lower95 = hdi_target$CI_low[2],
                    lower80 = hdi_target$CI_low[1],
                    upper80 = hdi_target$CI_high[1],
                    upper95 = hdi_target$CI_high[2]) %>%
  mutate(spec_var = "spec_index", elevation = "low", category = "spec_main") %>% 
  bind_rows(d_hdi, .)

hdi_target <- hdi(l_fit$b_spec_index + l_fit$b_spec_index_elevation, ci = c(.8, .95))
d_hdi <- data.frame(mean = mean(l_fit$b_spec_index + l_fit$b_spec_index_elevation), 
                    lower95 = hdi_target$CI_low[2],
                    lower80 = hdi_target$CI_low[1],
                    upper80 = hdi_target$CI_high[1],
                    upper95 = hdi_target$CI_high[2]) %>%
  mutate(spec_var = "spec_index", elevation = "high", category = "spec_main") %>% 
  bind_rows(d_hdi, .)

# land use x specialisation ----------------------------------------------------.

for (lu_var_i in c("AgrArea_prop", "LSUGL", "IAR")) {
  
  hdi_target <- hdi(l_fit[[paste0("b_", lu_var_i, "_spec_index")]], ci = c(.8, .95))
  d_hdi <- data.frame(mean = mean(l_fit[[paste0("b_", lu_var_i, "_spec_index")]]), 
                      lower95 = hdi_target$CI_low[2],
                      lower80 = hdi_target$CI_low[1],
                      upper80 = hdi_target$CI_high[1],
                      upper95 = hdi_target$CI_high[2]) %>%
    mutate(lu_var = lu_var_i, spec_var = "spec_index", category = "land_use_spec") %>% 
    bind_rows(d_hdi, .)
}

# elevation --------------------------------------------------------------------.

hdi_target <- hdi(l_fit$b_elevation, ci = c(.8, .95))
d_hdi <- data.frame(mean = mean(l_fit$b_elevation), 
                    lower95 = hdi_target$CI_low[2],
                    lower80 = hdi_target$CI_low[1],
                    upper80 = hdi_target$CI_high[1],
                    upper95 = hdi_target$CI_high[2]) %>%
  mutate(category = "elevation_main") %>% 
  bind_rows(d_hdi, .)

# group intercepts -------------------------------------------------------------.

for (gr_i in 1:nlevels(d_trend_comb_z$group)) {
  hdi_target <- hdi(l_fit$mu_a[, gr_i], ci = c(.8, .95))
  d_hdi <- data.frame(mean = mean(l_fit$mu_a[, gr_i]),
                      lower95 = hdi_target$CI_low[2],
                      lower80 = hdi_target$CI_low[1],
                      upper80 = hdi_target$CI_high[1],
                      upper95 = hdi_target$CI_high[2]) %>% 
    mutate(category = "group_intercept", group = paste0("Intercept (", 
                                                        v_labels_groups[levels(d_trend_comb_z$group)[gr_i]],
                                                        ")")) %>% 
    bind_rows(d_hdi, .)
}

# time interval ----------------------------------------------------------------.

v_intervals <- c("1985-1990", "1990-1995",
                 "1995-2000", "2000-2005", 
                 "2005-2010", "2010-2015",
                 "2015-2020")

for (int_i in 1:7) {
  hdi_target <- hdi(l_fit$b_int[, int_i], ci = c(.8, .95))
  d_hdi <- data.frame(mean = mean(l_fit$b_int[, int_i]),
                      lower95 = hdi_target$CI_low[2],
                      lower80 = hdi_target$CI_low[1],
                      upper80 = hdi_target$CI_high[1],
                      upper95 = hdi_target$CI_high[2]) %>% 
    mutate(category = "time_interval", int = paste("Interval", v_intervals[int_i])) %>% 
    bind_rows(d_hdi, .)
}

# final edits ------------------------------------------------------------------.


d_hdi <- d_hdi %>% 
  mutate(cl_var = factor(cl_var,
                         levels = c("T_all",
                                    "T_seas",
                                    "P"),
                         labels = c("T. mean",
                                    "T. seasonality",
                                    "P. summer")),
         lu_var = factor(lu_var,
                         levels = c("AgrArea_prop",
                                    "LSUGL",
                                    "IAR"),
                         labels = c("Agr. area",
                                    "Grassland-use int.",
                                    "Crop-use int.")),
         trT_var = factor(trT_var,
                          levels = c("Tniche"),
                          labels = c("Temp. niche")),
         spec_var = factor(spec_var,
                           levels = c("spec_index"),
                           labels = c("Specialisation")),
         int = factor(int),
         group = factor(group),
         mean = mean * scaling_parameters$trend_mean_sqrt$sd,
         lower95 = lower95 * scaling_parameters$trend_mean_sqrt$sd,
         lower80 = lower80 * scaling_parameters$trend_mean_sqrt$sd,
         upper80 = upper80 * scaling_parameters$trend_mean_sqrt$sd,
         upper95 = upper95 * scaling_parameters$trend_mean_sqrt$sd) %>% # change to original level
  rowwise() %>% 
  mutate(var = ifelse(is.na(elevation), paste(c(cl_var, lu_var, trT_var, spec_var, int, group), collapse = " × "),
                      paste0(paste(c(cl_var, lu_var, trT_var, spec_var, int, group), collapse = " × "), " (", elevation, ")")),
         var = gsub(" × NA", "", var),
         var = gsub("NA × ", "", var),
         var= ifelse(category == "elevation_main", "Elevation (high)", var)) %>% 
  ungroup()

d_hdi %>%
  select(var, lower95, lower80, mean, upper80, upper95)


# model coefficients plot ######################################################
################################################################################.

# exemplary for model version 1

l_fit <- extract(fit1, pars = c("b_T_all", "b_T_seas", "b_P",  
                               "b_AgrArea_prop", "b_LSUGL", "b_IAR", 
                               "b_Tniche", "b_spec_index", 
                               "b_elevation",
                               "b_T_all_Tniche", "b_T_seas_Tniche", "b_P_Tniche",
                               "b_AgrArea_prop_spec_index", "b_LSUGL_spec_index", "b_IAR_spec_index",
                               "b_T_all_AgrArea_prop", "b_T_seas_AgrArea_prop", "b_P_AgrArea_prop", 
                               "b_T_all_LSUGL", "b_T_seas_LSUGL", "b_P_LSUGL",
                               "b_T_all_IAR", "b_T_seas_IAR", "b_P_IAR",
                               "b_Tniche_elevation", "b_spec_index_elevation", 
                               "b_T_all_elevation", "b_T_seas_elevation", "b_P_elevation", 
                               "b_AgrArea_prop_elevation", "b_LSUGL_elevation"))

d_density <- data.frame()
d_mean <- data.frame()

# threshold_y <- .01 # exclude low densities

# climate main -----------------------------------------------------------------.

for (cl_var_i in c("T_all", "T_seas", "P")) {
  
  hdi_target <- hdi(l_fit[[paste0("b_", cl_var_i)]], ci = .95)
  d_density <- density(l_fit[[paste0("b_", cl_var_i)]], 
                       from = hdi_target$CI_low, to = hdi_target$CI_high) %>% 
    as.data.frame() %>% 
    mutate(cl_var = cl_var_i, elevation = "low", category = "climate_main",
           xwidth = mean(diff(x))) %>% 
    bind_rows(d_density, .)
  d_mean <- data.frame(mean = mean(l_fit[[paste0("b_", cl_var_i)]])) %>% 
    mutate(cl_var = cl_var_i, elevation = "low", category = "climate_main") %>% 
    bind_rows(d_mean, .)
  
  hdi_target <- hdi(l_fit[[paste0("b_", cl_var_i)]] + l_fit[[paste0("b_", cl_var_i, "_elevation")]], ci = .95)
  d_density <- density(l_fit[[paste0("b_", cl_var_i)]] + l_fit[[paste0("b_", cl_var_i, "_elevation")]], 
                       from = hdi_target$CI_low, to = hdi_target$CI_high) %>% 
    as.data.frame() %>% 
    mutate(cl_var = cl_var_i, elevation = "high", category = "climate_main",
           xwidth = mean(diff(x))) %>% 
    bind_rows(d_density, .)
  d_mean <- data.frame(mean = mean(l_fit[[paste0("b_", cl_var_i)]] + l_fit[[paste0("b_", cl_var_i, "_elevation")]])) %>% 
    mutate(cl_var = cl_var_i, elevation = "high", category = "climate_main") %>% 
    bind_rows(d_mean, .)
  
}

# land use main ----------------------------------------------------------------.

for (lu_var_i in c("AgrArea_prop", "LSUGL")) {
  
  hdi_target <- hdi(l_fit[[paste0("b_", lu_var_i)]], ci = .95)
  d_density <- density(l_fit[[paste0("b_", lu_var_i)]], 
                       from = hdi_target$CI_low, to = hdi_target$CI_high) %>% 
    as.data.frame() %>% 
    mutate(lu_var = lu_var_i, elevation = "low", category = "land_use_main",
           xwidth = mean(diff(x))) %>% 
    bind_rows(d_density, .)
  d_mean <- data.frame(mean = mean(l_fit[[paste0("b_", lu_var_i)]])) %>% 
    mutate(lu_var = lu_var_i, elevation = "low", category = "land_use_main") %>%
    bind_rows(d_mean, .)
  
  hdi_target <- hdi(l_fit[[paste0("b_", lu_var_i)]] + l_fit[[paste0("b_", lu_var_i, "_elevation")]], ci = .95)
  d_density <- density(l_fit[[paste0("b_", lu_var_i)]] + l_fit[[paste0("b_", lu_var_i, "_elevation")]], 
                       from = hdi_target$CI_low, to = hdi_target$CI_high) %>% 
    as.data.frame() %>% 
    mutate(lu_var = lu_var_i, elevation = "high", category = "land_use_main",
           xwidth = mean(diff(x))) %>% 
    bind_rows(d_density, .)
  d_mean <- data.frame(mean = mean(l_fit[[paste0("b_", lu_var_i)]] + l_fit[[paste0("b_", lu_var_i, "_elevation")]])) %>% 
    mutate(lu_var = lu_var_i, elevation = "high", category = "land_use_main") %>%
    bind_rows(d_mean, .)
  
}

# no altitude interaction for IAR
hdi_target <- hdi(l_fit[["b_IAR"]], ci = .95)
d_density <- density(l_fit[["b_IAR"]], 
                     from = hdi_target$CI_low, to = hdi_target$CI_high) %>% 
  as.data.frame() %>% 
  mutate(lu_var = "IAR", elevation = NA, category = "land_use_main",
         xwidth = mean(diff(x))) %>% 
  bind_rows(d_density, .)
d_mean <- data.frame(mean = mean(l_fit[["b_IAR"]])) %>% 
  mutate(lu_var = "IAR", elevation = NA, category = "land_use_main") %>%
  bind_rows(d_mean, .)


# climate x land use -----------------------------------------------------------.

for (cl_var_i in c("T_all", "T_seas", "P")) {
  for (lu_var_i in c("AgrArea_prop", "LSUGL", "IAR")) {
    
    hdi_target <- hdi(l_fit[[paste0("b_", cl_var_i, "_", lu_var_i)]], ci = .95)
    d_density <- density(l_fit[[paste0("b_", cl_var_i, "_", lu_var_i)]], 
                         from = hdi_target$CI_low, to = hdi_target$CI_high) %>% 
      as.data.frame() %>% 
      mutate(cl_var = cl_var_i, lu_var = lu_var_i, category = "climate_land_use",
             xwidth = mean(diff(x))) %>% 
      bind_rows(d_density, .)
    d_mean <- data.frame(mean = mean(l_fit[[paste0("b_", cl_var_i, "_", lu_var_i)]])) %>% 
      mutate(cl_var = cl_var_i, lu_var = lu_var_i, category = "climate_land_use") %>%
      bind_rows(d_mean, .)
    
    
  }
}

# T niche main -----------------------------------------------------------------.

hdi_target <- hdi(l_fit$b_Tniche, ci = .95)
d_density <- density(l_fit$b_Tniche, 
                     from = hdi_target$CI_low, to = hdi_target$CI_high) %>% 
  as.data.frame() %>% 
  mutate(trT_var = "Tniche",  elevation = "low", category = "Tniche_main",
         xwidth = mean(diff(x))) %>% 
  bind_rows(d_density, .)
d_mean <- data.frame(mean = mean(l_fit$b_Tniche)) %>% 
  mutate(trT_var = "Tniche",  elevation = "low", category = "Tniche_main") %>%
  bind_rows(d_mean, .)

hdi_target <- hdi(l_fit$b_Tniche + l_fit$b_Tniche_elevation, ci = .95)
d_density <- density(l_fit$b_Tniche + l_fit$b_Tniche_elevation, 
                     from = hdi_target$CI_low, to = hdi_target$CI_high) %>% 
  as.data.frame() %>% 
  mutate(trT_var = "Tniche",  elevation = "high", category = "Tniche_main",
         xwidth = mean(diff(x))) %>% 
  bind_rows(d_density, .)
d_mean <- data.frame(mean = mean(l_fit$b_Tniche + l_fit$b_Tniche_elevation)) %>% 
  mutate(trT_var = "Tniche",  elevation = "high", category = "Tniche_main") %>%
  bind_rows(d_mean, .)


# T niche x climate ------------------------------------------------------------.

for (cl_var_i in c("T_all", "T_seas", "P")) {
  
  hdi_target <- hdi(l_fit[[paste0("b_", cl_var_i, "_Tniche")]], ci = .95)
  d_density <- density(l_fit[[paste0("b_", cl_var_i, "_Tniche")]], 
                       from = hdi_target$CI_low, to = hdi_target$CI_high) %>% 
    as.data.frame() %>% 
    mutate(cl_var = cl_var_i, trT_var = "Tniche", category = "climate_Tniche",
           xwidth = mean(diff(x))) %>% 
    bind_rows(d_density, .)
  d_mean <- data.frame(mean = mean(l_fit[[paste0("b_", cl_var_i, "_Tniche")]])) %>% 
    mutate(cl_var = cl_var_i, trT_var = "Tniche", category = "climate_Tniche") %>%
    bind_rows(d_mean, .)
  
}

# specialisation ---------------------------------------------------------------.

hdi_target <- hdi(l_fit$b_spec_index, ci = .95)
d_density <- density(l_fit$b_spec_index, 
                     from = hdi_target$CI_low, to = hdi_target$CI_high) %>% 
  as.data.frame() %>% 
  mutate(spec_var = "spec_index", elevation = "low", category = "spec_main",
         xwidth = mean(diff(x))) %>% 
  bind_rows(d_density, .)
d_mean <- data.frame(mean = mean(l_fit$b_spec_index)) %>% 
  mutate(spec_var = "spec_index", elevation = "low", category = "spec_main") %>%
  bind_rows(d_mean, .)

hdi_target <- hdi(l_fit$b_spec_index + l_fit$b_spec_index_elevation, ci = .95)
d_density <- density(l_fit$b_spec_index + l_fit$b_spec_index_elevation, 
                     from = hdi_target$CI_low, to = hdi_target$CI_high) %>% 
  as.data.frame() %>% 
  mutate(spec_var = "spec_index", elevation = "high", category = "spec_main",
         xwidth = mean(diff(x))) %>% 
  bind_rows(d_density, .)
d_mean <- data.frame(mean = mean(l_fit$b_spec_index + l_fit$b_spec_index_elevation)) %>% 
  mutate(spec_var = "spec_index", elevation = "high", category = "spec_main") %>%
  bind_rows(d_mean, .)

# land use x specialisation ----------------------------------------------------.

for (lu_var_i in c("AgrArea_prop", "LSUGL", "IAR")) {
  
  hdi_target <- hdi(l_fit[[paste0("b_", lu_var_i, "_spec_index")]], ci = .95)
  d_density <- density(l_fit[[paste0("b_", lu_var_i, "_spec_index")]], 
                       from = hdi_target$CI_low, to = hdi_target$CI_high) %>% 
    as.data.frame() %>% 
    mutate(lu_var = lu_var_i, spec_var = "spec_index", category = "land_use_spec",
           xwidth = mean(diff(x))) %>% 
    bind_rows(d_density, .)
  d_mean <- data.frame(mean = mean(l_fit[[paste0("b_", lu_var_i, "_spec_index")]])) %>% 
    mutate(lu_var = lu_var_i, spec_var = "spec_index", category = "land_use_spec") %>%
    bind_rows(d_mean, .)
}

# elevation --------------------------------------------------------------------.

hdi_target <- hdi(l_fit$b_elevation, ci = .95)
d_density <- density(l_fit$b_elevation, 
                     from = hdi_target$CI_low, to = hdi_target$CI_high) %>% 
  as.data.frame() %>% 
  mutate(category = "elevation_main",
         xwidth = mean(diff(x))) %>% 
  bind_rows(d_density, .)
d_mean <- data.frame(mean = mean(l_fit$b_elevation)) %>% 
  mutate(category = "elevation_main") %>%
  bind_rows(d_mean, .)


# final edits ------------------------------------------------------------------.

d_density <- d_density %>% 
  mutate(cl_var = factor(cl_var, 
                         levels = c("T_all",
                                    "T_seas",
                                    "P")),
         lu_var = factor(lu_var, 
                         levels = c("AgrArea_prop", 
                                    "LSUGL",
                                    "IAR")),
         trT_var = factor(trT_var, 
                          levels = c("Tniche")),
         spec_var = factor(spec_var,
                           levels = c("spec_index")),
         x = x * scaling_parameters$trend_mean_sqrt$sd,
         xwidth = xwidth * scaling_parameters$trend_mean_sqrt$sd) # change to original level

d_mean <- d_mean %>% 
  mutate(cl_var = factor(cl_var, 
                         levels = c("T_all",
                                    "T_seas",
                                    "P")),
         lu_var = factor(lu_var, 
                         levels = c("AgrArea_prop", 
                                    "LSUGL",
                                    "IAR")),
         trT_var = factor(trT_var, 
                          levels = c("Tniche")),
         spec_var = factor(spec_var,
                           levels = c("spec_index")),
         mean = mean * scaling_parameters$trend_mean_sqrt$sd)

# plotting squares -------------------------------------------------------------.

scale_height <- 1

range_effect <- c(min(d_density$x),
                  max(d_density$x))
range_effect <- c(-max(abs(range_effect)), # to make it symmetric
                  max(abs(range_effect)))

p1 <-
  d_density %>%
  left_join(d_mean,
            c("cl_var", "elevation", "category", "lu_var", "trT_var", "spec_var")) %>% 
  filter(category == "climate_main") %>%
  mutate(sign_high = ifelse(elevation == "high", 1, -1)) %>% 
  group_by(cl_var, elevation) %>% 
  mutate(x0 = ifelse((any(x > 0) & any(x < 0) & y == max(y)), 0, NA),
         label = formatC(mean, format = "f", digits = 3)) %>% 
  ungroup() %>% 
  mutate(x_cat = cut(x, breaks = c(-Inf, -.015, -.01, -.005, 0, .005, .01, .015, Inf))) %>% 
  ggplot() +
  geom_line(aes(x = x, y = y * sign_high), col = NA) +
  geom_tile(aes(x = x, fill = x_cat, y = y/2 * sign_high * scale_height, height = y * scale_height, width = xwidth)) +
  geom_line(aes(x = x, y = y * sign_high * scale_height)) +
  geom_segment(aes(x = x0, xend = x0, y = 0, yend = y * sign_high * scale_height), lty = 2, size = .5) +
  geom_text(data=function(x) {x %>% group_by(cl_var, elevation) %>% 
      filter(abs(x - median(x)) == min(abs(x - median(x)))) %>% 
      select(-x) %>% distinct()}, aes(x = mean, label = label, y = y/4 * sign_high * scale_height), col = 1) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(8, "RdYlBu"), drop = F) +
  theme_nothing() +
  theme(panel.spacing.y = unit(-.2, "cm"),
        panel.background = element_rect(fill = "grey90", colour = NA)) +
  facet_wrap(elevation ~ cl_var, scales = "free", nrow = 2)


p2a <-
  d_density %>% 
  left_join(d_mean,
            c("cl_var", "elevation", "category", "lu_var", "trT_var", "spec_var")) %>% 
  filter(category == "land_use_main",
         lu_var != "IAR") %>% 
  mutate(sign_high = ifelse(elevation == "high", 1, -1)) %>% 
  group_by(lu_var, elevation) %>% 
  mutate(x0 = ifelse((any(x > 0) & any(x < 0) & y == max(y)), 0, NA),
         label = formatC(mean, format = "f", digits = 3)) %>% 
  ungroup() %>% 
  mutate(x_cat = cut(x, breaks = c(-Inf, -.015, -.01, -.005, 0, .005, .01, .015, Inf))) %>% 
  ggplot() +
  geom_line(aes(x = x, y = y * sign_high), col = NA) +
  geom_tile(aes(x = x, fill = x_cat, y = y/2 * sign_high * scale_height, height = y * scale_height, width = xwidth)) +
  geom_line(aes(x = x, y = y * sign_high * scale_height)) +
  geom_segment(aes(x = x0, xend = x0, y = 0, yend = y * sign_high * scale_height), lty = 2, size = .5) +
  geom_text(data=function(x) {x %>% group_by(lu_var, elevation) %>% 
      filter(abs(x - median(x)) == min(abs(x - median(x)))) %>% 
      select(-x) %>% distinct()}, aes(x = mean, label = label, y = y/4 * sign_high * scale_height), col = 1) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(8, "RdYlBu"), drop = F) +
  theme_nothing() +
  theme(panel.spacing.y = unit(-.2, "cm"),
        panel.background = element_rect(fill = "grey90", colour = NA)) +
  facet_wrap(elevation ~ lu_var, scales = "free", nrow = 2)

p2b <- 
  d_density %>% 
  left_join(d_mean,
            c("cl_var", "elevation", "category", "lu_var", "trT_var", "spec_var")) %>% 
  filter(category == "land_use_main",
         lu_var == "IAR") %>% 
  mutate(x0 = ifelse((any(x > 0) & any(x < 0) & y == max(y)), 0, NA),
         x0 = as.numeric(x0),
         label = formatC(mean, format = "f", digits = 3)) %>% 
  mutate(x_cat = cut(x, breaks = c(-Inf,  -.015, -.01, -.005, 0, .005, .01, .015, Inf))) %>% 
  ggplot() +
  geom_tile(aes(x = x, fill = x_cat, y = y/2, height = y, width = xwidth)) +
  geom_line(aes(x = x, y = y)) +
  geom_segment(aes(x = x0, xend = x0, y = 0, yend = y), lty = 2, size = .5) +
  geom_text(data=function(x) {x %>%
      filter(abs(x - median(x)) == min(abs(x - median(x)))) %>% 
      select(-x) %>% distinct()}, aes(x = mean, label = label, y = y/4 * scale_height), col = 1) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(8, "RdYlBu"), drop = F) +  
  theme_nothing() +
  theme(panel.background = element_rect(fill = "grey90", colour = NA))


p3 <-
  d_density %>% 
  left_join(d_mean,
            c("cl_var", "elevation", "category", "lu_var", "trT_var", "spec_var")) %>% 
  filter(category == "Tniche_main") %>% 
  mutate(sign_high = ifelse(elevation == "high", 1, -1)) %>% 
  group_by(trT_var, elevation) %>% 
  mutate(x0 = ifelse((any(x > 0) & any(x < 0) & y == max(y)), 0, NA),
         x0 = as.numeric(x0),
         label = formatC(mean, format = "f", digits = 3)) %>% 
  ungroup() %>% 
  mutate(x_cat = cut(x, breaks = c(-Inf, -.015, -.01, -.005, 0, .005, .01, .015,  Inf))) %>% 
  ggplot() +
  geom_tile(aes(x = x, fill = x_cat, y = y/2 * sign_high, height = y, width = xwidth)) +
  geom_line(aes(x = x, y = y * sign_high)) +
  geom_segment(aes(x = x0, xend = x0, y = 0, yend = y * sign_high), lty = 2, size = .5) +
  geom_text(data=function(x) {x %>% group_by(trT_var, elevation) %>% 
      filter(abs(x - median(x)) == min(abs(x - median(x)))) %>% 
      select(-x) %>% distinct()}, aes(x = mean, label = label, y = y/4 * sign_high * scale_height), col = 1) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(8, "RdYlBu"), drop = F) +
  theme_nothing() +
  theme(panel.spacing.y = unit(-.1, "cm"),
        panel.background = element_rect(fill = "grey90", colour = NA)) +
  facet_wrap(elevation ~ trT_var, scales = "free", nrow = 2)

p4 <-
  d_density %>% 
  left_join(d_mean,
            c("cl_var", "elevation", "category", "lu_var", "trT_var", "spec_var")) %>% 
  filter(category == "spec_main") %>% 
  mutate(sign_high = ifelse(elevation == "high", 1, -1)) %>% 
  group_by(elevation) %>% 
  mutate(x_double = median(x) + (x - median(x)) * 2) %>% 
  group_by(spec_var, elevation) %>% 
  mutate(x0 = ifelse((any(x > 0) & any(x < 0) & y == max(y)), 0, NA),
         x0 = as.numeric(x0),
         label = formatC(mean, format = "f", digits = 3)) %>% 
  ungroup() %>%
  mutate(x_cat = cut(x, breaks = c(-Inf, -.015, -.01, -.005, 0, .005, .01, .015, Inf))) %>% 
  ggplot() +
  geom_tile(aes(x = x, fill = x_cat, y = y/2 * sign_high, height = y, width = xwidth)) +
  geom_line(aes(x = x, y = y * sign_high)) +
  geom_segment(aes(x = x0, xend = x0, y = 0, yend = y * sign_high), lty = 2, size = .5) +
  geom_text(data=function(x) {x %>% group_by(spec_var, elevation) %>% 
      filter(abs(x - median(x)) == min(abs(x - median(x)))) %>% 
      select(-x) %>% distinct()}, aes(x = mean, label = label, y = y/4 * sign_high * scale_height), col = 1) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(8, "RdYlBu"), drop = F) +
  theme_nothing() +
  theme(panel.spacing.y = unit(-.1, "cm"),
        panel.background = element_rect(fill = "grey90", colour = NA)) +
  facet_wrap(elevation ~ spec_var, scales = "free", nrow = 2)

p5 <-
  d_density %>% 
  left_join(d_mean,
            c("cl_var", "elevation", "category", "lu_var", "trT_var", "spec_var")) %>% 
  filter(category == "climate_land_use") %>% 
  group_by(lu_var, cl_var) %>% 
  mutate(x0 = ifelse((any(x > 0) & any(x < 0) & y == max(y)), 0, NA),
         x0 = as.numeric(x0),
         label = formatC(mean, format = "f", digits = 3)) %>% 
  ungroup() %>%
  mutate(x_cat = cut(x, breaks = c(-Inf, -.015, -.01, -.005, 0, .005, .01, .015, Inf))) %>% 
  ggplot() +
  geom_tile(aes(x = x, fill = x_cat, y = y/2, height = y, width = xwidth)) +
  geom_line(aes(x = x, y = y)) +
  geom_segment(aes(x = x0, xend = x0, y = 0, yend = y), lty = 2, size = .5) +
  geom_text(data=function(x) {x %>% group_by(lu_var, cl_var) %>% 
      filter(abs(x - median(x)) == min(abs(x - median(x)))) %>% 
      select(-x) %>% distinct()}, aes(x = mean, label = label, y = y/4 * scale_height), col = 1) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(8, "RdYlBu"), drop = F) +
  theme_nothing() +
  theme(panel.background = element_rect(fill = "grey90", colour = NA)) +
  facet_wrap(lu_var ~ cl_var, scales = "free", nrow = 3)

p6 <-
  d_density %>% 
  left_join(d_mean,
            c("cl_var", "elevation", "category", "lu_var", "trT_var", "spec_var")) %>% 
  filter(category == "climate_Tniche") %>% 
  group_by(trT_var, cl_var) %>% 
  mutate(x0 = ifelse((any(x > 0) & any(x < 0) & y == max(y)), 0, NA),
         x0 = as.numeric(x0),
         label = formatC(mean, format = "f", digits = 3)) %>% 
  ungroup() %>% 
  mutate(x_cat = cut(x, breaks = c(-Inf, -.015, -.01, -.005, 0, .005, .01, .015, Inf))) %>% 
  ggplot() +
  geom_tile(aes(x = x, fill = x_cat, y = y/2, height = y, width = xwidth)) +
  geom_line(aes(x = x, y = y)) +
  geom_segment(aes(x = x0, xend = x0, y = 0, yend = y), lty = 2, size = .5) +
  geom_text(data=function(x) {x %>% group_by(cl_var, elevation) %>% 
      filter(abs(x - median(x)) == min(abs(x - median(x)))) %>% 
      select(-x) %>% distinct()}, aes(x = mean, label = label, y = y/4 * scale_height), col = 1) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(8, "RdYlBu"), drop = F) +
  theme_nothing() +
  theme(panel.background = element_rect(fill = "grey90", colour = NA)) +
  facet_wrap(trT_var ~ cl_var, scales = "free", nrow = 1)


p7 <-
  d_density %>% 
  left_join(d_mean,
            c("cl_var", "elevation", "category", "lu_var", "trT_var", "spec_var")) %>% 
  filter(category == "land_use_spec") %>%
  group_by(spec_var, lu_var) %>% 
  mutate(x0 = ifelse((any(x > 0) & any(x < 0) & y == max(y)), 0, NA),
         x0 = as.numeric(x0),
         label = formatC(mean, format = "f", digits = 3)) %>% 
  ungroup() %>% 
  mutate(x_cat = cut(x, breaks = c(-Inf, -.015, -.01, -.005, 0, .005, .01, .015, Inf))) %>% 
  ggplot() +
  geom_tile(aes(x = x, fill = x_cat, y = y/2, height = y, width = xwidth)) +
  geom_line(aes(x = x, y = y)) +
  geom_segment(aes(x = x0, xend = x0, y = 0, yend = y), lty = 2, size = .5) +
  geom_text(data=function(x) {x %>% group_by(spec_var, lu_var) %>% 
      filter(abs(x - median(x)) == min(abs(x - median(x)))) %>% 
      select(-x) %>% distinct()}, aes(x = mean, label = label, y = y/4 * scale_height), col = 1) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(8, "RdYlBu"), drop = F) +
  theme_nothing() +
  theme(panel.background = element_rect(fill = "grey90", colour = NA)) +
  facet_wrap(spec_var ~ lu_var, scales = "free", nrow = 1)

p8 <-
  d_density %>% 
  left_join(d_mean,
            c("cl_var", "elevation", "category", "lu_var", "trT_var", "spec_var")) %>% 
  filter(category == "elevation_main") %>% 
  mutate(x0 = ifelse((any(x > 0) & any(x < 0) & y == max(y)), 0, NA),
         x0 = as.numeric(x0),
         label = formatC(mean, format = "f", digits = 3)) %>% 
  mutate(x_cat = cut(x, breaks = c(-Inf, -.015, -.01, -.005, 0, .005, .01, .015, Inf))) %>% 
  ggplot() +
  geom_tile(aes(x = x, fill = x_cat, y = y/2, height = y, width = xwidth)) +
  geom_line(aes(x = x, y = y)) +
  geom_segment(aes(x = x0, xend = x0, y = 0, yend = y), lty = 2, size = .5) +
  geom_text(data=function(x) {x %>% 
      filter(abs(x - median(x)) == min(abs(x - median(x)))) %>% 
      select(-x) %>% distinct()}, aes(x = mean, label = label, y = y/4 * scale_height), col = 1) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(8, "RdYlBu"), drop = F) +
  theme_nothing() +
  theme(panel.background = element_rect(fill = "grey90", colour = NA))

p_empty <- ggplot() + theme_nothing()  

fct <- 4

ggplot_labeller <- function(labels, angle, fontsize){
  
  p <- ggplot(data.frame(x = 0, y = 0))
  
  if (angle == 0){
    for (i in seq_len(length(labels))){
      p <- p +
        annotate(geom = "text", label = labels[i], x = i, y = 0, 
                 size = fontsize/ggplot2:::.pt, angle = angle, hjust = "middle", 
                 vjust = "center",
                 parse = F)
      
    }
    p +
      theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "grey80", colour = NA),
        plot.margin = margin(0, 0, -.1, -.1, "cm"),
        plot.background = element_rect(fill = NA, colour = NA)) +
      coord_cartesian(xlim = c(.55, length(labels) + .45), expand = F)
  } else if (angle == -90){
    for (i in seq_len(length(labels))){
      p <- p +
        annotate(geom = "text", label = labels[i], x = 0, y = length(labels) - i + 1, 
                 size = fontsize/ggplot2:::.pt, angle = angle, hjust = "middle", 
                 vjust = "center",
                 parse = F)
      
    }
    p +
      theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "grey80", colour = NA),
        plot.margin = margin(0, 0, -.1, -.1, "cm"),
        plot.background = element_rect(fill = NA, colour = NA)) +
      coord_cartesian(ylim = c(.55, length(labels) + .45), expand = F)
  }
}

plot_grid(
  plot_grid(
    plot_grid(p_empty, 
              plot_grid(ggplot_labeller(c(expression(Delta*"T."~"mean"), 
                                          expression(Delta*"T."~"seasonality"),
                                          expression(Delta*"P."~"summer")),
                                        0, 7),
                        p1, nrow = 2, rel_heights = c(1, 2 * fct)),
              p_empty, ncol = 1, rel_heights = c(fct/2, 2 * fct, fct/2)),
    
    p_empty,
    
    plot_grid(plot_grid(ggplot_labeller(c(expression(Delta*"T."~"mean"), 
                                          expression(Delta*"T."~"seasonality"),
                                          expression(Delta*"P."~"summer")),
                                        0, 7),
                        p5, nrow = 2, rel_heights = c(1, 3 * fct)),
              plot_grid(p_empty, ggplot_labeller(c(expression(Delta*"Agr."~"area"), 
                                                   expression(Delta*"Grassland"~"int."),
                                                   expression(Delta*"Crop"~"int.")),
                                                 -90, 7),
                        ncol = 1,
                        rel_heights = c(1, 3*fct)),
              rel_widths = c(3 * fct, 1)),
    
    plot_grid(p_empty,
              plot_grid(ggplot_labeller(c(expression(Delta*"Agr."~"area"),
                                          expression(Delta*"Grassland"~"int."),
                                          expression(Delta*"Crop"~"int.")),
                                        0, 7),
                        plot_grid(p2a, 
                                  plot_grid(p_empty, p2b + theme(plot.margin = margin(0, 0, 0, .2, "cm")), 
                                            p_empty, ncol = 1, rel_heights = c(fct/2, fct, fct/2)),
                                  nrow = 1, rel_widths = c(2 * fct, fct)),
                        nrow = 2, rel_heights = c(1, 2 * fct)),
              p_empty, ncol = 1, rel_heights = c(fct/2, 2 * fct, fct/2)), 
    
    p_empty,
    
    nrow = 1, rel_widths = c(3 *fct, 1, 3 * fct + 1, 3 * fct, 1), 
    rel_heights = c(2 * fct + 1, fct + 1),
    scale = .9
  ),
  
  plot_grid(  
    plot_grid(plot_grid(ggplot_labeller(c(expression(Delta*"T."~"mean"), 
                                          expression(Delta*"T."~"seasonality"),
                                          expression(Delta*"P."~"summer")),
                                        0, 7),
                        p6, nrow = 2, rel_heights = c(1, 1 * fct)),
              plot_grid(p_empty, ggplot_labeller(c("Temp. niche"),
                                                 -90, 7),
                        ncol = 1,
                        rel_heights = c(1, fct)),
              rel_widths = c(3 * fct, 1)),
    
    plot_grid(p_empty, plot_grid(ggplot_labeller("Elevation", 0, 7),
                                 p8, ncol = 1, rel_heights = c(1, fct)), 
              p_empty, nrow = 1),
    
    p_empty,
    
    plot_grid(plot_grid(ggplot_labeller(c(expression(Delta*"Agr."~"area"),
                                          expression(Delta*"Grassland"~"int."),
                                          expression(Delta*"Crop"~"int.")),
                                        0, 7),
                        p7, nrow = 2, rel_heights = c(1, 1 * fct)),
              plot_grid(p_empty, ggplot_labeller("Specialisation", -90, 7),
                        ncol = 1,
                        rel_heights = c(1, fct)),
              rel_widths = c(3 * fct, 1)),
    
    nrow = 1, rel_widths = c(3 *fct + 1, 3 * fct, 1, 3 * fct + 1), 
    scale = .9
  ),
  
  plot_grid(
    plot_grid(p_empty, plot_grid(ggplot_labeller("Temp. niche", 0, 7),
                                 p3, ncol = 1, rel_heights = c(1, 2 * fct)), 
              p_empty, nrow = 1),
    
    p_empty,
    
    p_empty,
    
    p_empty,
    
    plot_grid(p_empty, plot_grid(ggplot_labeller2("Specialisation", 0, 7),
                                 p4, ncol = 1, rel_heights = c(1, 2 * fct)), 
              p_empty, nrow = 1),
    
    p_empty, 
    
    nrow = 1, rel_widths = c(3 *fct, 1, 3 * fct, 1, 3 * fct, 1), scale = .9
    
  ),
  nrow = 3, rel_heights = c(3 * fct + 1, 1 * fct + 1, 2 * fct + 1))



# plot effect size -------------------------------------------------------------.

trend_mean_sqrt <- d_trend_comb_z %>% 
  mutate(trend_mean_sqrt = trend_mean_sqrt * scaling_parameters$trend_mean_sqrt$sd +
           scaling_parameters$trend_mean_sqrt$mean) %>% 
  group_by(species) %>% 
  summarise(trend_mean_sqrt = mean(trend_mean_sqrt)) %>% 
  {mean(.$trend_mean_sqrt)}

d_gradient <- data.frame(upper = c(-.015, -.01, -.005, 0, .005, .01, .015, range_effect[2]),
                         lower = c(range_effect[1], -.015, -.01, -.005, 0, .005, .01, .015),
                         trend_mean_sqrt = trend_mean_sqrt,
                         z_mean = mean(l_trends$z_mean_reg$occ_m)) %>% 
  mutate(effect = (upper + lower) / 2,
         upper_trend_mean_sqrt = trend_mean_sqrt + upper, 
         upper_sign_trend = sign(upper_trend_mean_sqrt),
         lower_trend_mean_sqrt = trend_mean_sqrt + lower, 
         lower_sign_trend = sign(lower_trend_mean_sqrt)) %>% 
  mutate(upper_trend_mean = upper_trend_mean_sqrt ^ 2 * upper_sign_trend,
         upper_z_end = z_mean + upper_trend_mean * 5,
         lower_trend_mean = lower_trend_mean_sqrt ^ 2 * lower_sign_trend,
         lower_z_end = z_mean + lower_trend_mean * 5) 

d_gradient_sel <- data.frame(effect = c(-.015, -.01, 0, .01, .015),
                             trend_mean_sqrt = trend_mean_sqrt,
                             z_mean = mean(l_trends$z_mean_reg$occ_m)) %>% 
  mutate(trend_mean_sqrt = trend_mean_sqrt + effect,
         sign_trend = sign(trend_mean_sqrt)) %>% 
  mutate(trend_mean = trend_mean_sqrt ^ 2 * sign_trend,
         z_end = z_mean + trend_mean * 5)

p_gradient2 <-
  d_gradient %>%
  mutate(start_x = 0.2, end_x = 1,
         effect_cat = cut(effect, breaks = c(-Inf, -.015, -.01, -.005, 0, .005, .01, .015, Inf))) %>% 
  ggplot() +
  geom_rect(aes(xmin = start_x, xmax = end_x, ymin = lower_z_end, ymax = upper_z_end, fill = effect_cat, col = NULL, size = 0)) +
  geom_segment(data = d_gradient_sel,
               aes(x = .1, xend = .3, y = z_end, yend = z_end)) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(8, "RdYlBu"), drop = F) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, -.25, "cm"),
        panel.background = element_rect(fill = NA, colour = NA),
        plot.background = element_rect(fill = NA, colour = NA)) +
  coord_cartesian(xlim = c(.1, 1))

yrange <- ggplot_build(p_gradient2)$layout$panel_scales_y[[1]]$range$range

p_gradient1 <-
  d_gradient_sel %>% 
  select(effect, z_mean, z_end) %>% 
  gather(var, value, -effect) %>% 
  mutate(x = ifelse(var == "z_mean", 0, 5)) %>% 
  ggplot(aes(x = x, y = value)) +
  geom_line(aes(group = effect), size = .6) +
  geom_text(aes(label = formatC(effect, digits = 3, format = "f"), x = 6.6, hjust = 1),
            data = function(y) {y[y$x == 5, ]}, size = 7 / ggplot2:::.pt) +
  scale_colour_gradient2(limits = range_effect) +
  coord_cartesian(xlim = c(0, 6.6)) +
  scale_x_continuous(breaks = c(0, 5)) +
  scale_y_continuous(limits = yrange) +
  annotate(geom = "segment", x = -Inf, xend = 5.1, y = -Inf, yend = -Inf, size = 1) +
  theme(axis.line.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "pt"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        axis.title.x = element_text(hjust = 2.5/6.4),
        panel.background = element_rect(fill = NA, colour = NA),
        plot.background = element_rect(fill = NA, colour = NA)) +
  xlab("Years") +
  ylab("Mean occupancy")

plot_grid(p_gradient1, p_gradient2,
          nrow = 1, rel_widths = c(7, 1), align = "h")

# Trait effect sizes ###########################################################
################################################################################,

range_std_T <- (range(d_temperature_niche$Tniche_mean) - scaling_parameters$Tniche_mean$mean) / scaling_parameters$Tniche_mean$sd

range_std_sp <- (range(d_spec$spec_index) - scaling_parameters$spec_index$mean) / scaling_parameters$spec_index$sd

d_grid <- expand.grid(Tniche_mean = seq(from = range_std_T[1], to = range_std_T[2], length.out = 200),
                      spec_index = seq(from = range_std_sp[1], to = range_std_sp[2], length.out = 200))

mean_std_T <- d_mean %>% 
  filter(category == "Tniche_main") %>% 
  summarise(mean = mean(mean))

mean_std_sp <- d_mean %>% 
  filter(category == "spec_main") %>% 
  summarise(mean = mean(mean))

d_grid <- d_grid %>% 
  mutate(change = trend_mean_sqrt + Tniche_mean * mean_std_T$mean + spec_index * mean_std_sp$mean,
         sign_change = sign(change),
         change = change ^ 2 * sign(change),
         change = change * 40,
         Tniche_mean = Tniche_mean * scaling_parameters$Tniche_mean$sd + scaling_parameters$Tniche_mean$mean,
         spec_index = spec_index * scaling_parameters$spec_index$sd + scaling_parameters$spec_index$mean)

d_traits <- d_temperature_niche %>% 
  left_join(d_spec, by = c("species", "group"))

d_hull <- d_traits %>% 
  slice(chull(Tniche_mean, spec_index)) %>% 
  select(Tniche_mean, spec_index)

sf_hull <- st_multipoint(as.matrix(d_hull)) %>% 
  st_cast("POLYGON")

d_grid_sub <-
  d_grid %>% 
  mutate(X = Tniche_mean,
         Y = spec_index) %>% 
  st_as_sf(coords = c("X", "Y")) %>% 
  st_filter(sf_hull, .pred = st_intersects()) %>% 
  as.data.frame() %>% 
  select(-geometry) 

d_grid_sub %>% 
  ggplot(aes(x = Tniche_mean,
             y = spec_index)) +
  geom_tile(aes(fill = change)) +
  geom_point(data = d_traits, alpha = .5) +
  scale_fill_gradient2() +
  labs(x = "Temperature niche [°C]",
       y = "Specialisation",
       fill = "Occupancy change")

# Scenario predictions (for model version 2) ###################################
################################################################################.

# calculate predictions for different scenarios
# compare predictions to observed trends through correlations

l_fit <- extract(fit2, pars = c("mu_a", "b_T_all", "b_T_seas", "b_P",  
                               "b_AgrArea_prop", "b_LSUGL", "b_IAR", 
                               "b_Tniche", "b_spec_index", 
                               "b_elevation",
                               "b_T_all_Tniche", "b_T_seas_Tniche", "b_P_Tniche",
                               "b_AgrArea_prop_spec_index", "b_LSUGL_spec_index", "b_IAR_spec_index",
                               "b_T_all_AgrArea_prop", "b_T_seas_AgrArea_prop", "b_P_AgrArea_prop", 
                               "b_T_all_LSUGL", "b_T_seas_LSUGL", "b_P_LSUGL",
                               "b_T_all_IAR", "b_T_seas_IAR", "b_P_IAR",
                               "b_Tniche_elevation", "b_spec_index_elevation", 
                               "b_T_all_elevation", "b_T_seas_elevation", "b_P_elevation", 
                               "b_AgrArea_prop_elevation", "b_LSUGL_elevation", "b_int",
                               "alpha_spec_int", "alpha_reg", 
                               "alpha_spec_s_T_all", "alpha_spec_s_T_seas",
                               "alpha_spec_s_P", "alpha_spec_s_AgrArea_prop",
                               "alpha_spec_s_LSUGL", "alpha_spec_s_IAR",
                               "alpha_spec_s_T_all_AgrArea_prop", "alpha_spec_s_T_seas_AgrArea_prop",
                               "alpha_spec_s_P_AgrArea_prop", "alpha_spec_s_T_all_LSUGL",
                               "alpha_spec_s_T_seas_LSUGL", "alpha_spec_s_P_LSUGL",
                               "alpha_spec_s_T_all_IAR", "alpha_spec_s_T_seas_IAR",
                               "alpha_spec_s_P_IAR"))

f_pred <- function(iter_i,
                   group, spec, reg,
                   elevation, year_start,
                   T_all, T_seas, P,
                   AgrArea_prop, LSUGL, IAR,
                   Tniche, spec_index) {
  
  pred_trends <- 
    # Intercept
    l_fit$mu_a[iter_i, group] +
    # random intercepts
    l_fit$alpha_spec_int[iter_i, spec] +
    l_fit$alpha_reg[iter_i, reg] +
    # climate main effects + random slopes
    (l_fit$b_T_all[iter_i] + l_fit$alpha_spec_s_T_all[iter_i, spec]) * T_all +
    (l_fit$b_T_seas[iter_i] + l_fit$alpha_spec_s_T_seas[iter_i, spec]) * T_seas +
    (l_fit$b_P[iter_i] + l_fit$alpha_spec_s_P[iter_i, spec]) * P +
    # land use main effects + random slopes
    (l_fit$b_AgrArea_prop[iter_i] + l_fit$alpha_spec_s_AgrArea_prop[iter_i, spec]) * AgrArea_prop +
    (l_fit$b_LSUGL[iter_i] + l_fit$alpha_spec_s_LSUGL[iter_i, spec]) * LSUGL +
    (l_fit$b_IAR[iter_i] + l_fit$alpha_spec_s_IAR[iter_i, spec]) * IAR +
    # trait main effects
    l_fit$b_Tniche[iter_i] * Tniche +
    l_fit$b_spec_index[iter_i] * spec_index +
    # elevation main effect
    l_fit$b_elevation[iter_i] * elevation +
    # Climate * trait interaction
    l_fit$b_T_all_Tniche[iter_i] * T_all  * Tniche +
    l_fit$b_T_seas_Tniche[iter_i] * T_seas * Tniche +
    l_fit$b_P_Tniche[iter_i] * P * Tniche +
    # Land use * trait interaction
    l_fit$b_AgrArea_prop_spec_index[iter_i] * AgrArea_prop * spec_index +
    l_fit$b_LSUGL_spec_index[iter_i] * LSUGL * spec_index +
    l_fit$b_IAR_spec_index[iter_i] * IAR * spec_index +
    # Climate * land use interactions + random slopes
    (l_fit$b_T_all_AgrArea_prop[iter_i] + l_fit$alpha_spec_s_T_all_AgrArea_prop[iter_i, spec]) * T_all * AgrArea_prop +
    (l_fit$b_T_seas_AgrArea_prop[iter_i] + l_fit$alpha_spec_s_T_seas_AgrArea_prop[iter_i, spec]) * T_seas * AgrArea_prop +
    (l_fit$b_P_AgrArea_prop[iter_i] + l_fit$alpha_spec_s_P_AgrArea_prop[iter_i, spec]) * P * AgrArea_prop +
    (l_fit$b_T_all_LSUGL[iter_i] + l_fit$alpha_spec_s_T_all_LSUGL[iter_i, spec]) * T_all * LSUGL +
    (l_fit$b_T_seas_LSUGL[iter_i] + l_fit$alpha_spec_s_T_seas_LSUGL[iter_i, spec]) * T_seas * LSUGL +
    (l_fit$b_P_LSUGL[iter_i] + l_fit$alpha_spec_s_P_LSUGL[iter_i, spec]) * P * LSUGL +
    (l_fit$b_T_all_IAR[iter_i] + l_fit$alpha_spec_s_T_all_IAR[iter_i, spec]) * T_all * IAR +
    (l_fit$b_T_seas_IAR[iter_i] + l_fit$alpha_spec_s_T_seas_IAR[iter_i, spec]) * T_seas * IAR +
    (l_fit$b_P_IAR[iter_i] + l_fit$alpha_spec_s_P_IAR[iter_i, spec]) * P * IAR +
    # elevation * trait interaction
    l_fit$b_Tniche_elevation[iter_i] * Tniche * elevation +
    l_fit$b_spec_index_elevation[iter_i] * spec_index * elevation +
    # elevation * climate interaction
    l_fit$b_T_all_elevation[iter_i] * T_all * elevation +
    l_fit$b_T_seas_elevation[iter_i] * T_seas * elevation +
    l_fit$b_P_elevation[iter_i] * P * elevation +
    # elevation * land use interaction
    l_fit$b_AgrArea_prop_elevation[iter_i] * AgrArea_prop * elevation +
    l_fit$b_LSUGL_elevation[iter_i] * LSUGL * elevation +
    # year interval effect
    year_start %*% l_fit$b_int[iter_i, ]
  
  d_trend_comb_agrsp_z %>% 
    mutate(trend_mean_sqrt = pred_trends) %>% 
    mutate(trend_mean = trend_mean_sqrt * scaling_parameters_agrsp$trend_mean_sqrt$sd +
             scaling_parameters_agrsp$trend_mean_sqrt$mean,
           sign_trend = sign(trend_mean),
           trend_mean = sign_trend * trend_mean^2) %>% 
    group_by(group, species, zone, elevation) %>% 
    summarise(trend_sum = sum(trend_mean),
              .groups = "drop") %>% 
    left_join(d_n_squares, by = c("group", "zone")) %>% 
    group_by(species) %>% 
    summarise(trend_all = sum(trend_sum * n_squares) / sum(d_n_squares$n_squares[d_n_squares$group == v_splist[unique(species)]]),
              .groups = "drop") %>% 
    mutate(trend_all = trend_all * 5,
           run = iter_i)
}

n_iter <- 4000

# Acutal data ##################################################################.

trend_actual <- d_trend_comb_agrsp_z %>% 
  mutate(trend_mean = trend_mean_sqrt * scaling_parameters_agrsp$trend_mean_sqrt$sd +
           scaling_parameters_agrsp$trend_mean_sqrt$mean,
         sign_trend = sign(trend_mean),
         trend_mean = sign_trend * trend_mean^2) %>% 
  group_by(group, species, zone, elevation) %>% 
  summarise(trend_sum = sum(trend_mean),
            .groups = "drop") %>% 
  left_join(d_n_squares, by = c("group", "zone")) %>% 
  group_by(species) %>% 
  summarise(trend_all = sum(trend_sum * n_squares) / sum(d_n_squares$n_squares[d_n_squares$group == v_splist[unique(species)]]),
            .groups = "drop") %>% 
  mutate(trend_all = trend_all * 5)


# Actual data prediction #######################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_agrsp_z$group),
                       spec = as.numeric(d_trend_comb_agrsp_z$species),
                       reg = as.numeric(d_trend_comb_agrsp_z$zone),
                       
                       T_all = d_trend_comb_agrsp_z$T_all_change_mean,
                       T_seas = d_trend_comb_agrsp_z$T_seas_change_mean,
                       P = d_trend_comb_agrsp_z$P_warmest_quarter_change_mean,
                       AgrArea_prop = d_trend_comb_agrsp_z$AgrArea_prop_change_mean,
                       LSUGL = d_trend_comb_agrsp_z$LSUGL_change_mean,
                       IAR = d_trend_comb_agrsp_z$IAR_change_mean,
                       
                       Tniche = d_trend_comb_agrsp_z$Tniche_mean,
                       spec_index = d_trend_comb_agrsp_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_agrsp_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_agrsp_z)[, -1]) %>% 
  bind_rows() %>% 
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "actual data")

# All zero #####################################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_agrsp_z$group),
                       spec = as.numeric(d_trend_comb_agrsp_z$species),
                       reg = as.numeric(d_trend_comb_agrsp_z$zone),
                       
                       T_all = - scaling_parameters_agrsp$T_all_change_mean$mean / 
                         scaling_parameters_agrsp$T_all_change_mean$sd,
                       T_seas = - scaling_parameters_agrsp$T_seas_change_mean$mean / 
                         scaling_parameters_agrsp$T_seas_change_mean$sd,
                       P = - scaling_parameters_agrsp$P_warmest_quarter_change_mean$mean / 
                         scaling_parameters_agrsp$P_warmest_quarter_change_mean$sd,
                       AgrArea_prop = - scaling_parameters_agrsp$AgrArea_prop_change_mean$mean / 
                         scaling_parameters_agrsp$AgrArea_prop_change_mean$sd,
                       LSUGL =  - scaling_parameters_agrsp$LSUGL_change_mean$mean / 
                         scaling_parameters_agrsp$LSUGL_change_mean$sd,
                       IAR =  - scaling_parameters_agrsp$IAR_change_mean$mean / 
                         scaling_parameters_agrsp$IAR_change_mean$sd,
                       
                       Tniche = d_trend_comb_agrsp_z$Tniche_mean,
                       spec_index = d_trend_comb_agrsp_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_agrsp_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_agrsp_z)[, -1]) %>% 
  bind_rows() %>% 
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "target changes",
         cl_var = "none",
         lu_var = "none") %>% 
  bind_rows(d_trend_pred, .)

# Only T_all ###################################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_agrsp_z$group),
                       spec = as.numeric(d_trend_comb_agrsp_z$species),
                       reg = as.numeric(d_trend_comb_agrsp_z$zone),
                       
                       T_all = d_trend_comb_agrsp_z$T_all_change_mean,
                       T_seas = - scaling_parameters_agrsp$T_seas_change_mean$mean / 
                         scaling_parameters_agrsp$T_seas_change_mean$sd,
                       P = - scaling_parameters_agrsp$P_warmest_quarter_change_mean$mean / 
                         scaling_parameters_agrsp$P_warmest_quarter_change_mean$sd,
                       AgrArea_prop = - scaling_parameters_agrsp$AgrArea_prop_change_mean$mean / 
                         scaling_parameters_agrsp$AgrArea_prop_change_mean$sd,
                       LSUGL =  - scaling_parameters_agrsp$LSUGL_change_mean$mean / 
                         scaling_parameters_agrsp$LSUGL_change_mean$sd,
                       IAR =  - scaling_parameters_agrsp$IAR_change_mean$mean / 
                         scaling_parameters_agrsp$IAR_change_mean$sd,
                       
                       Tniche = d_trend_comb_agrsp_z$Tniche_mean,
                       spec_index = d_trend_comb_agrsp_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_agrsp_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_agrsp_z)[, -1]) %>% 
  bind_rows() %>% 
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "target changes",
         cl_var = "T_all",
         lu_var = "none") %>% 
  bind_rows(d_trend_pred, .)

# Only T_seas ##################################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_agrsp_z$group),
                       spec = as.numeric(d_trend_comb_agrsp_z$species),
                       reg = as.numeric(d_trend_comb_agrsp_z$zone),
                       
                       T_all = - scaling_parameters_agrsp$T_all_change_mean$mean / 
                         scaling_parameters_agrsp$T_all_change_mean$sd,
                       T_seas = d_trend_comb_agrsp_z$T_seas_change_mean,
                       P = - scaling_parameters_agrsp$P_warmest_quarter_change_mean$mean / 
                         scaling_parameters_agrsp$P_warmest_quarter_change_mean$sd,
                       AgrArea_prop = - scaling_parameters_agrsp$AgrArea_prop_change_mean$mean / 
                         scaling_parameters_agrsp$AgrArea_prop_change_mean$sd,
                       LSUGL =  - scaling_parameters_agrsp$LSUGL_change_mean$mean / 
                         scaling_parameters_agrsp$LSUGL_change_mean$sd,
                       IAR =  - scaling_parameters_agrsp$IAR_change_mean$mean / 
                         scaling_parameters_agrsp$IAR_change_mean$sd,
                       
                       Tniche = d_trend_comb_agrsp_z$Tniche_mean,
                       spec_index = d_trend_comb_agrsp_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_agrsp_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_agrsp_z)[, -1]) %>% 
  bind_rows() %>% 
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "target changes",
         cl_var = "T_seas",
         lu_var = "none") %>% 
  bind_rows(d_trend_pred, .)

# Only P #######################################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_agrsp_z$group),
                       spec = as.numeric(d_trend_comb_agrsp_z$species),
                       reg = as.numeric(d_trend_comb_agrsp_z$zone),
                       
                       T_all = - scaling_parameters_agrsp$T_all_change_mean$mean / 
                         scaling_parameters_agrsp$T_all_change_mean$sd,
                       T_seas = - scaling_parameters_agrsp$T_seas_change_mean$mean / 
                         scaling_parameters_agrsp$T_seas_change_mean$sd,
                       P = d_trend_comb_agrsp_z$P_warmest_quarter_change_mean,
                       AgrArea_prop = - scaling_parameters_agrsp$AgrArea_prop_change_mean$mean / 
                         scaling_parameters_agrsp$AgrArea_prop_change_mean$sd,
                       LSUGL =  - scaling_parameters_agrsp$LSUGL_change_mean$mean / 
                         scaling_parameters_agrsp$LSUGL_change_mean$sd,
                       IAR =  - scaling_parameters_agrsp$IAR_change_mean$mean / 
                         scaling_parameters_agrsp$IAR_change_mean$sd,
                       
                       Tniche = d_trend_comb_agrsp_z$Tniche_mean,
                       spec_index = d_trend_comb_agrsp_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_agrsp_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_agrsp_z)[, -1]) %>% 
  bind_rows() %>% 
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "target changes",
         cl_var = "P",
         lu_var = "none") %>% 
  bind_rows(d_trend_pred, .)

# Only AgrArea_prop #################################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_agrsp_z$group),
                       spec = as.numeric(d_trend_comb_agrsp_z$species),
                       reg = as.numeric(d_trend_comb_agrsp_z$zone),
                       
                       T_all = - scaling_parameters_agrsp$T_all_change_mean$mean / 
                         scaling_parameters_agrsp$T_all_change_mean$sd,
                       T_seas = - scaling_parameters_agrsp$T_seas_change_mean$mean / 
                         scaling_parameters_agrsp$T_seas_change_mean$sd,
                       P = - scaling_parameters_agrsp$P_warmest_quarter_change_mean$mean / 
                         scaling_parameters_agrsp$P_warmest_quarter_change_mean$sd,
                       AgrArea_prop = d_trend_comb_agrsp_z$AgrArea_prop_change_mean,
                       LSUGL =  - scaling_parameters_agrsp$LSUGL_change_mean$mean / 
                         scaling_parameters_agrsp$LSUGL_change_mean$sd,
                       IAR =  - scaling_parameters_agrsp$IAR_change_mean$mean / 
                         scaling_parameters_agrsp$IAR_change_mean$sd,
                       
                       Tniche = d_trend_comb_agrsp_z$Tniche_mean,
                       spec_index = d_trend_comb_agrsp_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_agrsp_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_agrsp_z)[, -1]) %>% 
  bind_rows() %>% 
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "target changes",
         lu_var = "AgrArea_prop",
         cl_var = "none") %>% 
  bind_rows(d_trend_pred, .)

# Only LSUGL ###################################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_agrsp_z$group),
                       spec = as.numeric(d_trend_comb_agrsp_z$species),
                       reg = as.numeric(d_trend_comb_agrsp_z$zone),
                       
                       T_all = - scaling_parameters_agrsp$T_all_change_mean$mean / 
                         scaling_parameters_agrsp$T_all_change_mean$sd,
                       T_seas = - scaling_parameters_agrsp$T_seas_change_mean$mean / 
                         scaling_parameters_agrsp$T_seas_change_mean$sd,
                       P = - scaling_parameters_agrsp$P_warmest_quarter_change_mean$mean / 
                         scaling_parameters_agrsp$P_warmest_quarter_change_mean$sd,
                       AgrArea_prop = - scaling_parameters_agrsp$AgrArea_prop_change_mean$mean / 
                         scaling_parameters_agrsp$AgrArea_prop_change_mean$sd,
                       LSUGL =  d_trend_comb_agrsp_z$LSUGL_change_mean,
                       IAR =  - scaling_parameters_agrsp$IAR_change_mean$mean / 
                         scaling_parameters_agrsp$IAR_change_mean$sd,
                       
                       Tniche = d_trend_comb_agrsp_z$Tniche_mean,
                       spec_index = d_trend_comb_agrsp_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_agrsp_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_agrsp_z)[, -1]) %>% 
  bind_rows() %>% 
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "target changes",
         lu_var  = "LSUGL",
         cl_var = "none") %>% 
  bind_rows(d_trend_pred, .)

# Only InsINnd #################################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_agrsp_z$group),
                       spec = as.numeric(d_trend_comb_agrsp_z$species),
                       reg = as.numeric(d_trend_comb_agrsp_z$zone),
                       
                       T_all = - scaling_parameters_agrsp$T_all_change_mean$mean / 
                         scaling_parameters_agrsp$T_all_change_mean$sd,
                       T_seas = - scaling_parameters_agrsp$T_seas_change_mean$mean / 
                         scaling_parameters_agrsp$T_seas_change_mean$sd,
                       P = - scaling_parameters_agrsp$P_warmest_quarter_change_mean$mean / 
                         scaling_parameters_agrsp$P_warmest_quarter_change_mean$sd,
                       AgrArea_prop = - scaling_parameters_agrsp$AgrArea_prop_change_mean$mean / 
                         scaling_parameters_agrsp$AgrArea_prop_change_mean$sd,
                       LSUGL =  - scaling_parameters_agrsp$LSUGL_change_mean$mean / 
                         scaling_parameters_agrsp$LSUGL_change_mean$sd,
                       IAR =  d_trend_comb_agrsp_z$IAR_change_mean,
                       
                       Tniche = d_trend_comb_agrsp_z$Tniche_mean,
                       spec_index = d_trend_comb_agrsp_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_agrsp_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_agrsp_z)[, -1]) %>% 
  bind_rows() %>% 
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "target changes",
         lu_var  = "IAR",
         cl_var = "none") %>% 
  bind_rows(d_trend_pred, .)


#  T_all  & AgrArea_prop ############################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_agrsp_z$group),
                       spec = as.numeric(d_trend_comb_agrsp_z$species),
                       reg = as.numeric(d_trend_comb_agrsp_z$zone),
                       
                       T_all = d_trend_comb_agrsp_z$T_all_change_mean,
                       T_seas = - scaling_parameters_agrsp$T_seas_change_mean$mean / 
                         scaling_parameters_agrsp$T_seas_change_mean$sd,
                       P = - scaling_parameters_agrsp$P_warmest_quarter_change_mean$mean / 
                         scaling_parameters_agrsp$P_warmest_quarter_change_mean$sd,
                       AgrArea_prop = d_trend_comb_agrsp_z$AgrArea_prop_change_mean,
                       LSUGL =  - scaling_parameters_agrsp$LSUGL_change_mean$mean / 
                         scaling_parameters_agrsp$LSUGL_change_mean$sd,
                       IAR =  - scaling_parameters_agrsp$IAR_change_mean$mean / 
                         scaling_parameters_agrsp$IAR_change_mean$sd,
                       
                       Tniche = d_trend_comb_agrsp_z$Tniche_mean,
                       spec_index = d_trend_comb_agrsp_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_agrsp_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_agrsp_z)[, -1]) %>% 
  bind_rows() %>% 
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "target changes",
         cl_var = "T_all",
         lu_var = "AgrArea_prop") %>% 
  bind_rows(d_trend_pred, .)

#  T_all  & LSUGL ##############################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_agrsp_z$group),
                       spec = as.numeric(d_trend_comb_agrsp_z$species),
                       reg = as.numeric(d_trend_comb_agrsp_z$zone),
                       
                       T_all = d_trend_comb_agrsp_z$T_all_change_mean,
                       T_seas = - scaling_parameters_agrsp$T_seas_change_mean$mean / 
                         scaling_parameters_agrsp$T_seas_change_mean$sd,
                       P = - scaling_parameters_agrsp$P_warmest_quarter_change_mean$mean / 
                         scaling_parameters_agrsp$P_warmest_quarter_change_mean$sd,
                       AgrArea_prop = - scaling_parameters_agrsp$AgrArea_prop_change_mean$mean / 
                         scaling_parameters_agrsp$AgrArea_prop_change_mean$sd,
                       LSUGL =  d_trend_comb_agrsp_z$LSUGL_change_mean,
                       IAR =  - scaling_parameters_agrsp$IAR_change_mean$mean / 
                         scaling_parameters_agrsp$IAR_change_mean$sd,
                       
                       Tniche = d_trend_comb_agrsp_z$Tniche_mean,
                       spec_index = d_trend_comb_agrsp_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_agrsp_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_agrsp_z)[, -1]) %>% 
  bind_rows() %>% 
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "target changes",
         cl_var = "T_all",
         lu_var = "LSUGL") %>% 
  bind_rows(d_trend_pred, .)

# T_all  & IAR ##############################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_agrsp_z$group),
                       spec = as.numeric(d_trend_comb_agrsp_z$species),
                       reg = as.numeric(d_trend_comb_agrsp_z$zone),
                       
                       T_all = d_trend_comb_agrsp_z$T_all_change_mean,
                       T_seas = - scaling_parameters_agrsp$T_seas_change_mean$mean / 
                         scaling_parameters_agrsp$T_seas_change_mean$sd,
                       P = - scaling_parameters_agrsp$P_warmest_quarter_change_mean$mean / 
                         scaling_parameters_agrsp$P_warmest_quarter_change_mean$sd,
                       AgrArea_prop = - scaling_parameters_agrsp$AgrArea_prop_change_mean$mean / 
                         scaling_parameters_agrsp$AgrArea_prop_change_mean$sd,
                       LSUGL =  - scaling_parameters_agrsp$LSUGL_change_mean$mean / 
                         scaling_parameters_agrsp$LSUGL_change_mean$sd,
                       IAR =  d_trend_comb_agrsp_z$IAR_change_mean,
                       
                       Tniche = d_trend_comb_agrsp_z$Tniche_mean,
                       spec_index = d_trend_comb_agrsp_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_agrsp_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_agrsp_z)[, -1]) %>% 
  bind_rows() %>% 
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "target changes",
         cl_var = "T_all",
         lu_var = "IAR") %>% 
  bind_rows(d_trend_pred, .)

# T_seas & AgrArea_prop #############################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_agrsp_z$group),
                       spec = as.numeric(d_trend_comb_agrsp_z$species),
                       reg = as.numeric(d_trend_comb_agrsp_z$zone),
                       
                       T_all = - scaling_parameters_agrsp$T_all_change_mean$mean / 
                         scaling_parameters_agrsp$T_all_change_mean$sd,
                       T_seas = d_trend_comb_agrsp_z$T_seas_change_mean,
                       P = - scaling_parameters_agrsp$P_warmest_quarter_change_mean$mean / 
                         scaling_parameters_agrsp$P_warmest_quarter_change_mean$sd,
                       AgrArea_prop = d_trend_comb_agrsp_z$AgrArea_prop_change_mean,
                       LSUGL =  - scaling_parameters_agrsp$LSUGL_change_mean$mean / 
                         scaling_parameters_agrsp$LSUGL_change_mean$sd,
                       IAR =  - scaling_parameters_agrsp$IAR_change_mean$mean / 
                         scaling_parameters_agrsp$IAR_change_mean$sd,
                       
                       Tniche = d_trend_comb_agrsp_z$Tniche_mean,
                       spec_index = d_trend_comb_agrsp_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_agrsp_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_agrsp_z)[, -1]) %>% 
  bind_rows() %>% 
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "target changes",
         cl_var = "T_seas",
         lu_var = "AgrArea_prop") %>% 
  bind_rows(d_trend_pred, .)


# T_seas & LSUGL ###############################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_agrsp_z$group),
                       spec = as.numeric(d_trend_comb_agrsp_z$species),
                       reg = as.numeric(d_trend_comb_agrsp_z$zone),
                       
                       T_all = - scaling_parameters_agrsp$T_all_change_mean$mean / 
                         scaling_parameters_agrsp$T_all_change_mean$sd,
                       T_seas = d_trend_comb_agrsp_z$T_seas_change_mean,
                       P = - scaling_parameters_agrsp$P_warmest_quarter_change_mean$mean / 
                         scaling_parameters_agrsp$P_warmest_quarter_change_mean$sd,
                       AgrArea_prop = - scaling_parameters_agrsp$AgrArea_prop_change_mean$mean / 
                         scaling_parameters_agrsp$AgrArea_prop_change_mean$sd,
                       LSUGL =  d_trend_comb_agrsp_z$LSUGL_change_mean,
                       IAR =  - scaling_parameters_agrsp$IAR_change_mean$mean / 
                         scaling_parameters_agrsp$IAR_change_mean$sd,
                       
                       Tniche = d_trend_comb_agrsp_z$Tniche_mean,
                       spec_index = d_trend_comb_agrsp_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_agrsp_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_agrsp_z)[, -1]) %>% 
  bind_rows() %>% 
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "target changes",
         cl_var = "T_seas",
         lu_var = "LSUGL") %>% 
  bind_rows(d_trend_pred, .)

# T_seas & IAR ##############################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_agrsp_z$group),
                       spec = as.numeric(d_trend_comb_agrsp_z$species),
                       reg = as.numeric(d_trend_comb_agrsp_z$zone),
                       
                       T_all = - scaling_parameters_agrsp$T_all_change_mean$mean / 
                         scaling_parameters_agrsp$T_all_change_mean$sd,
                       T_seas = d_trend_comb_agrsp_z$T_seas_change_mean,
                       P = - scaling_parameters_agrsp$P_warmest_quarter_change_mean$mean / 
                         scaling_parameters_agrsp$P_warmest_quarter_change_mean$sd,
                       AgrArea_prop = - scaling_parameters_agrsp$AgrArea_prop_change_mean$mean / 
                         scaling_parameters_agrsp$AgrArea_prop_change_mean$sd,
                       LSUGL =  - scaling_parameters_agrsp$LSUGL_change_mean$mean / 
                         scaling_parameters_agrsp$LSUGL_change_mean$sd,
                       IAR =  d_trend_comb_agrsp_z$IAR_change_mean,
                       
                       Tniche = d_trend_comb_agrsp_z$Tniche_mean,
                       spec_index = d_trend_comb_agrsp_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_agrsp_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_agrsp_z)[, -1]) %>% 
  bind_rows() %>% 
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "target changes",
         cl_var = "T_seas",
         lu_var = "IAR") %>% 
  bind_rows(d_trend_pred, .)

# P & AgrArea_prop ##################################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_agrsp_z$group),
                       spec = as.numeric(d_trend_comb_agrsp_z$species),
                       reg = as.numeric(d_trend_comb_agrsp_z$zone),
                       
                       T_all = - scaling_parameters_agrsp$T_all_change_mean$mean / 
                         scaling_parameters_agrsp$T_all_change_mean$sd,
                       T_seas = - scaling_parameters_agrsp$T_seas_change_mean$mean / 
                         scaling_parameters_agrsp$T_seas_change_mean$sd,
                       P = d_trend_comb_agrsp_z$P_warmest_quarter_change_mean,
                       AgrArea_prop = d_trend_comb_agrsp_z$AgrArea_prop_change_mean,
                       LSUGL =  - scaling_parameters_agrsp$LSUGL_change_mean$mean / 
                         scaling_parameters_agrsp$LSUGL_change_mean$sd,
                       IAR =  - scaling_parameters_agrsp$IAR_change_mean$mean / 
                         scaling_parameters_agrsp$IAR_change_mean$sd,
                       
                       Tniche = d_trend_comb_agrsp_z$Tniche_mean,
                       spec_index = d_trend_comb_agrsp_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_agrsp_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_agrsp_z)[, -1]) %>% 
  bind_rows() %>% 
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "target changes",
         cl_var = "P",
         lu_var = "AgrArea_prop") %>% 
  bind_rows(d_trend_pred, .)

# P & LSUGL ####################################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_agrsp_z$group),
                       spec = as.numeric(d_trend_comb_agrsp_z$species),
                       reg = as.numeric(d_trend_comb_agrsp_z$zone),
                       
                       T_all = - scaling_parameters_agrsp$T_all_change_mean$mean / 
                         scaling_parameters_agrsp$T_all_change_mean$sd,
                       T_seas = - scaling_parameters_agrsp$T_seas_change_mean$mean / 
                         scaling_parameters_agrsp$T_seas_change_mean$sd,
                       P = d_trend_comb_agrsp_z$P_warmest_quarter_change_mean,
                       AgrArea_prop = - scaling_parameters_agrsp$AgrArea_prop_change_mean$mean / 
                         scaling_parameters_agrsp$AgrArea_prop_change_mean$sd,
                       LSUGL =  d_trend_comb_agrsp_z$LSUGL_change_mean,
                       IAR =  - scaling_parameters_agrsp$IAR_change_mean$mean / 
                         scaling_parameters_agrsp$IAR_change_mean$sd,
                       
                       Tniche = d_trend_comb_agrsp_z$Tniche_mean,
                       spec_index = d_trend_comb_agrsp_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_agrsp_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_agrsp_z)[, -1]) %>% 
  bind_rows() %>%
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "target changes",
         cl_var = "P",
         lu_var = "LSUGL") %>% 
  bind_rows(d_trend_pred, .)

# P & IAR ###################################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_agrsp_z$group),
                       spec = as.numeric(d_trend_comb_agrsp_z$species),
                       reg = as.numeric(d_trend_comb_agrsp_z$zone),
                       
                       T_all = - scaling_parameters_agrsp$T_all_change_mean$mean / 
                         scaling_parameters_agrsp$T_all_change_mean$sd,
                       T_seas = - scaling_parameters_agrsp$T_seas_change_mean$mean / 
                         scaling_parameters_agrsp$T_seas_change_mean$sd,
                       P = d_trend_comb_agrsp_z$P_warmest_quarter_change_mean,
                       AgrArea_prop = - scaling_parameters_agrsp$AgrArea_prop_change_mean$mean / 
                         scaling_parameters_agrsp$AgrArea_prop_change_mean$sd,
                       LSUGL =  - scaling_parameters_agrsp$LSUGL_change_mean$mean / 
                         scaling_parameters_agrsp$LSUGL_change_mean$sd,
                       IAR =  d_trend_comb_agrsp_z$IAR_change_mean,
                       
                       Tniche = d_trend_comb_agrsp_z$Tniche_mean,
                       spec_index = d_trend_comb_agrsp_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_agrsp_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_agrsp_z)[, -1]) %>% 
  bind_rows() %>% 
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "target changes",
         cl_var = "P",
         lu_var = "IAR") %>% 
  bind_rows(d_trend_pred, .)

# summarise and plot scenarios -------------------------------------------------.
# ------------------------------------------------------------------------------.

d_pred_actual_data <- d_trend_pred %>% 
  filter(type == "actual data") %>% 
  mutate(R2 = cor^2) %>% 
  summarise(mean = mean(R2),
            lower80 = hdi(R2, ci = .8)$CI_low,
            upper80 = hdi(R2, ci = .8)$CI_high,
            lower95 = hdi(R2, ci = .95)$CI_low,
            upper95 = hdi(R2, ci = .95)$CI_high,
            .groups = "drop") %>% 
  mutate(cl_var = "All together", lu_var = "No change") # "No change" not true here, but needed for accurate plotting


d_trend_pred_agg <- d_trend_pred %>% 
  filter(type == "target changes") %>% 
  mutate(R2 = cor^2) %>% 
  group_by(cl_var, lu_var) %>% 
  summarise(mean = mean(R2),
            lower80 = hdi(R2, ci = .8)$CI_low,
            upper80 = hdi(R2, ci = .8)$CI_high,
            lower95 = hdi(R2, ci = .95)$CI_low,
            upper95 = hdi(R2, ci = .95)$CI_high,
            .groups = "drop") %>% 
  add_row(mean = NA,
          cl_var = "All together",
          lu_var = "none") %>%
  mutate(cl_var = factor(cl_var, levels = c("none", "T_all", "T_seas", "P", "All together"),
                         labels = c("No change", "Annual mean\ntemperature", "Temperature\nseasonality", 
                                    "Summer\nprecipitation", "All together")),
         lu_var = factor(lu_var, levels = c("none", "AgrArea_prop", "LSUGL", "IAR"),
                         labels = c("No change", "Agricultural\narea", "Grassland-use\nintensity", "Crop-use\nintensity")))

d_trend_pred_agg %>%   
  ggplot(aes(x = cl_var, y = mean, fill = lu_var)) +
  geom_point(col = NA) +
  annotate(geom = "rect", xmin = 4.5, xmax = 5.5, ymin = -Inf, ymax = Inf,
           fill = "grey90") +
  annotate(geom = "segment", x = 4.5, xend = 4.5, y = -.042, 
           yend = Inf,
           col = "grey20") +
  geom_col(position = position_dodge(width = .8), width = .8) +
  geom_col(data = d_pred_actual_data,
           fill = "#45773C", width = .8/4) +
  geom_hline(yintercept = 0) +
  geom_linerange(aes(ymin = lower80, ymax = upper80),
                 position = position_dodge(width = .8), size = 1) +
  geom_linerange(aes(ymin = lower95, ymax = upper95),
                 position = position_dodge(width = .8), size = .5) +
  geom_linerange(data = d_pred_actual_data, aes(ymin = lower80, ymax = upper80),
                 size = 1) +
  geom_linerange(data = d_pred_actual_data, aes(ymin = lower95, ymax = upper95),
                 size = .5) +
  scale_fill_manual(values = c("grey30", "#5158BB", "#F26DF9", "#EB4B98"),
                    name = "Land-use change") +
  xlab("Climate change") +
  ylab(expression(Mean~italic(R^2))) +
  coord_cartesian(clip = "off",
                  ylim = c(0, max(c(d_trend_pred_agg$upper95,
                                    d_pred_actual_data$upper95), na.rm = T))) +
  theme(axis.title.x = element_text(hjust = 2.5 / 6.6),
        legend.text = element_text(margin = margin(t = 5, b = 5, unit = "pt")),
        legend.key.height = unit(30, "pt"))

# Scenario predictions (for model version 3) ###################################
################################################################################.

# calculate predictions for different scenarios
# compare predictions to observed trends through correlations

l_fit <- extract(fit3, pars = c("mu_a", "b_T_all", "b_T_seas", "b_P",  
                                "b_AgrArea_prop", "b_LSUGL", "b_IAR", 
                                "b_Tniche", "b_spec_index", 
                                "b_elevation",
                                "b_T_all_Tniche", "b_T_seas_Tniche", "b_P_Tniche",
                                "b_AgrArea_prop_spec_index", "b_LSUGL_spec_index", "b_IAR_spec_index",
                                "b_T_all_AgrArea_prop", "b_T_seas_AgrArea_prop", "b_P_AgrArea_prop", 
                                "b_T_all_LSUGL", "b_T_seas_LSUGL", "b_P_LSUGL",
                                "b_T_all_IAR", "b_T_seas_IAR", "b_P_IAR",
                                "b_Tniche_elevation", "b_spec_index_elevation", 
                                "b_T_all_elevation", "b_T_seas_elevation", "b_P_elevation", 
                                "b_AgrArea_prop_elevation", "b_LSUGL_elevation", "b_int",
                                "alpha_spec_int", "alpha_reg", 
                                "alpha_spec_s_T_all", "alpha_spec_s_T_seas",
                                "alpha_spec_s_P", "alpha_spec_s_AgrArea_prop",
                                "alpha_spec_s_LSUGL", "alpha_spec_s_IAR",
                                "alpha_spec_s_T_all_AgrArea_prop", "alpha_spec_s_T_seas_AgrArea_prop",
                                "alpha_spec_s_P_AgrArea_prop", "alpha_spec_s_T_all_LSUGL",
                                "alpha_spec_s_T_seas_LSUGL", "alpha_spec_s_P_LSUGL",
                                "alpha_spec_s_T_all_IAR", "alpha_spec_s_T_seas_IAR",
                                "alpha_spec_s_P_IAR"))

f_pred <- function(iter_i,
                   group, spec, reg,
                   elevation, year_start,
                   T_all, T_seas, P,
                   AgrArea_prop, LSUGL, IAR,
                   Tniche, spec_index) {
  
  pred_trends <- 
    # Intercept
    l_fit$mu_a[iter_i, group] +
    # random intercepts
    l_fit$alpha_spec_int[iter_i, spec] +
    l_fit$alpha_reg[iter_i, reg] +
    # climate main effects + random slopes
    (l_fit$b_T_all[iter_i] + l_fit$alpha_spec_s_T_all[iter_i, spec]) * T_all +
    (l_fit$b_T_seas[iter_i] + l_fit$alpha_spec_s_T_seas[iter_i, spec]) * T_seas +
    (l_fit$b_P[iter_i] + l_fit$alpha_spec_s_P[iter_i, spec]) * P +
    # land use main effects + random slopes
    (l_fit$b_AgrArea_prop[iter_i] + l_fit$alpha_spec_s_AgrArea_prop[iter_i, spec]) * AgrArea_prop +
    (l_fit$b_LSUGL[iter_i] + l_fit$alpha_spec_s_LSUGL[iter_i, spec]) * LSUGL +
    (l_fit$b_IAR[iter_i] + l_fit$alpha_spec_s_IAR[iter_i, spec]) * IAR +
    # trait main effects
    l_fit$b_Tniche[iter_i] * Tniche +
    l_fit$b_spec_index[iter_i] * spec_index +
    # elevation main effect
    l_fit$b_elevation[iter_i] * elevation +
    # Climate * trait interaction
    l_fit$b_T_all_Tniche[iter_i] * T_all  * Tniche +
    l_fit$b_T_seas_Tniche[iter_i] * T_seas * Tniche +
    l_fit$b_P_Tniche[iter_i] * P * Tniche +
    # Land use * trait interaction
    l_fit$b_AgrArea_prop_spec_index[iter_i] * AgrArea_prop * spec_index +
    l_fit$b_LSUGL_spec_index[iter_i] * LSUGL * spec_index +
    l_fit$b_IAR_spec_index[iter_i] * IAR * spec_index +
    # Climate * land use interactions + random slopes
    (l_fit$b_T_all_AgrArea_prop[iter_i] + l_fit$alpha_spec_s_T_all_AgrArea_prop[iter_i, spec]) * T_all * AgrArea_prop +
    (l_fit$b_T_seas_AgrArea_prop[iter_i] + l_fit$alpha_spec_s_T_seas_AgrArea_prop[iter_i, spec]) * T_seas * AgrArea_prop +
    (l_fit$b_P_AgrArea_prop[iter_i] + l_fit$alpha_spec_s_P_AgrArea_prop[iter_i, spec]) * P * AgrArea_prop +
    (l_fit$b_T_all_LSUGL[iter_i] + l_fit$alpha_spec_s_T_all_LSUGL[iter_i, spec]) * T_all * LSUGL +
    (l_fit$b_T_seas_LSUGL[iter_i] + l_fit$alpha_spec_s_T_seas_LSUGL[iter_i, spec]) * T_seas * LSUGL +
    (l_fit$b_P_LSUGL[iter_i] + l_fit$alpha_spec_s_P_LSUGL[iter_i, spec]) * P * LSUGL +
    (l_fit$b_T_all_IAR[iter_i] + l_fit$alpha_spec_s_T_all_IAR[iter_i, spec]) * T_all * IAR +
    (l_fit$b_T_seas_IAR[iter_i] + l_fit$alpha_spec_s_T_seas_IAR[iter_i, spec]) * T_seas * IAR +
    (l_fit$b_P_IAR[iter_i] + l_fit$alpha_spec_s_P_IAR[iter_i, spec]) * P * IAR +
    # elevation * trait interaction
    l_fit$b_Tniche_elevation[iter_i] * Tniche * elevation +
    l_fit$b_spec_index_elevation[iter_i] * spec_index * elevation +
    # elevation * climate interaction
    l_fit$b_T_all_elevation[iter_i] * T_all * elevation +
    l_fit$b_T_seas_elevation[iter_i] * T_seas * elevation +
    l_fit$b_P_elevation[iter_i] * P * elevation +
    # elevation * land use interaction
    l_fit$b_AgrArea_prop_elevation[iter_i] * AgrArea_prop * elevation +
    l_fit$b_LSUGL_elevation[iter_i] * LSUGL * elevation +
    # year interval effect
    year_start %*% l_fit$b_int[iter_i, ]
  
  d_trend_comb_z %>% 
    mutate(trend_mean_sqrt = pred_trends) %>% 
    mutate(trend_mean = trend_mean_sqrt * scaling_parameters$trend_mean_sqrt$sd +
             scaling_parameters$trend_mean_sqrt$mean,
           sign_trend = sign(trend_mean),
           trend_mean = sign_trend * trend_mean^2) %>% 
    group_by(group, species, zone, elevation) %>% 
    summarise(trend_sum = sum(trend_mean),
              .groups = "drop") %>% 
    left_join(d_n_squares, by = c("group", "zone")) %>% 
    group_by(species) %>% 
    summarise(trend_all = sum(trend_sum * n_squares) / sum(d_n_squares$n_squares[d_n_squares$group == v_splist[unique(species)]]),
              .groups = "drop") %>% 
    mutate(trend_all = trend_all * 5,
           run = iter_i)
}

n_iter <- 4000

# Acutal data ##################################################################.

trend_actual <- d_trend_comb_z %>% 
  mutate(trend_mean = trend_mean_sqrt * scaling_parameters$trend_mean_sqrt$sd +
           scaling_parameters$trend_mean_sqrt$mean,
         sign_trend = sign(trend_mean),
         trend_mean = sign_trend * trend_mean^2) %>% 
  group_by(group, species, zone, elevation) %>% 
  summarise(trend_sum = sum(trend_mean),
            .groups = "drop") %>% 
  left_join(d_n_squares, by = c("group", "zone")) %>% 
  group_by(species) %>% 
  summarise(trend_all = sum(trend_sum * n_squares) / sum(d_n_squares$n_squares[d_n_squares$group == v_splist[unique(species)]]),
            .groups = "drop") %>% 
  mutate(trend_all = trend_all * 5)


# Actual data prediction #######################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_z$group),
                       spec = as.numeric(d_trend_comb_z$species),
                       reg = as.numeric(d_trend_comb_z$zone),
                       
                       T_all = d_trend_comb_z$T_all_change_mean,
                       T_seas = d_trend_comb_z$T_seas_change_mean,
                       P = d_trend_comb_z$P_warmest_quarter_change_mean,
                       AgrArea_prop = d_trend_comb_z$AgrArea_prop_change_mean,
                       LSUGL = d_trend_comb_z$LSUGL_change_mean,
                       IAR = d_trend_comb_z$IAR_change_mean,
                       
                       Tniche = d_trend_comb_z$Tniche_mean,
                       spec_index = d_trend_comb_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_z)[, -1]) %>% 
  bind_rows() %>% 
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "actual data")

# All zero #####################################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_z$group),
                       spec = as.numeric(d_trend_comb_z$species),
                       reg = as.numeric(d_trend_comb_z$zone),
                       
                       T_all = - scaling_parameters$T_all_change_mean$mean / 
                         scaling_parameters$T_all_change_mean$sd,
                       T_seas = - scaling_parameters$T_seas_change_mean$mean / 
                         scaling_parameters$T_seas_change_mean$sd,
                       P = - scaling_parameters$P_warmest_quarter_change_mean$mean / 
                         scaling_parameters$P_warmest_quarter_change_mean$sd,
                       AgrArea_prop = - scaling_parameters$AgrArea_prop_change_mean$mean / 
                         scaling_parameters$AgrArea_prop_change_mean$sd,
                       LSUGL =  - scaling_parameters$LSUGL_change_mean$mean / 
                         scaling_parameters$LSUGL_change_mean$sd,
                       IAR =  - scaling_parameters$IAR_change_mean$mean / 
                         scaling_parameters$IAR_change_mean$sd,
                       
                       Tniche = d_trend_comb_z$Tniche_mean,
                       spec_index = d_trend_comb_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_z)[, -1]) %>% 
  bind_rows() %>% 
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "target changes",
         cl_var = "none",
         lu_var = "none") %>% 
  bind_rows(d_trend_pred, .)

# Only T_all ###################################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_z$group),
                       spec = as.numeric(d_trend_comb_z$species),
                       reg = as.numeric(d_trend_comb_z$zone),
                       
                       T_all = d_trend_comb_z$T_all_change_mean,
                       T_seas = - scaling_parameters$T_seas_change_mean$mean / 
                         scaling_parameters$T_seas_change_mean$sd,
                       P = - scaling_parameters$P_warmest_quarter_change_mean$mean / 
                         scaling_parameters$P_warmest_quarter_change_mean$sd,
                       AgrArea_prop = - scaling_parameters$AgrArea_prop_change_mean$mean / 
                         scaling_parameters$AgrArea_prop_change_mean$sd,
                       LSUGL =  - scaling_parameters$LSUGL_change_mean$mean / 
                         scaling_parameters$LSUGL_change_mean$sd,
                       IAR =  - scaling_parameters$IAR_change_mean$mean / 
                         scaling_parameters$IAR_change_mean$sd,
                       
                       Tniche = d_trend_comb_z$Tniche_mean,
                       spec_index = d_trend_comb_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_z)[, -1]) %>% 
  bind_rows() %>% 
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "target changes",
         cl_var = "T_all",
         lu_var = "none") %>% 
  bind_rows(d_trend_pred, .)

# Only T_seas ##################################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_z$group),
                       spec = as.numeric(d_trend_comb_z$species),
                       reg = as.numeric(d_trend_comb_z$zone),
                       
                       T_all = - scaling_parameters$T_all_change_mean$mean / 
                         scaling_parameters$T_all_change_mean$sd,
                       T_seas = d_trend_comb_z$T_seas_change_mean,
                       P = - scaling_parameters$P_warmest_quarter_change_mean$mean / 
                         scaling_parameters$P_warmest_quarter_change_mean$sd,
                       AgrArea_prop = - scaling_parameters$AgrArea_prop_change_mean$mean / 
                         scaling_parameters$AgrArea_prop_change_mean$sd,
                       LSUGL =  - scaling_parameters$LSUGL_change_mean$mean / 
                         scaling_parameters$LSUGL_change_mean$sd,
                       IAR =  - scaling_parameters$IAR_change_mean$mean / 
                         scaling_parameters$IAR_change_mean$sd,
                       
                       Tniche = d_trend_comb_z$Tniche_mean,
                       spec_index = d_trend_comb_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_z)[, -1]) %>% 
  bind_rows() %>% 
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "target changes",
         cl_var = "T_seas",
         lu_var = "none") %>% 
  bind_rows(d_trend_pred, .)

# Only P #######################################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_z$group),
                       spec = as.numeric(d_trend_comb_z$species),
                       reg = as.numeric(d_trend_comb_z$zone),
                       
                       T_all = - scaling_parameters$T_all_change_mean$mean / 
                         scaling_parameters$T_all_change_mean$sd,
                       T_seas = - scaling_parameters$T_seas_change_mean$mean / 
                         scaling_parameters$T_seas_change_mean$sd,
                       P = d_trend_comb_z$P_warmest_quarter_change_mean,
                       AgrArea_prop = - scaling_parameters$AgrArea_prop_change_mean$mean / 
                         scaling_parameters$AgrArea_prop_change_mean$sd,
                       LSUGL =  - scaling_parameters$LSUGL_change_mean$mean / 
                         scaling_parameters$LSUGL_change_mean$sd,
                       IAR =  - scaling_parameters$IAR_change_mean$mean / 
                         scaling_parameters$IAR_change_mean$sd,
                       
                       Tniche = d_trend_comb_z$Tniche_mean,
                       spec_index = d_trend_comb_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_z)[, -1]) %>% 
  bind_rows() %>% 
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "target changes",
         cl_var = "P",
         lu_var = "none") %>% 
  bind_rows(d_trend_pred, .)

# Only AgrArea_prop #################################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_z$group),
                       spec = as.numeric(d_trend_comb_z$species),
                       reg = as.numeric(d_trend_comb_z$zone),
                       
                       T_all = - scaling_parameters$T_all_change_mean$mean / 
                         scaling_parameters$T_all_change_mean$sd,
                       T_seas = - scaling_parameters$T_seas_change_mean$mean / 
                         scaling_parameters$T_seas_change_mean$sd,
                       P = - scaling_parameters$P_warmest_quarter_change_mean$mean / 
                         scaling_parameters$P_warmest_quarter_change_mean$sd,
                       AgrArea_prop = d_trend_comb_z$AgrArea_prop_change_mean,
                       LSUGL =  - scaling_parameters$LSUGL_change_mean$mean / 
                         scaling_parameters$LSUGL_change_mean$sd,
                       IAR =  - scaling_parameters$IAR_change_mean$mean / 
                         scaling_parameters$IAR_change_mean$sd,
                       
                       Tniche = d_trend_comb_z$Tniche_mean,
                       spec_index = d_trend_comb_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_z)[, -1]) %>% 
  bind_rows() %>% 
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "target changes",
         lu_var = "AgrArea_prop",
         cl_var = "none") %>% 
  bind_rows(d_trend_pred, .)

# Only LSUGL ###################################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_z$group),
                       spec = as.numeric(d_trend_comb_z$species),
                       reg = as.numeric(d_trend_comb_z$zone),
                       
                       T_all = - scaling_parameters$T_all_change_mean$mean / 
                         scaling_parameters$T_all_change_mean$sd,
                       T_seas = - scaling_parameters$T_seas_change_mean$mean / 
                         scaling_parameters$T_seas_change_mean$sd,
                       P = - scaling_parameters$P_warmest_quarter_change_mean$mean / 
                         scaling_parameters$P_warmest_quarter_change_mean$sd,
                       AgrArea_prop = - scaling_parameters$AgrArea_prop_change_mean$mean / 
                         scaling_parameters$AgrArea_prop_change_mean$sd,
                       LSUGL =  d_trend_comb_z$LSUGL_change_mean,
                       IAR =  - scaling_parameters$IAR_change_mean$mean / 
                         scaling_parameters$IAR_change_mean$sd,
                       
                       Tniche = d_trend_comb_z$Tniche_mean,
                       spec_index = d_trend_comb_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_z)[, -1]) %>% 
  bind_rows() %>% 
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "target changes",
         lu_var  = "LSUGL",
         cl_var = "none") %>% 
  bind_rows(d_trend_pred, .)

# Only InsINnd #################################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_z$group),
                       spec = as.numeric(d_trend_comb_z$species),
                       reg = as.numeric(d_trend_comb_z$zone),
                       
                       T_all = - scaling_parameters$T_all_change_mean$mean / 
                         scaling_parameters$T_all_change_mean$sd,
                       T_seas = - scaling_parameters$T_seas_change_mean$mean / 
                         scaling_parameters$T_seas_change_mean$sd,
                       P = - scaling_parameters$P_warmest_quarter_change_mean$mean / 
                         scaling_parameters$P_warmest_quarter_change_mean$sd,
                       AgrArea_prop = - scaling_parameters$AgrArea_prop_change_mean$mean / 
                         scaling_parameters$AgrArea_prop_change_mean$sd,
                       LSUGL =  - scaling_parameters$LSUGL_change_mean$mean / 
                         scaling_parameters$LSUGL_change_mean$sd,
                       IAR =  d_trend_comb_z$IAR_change_mean,
                       
                       Tniche = d_trend_comb_z$Tniche_mean,
                       spec_index = d_trend_comb_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_z)[, -1]) %>% 
  bind_rows() %>% 
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "target changes",
         lu_var  = "IAR",
         cl_var = "none") %>% 
  bind_rows(d_trend_pred, .)


#  T_all  & AgrArea_prop ############################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_z$group),
                       spec = as.numeric(d_trend_comb_z$species),
                       reg = as.numeric(d_trend_comb_z$zone),
                       
                       T_all = d_trend_comb_z$T_all_change_mean,
                       T_seas = - scaling_parameters$T_seas_change_mean$mean / 
                         scaling_parameters$T_seas_change_mean$sd,
                       P = - scaling_parameters$P_warmest_quarter_change_mean$mean / 
                         scaling_parameters$P_warmest_quarter_change_mean$sd,
                       AgrArea_prop = d_trend_comb_z$AgrArea_prop_change_mean,
                       LSUGL =  - scaling_parameters$LSUGL_change_mean$mean / 
                         scaling_parameters$LSUGL_change_mean$sd,
                       IAR =  - scaling_parameters$IAR_change_mean$mean / 
                         scaling_parameters$IAR_change_mean$sd,
                       
                       Tniche = d_trend_comb_z$Tniche_mean,
                       spec_index = d_trend_comb_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_z)[, -1]) %>% 
  bind_rows() %>% 
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "target changes",
         cl_var = "T_all",
         lu_var = "AgrArea_prop") %>% 
  bind_rows(d_trend_pred, .)

#  T_all  & LSUGL ##############################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_z$group),
                       spec = as.numeric(d_trend_comb_z$species),
                       reg = as.numeric(d_trend_comb_z$zone),
                       
                       T_all = d_trend_comb_z$T_all_change_mean,
                       T_seas = - scaling_parameters$T_seas_change_mean$mean / 
                         scaling_parameters$T_seas_change_mean$sd,
                       P = - scaling_parameters$P_warmest_quarter_change_mean$mean / 
                         scaling_parameters$P_warmest_quarter_change_mean$sd,
                       AgrArea_prop = - scaling_parameters$AgrArea_prop_change_mean$mean / 
                         scaling_parameters$AgrArea_prop_change_mean$sd,
                       LSUGL =  d_trend_comb_z$LSUGL_change_mean,
                       IAR =  - scaling_parameters$IAR_change_mean$mean / 
                         scaling_parameters$IAR_change_mean$sd,
                       
                       Tniche = d_trend_comb_z$Tniche_mean,
                       spec_index = d_trend_comb_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_z)[, -1]) %>% 
  bind_rows() %>% 
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "target changes",
         cl_var = "T_all",
         lu_var = "LSUGL") %>% 
  bind_rows(d_trend_pred, .)

# T_all  & IAR ##############################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_z$group),
                       spec = as.numeric(d_trend_comb_z$species),
                       reg = as.numeric(d_trend_comb_z$zone),
                       
                       T_all = d_trend_comb_z$T_all_change_mean,
                       T_seas = - scaling_parameters$T_seas_change_mean$mean / 
                         scaling_parameters$T_seas_change_mean$sd,
                       P = - scaling_parameters$P_warmest_quarter_change_mean$mean / 
                         scaling_parameters$P_warmest_quarter_change_mean$sd,
                       AgrArea_prop = - scaling_parameters$AgrArea_prop_change_mean$mean / 
                         scaling_parameters$AgrArea_prop_change_mean$sd,
                       LSUGL =  - scaling_parameters$LSUGL_change_mean$mean / 
                         scaling_parameters$LSUGL_change_mean$sd,
                       IAR =  d_trend_comb_z$IAR_change_mean,
                       
                       Tniche = d_trend_comb_z$Tniche_mean,
                       spec_index = d_trend_comb_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_z)[, -1]) %>% 
  bind_rows() %>% 
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "target changes",
         cl_var = "T_all",
         lu_var = "IAR") %>% 
  bind_rows(d_trend_pred, .)

# T_seas & AgrArea_prop #############################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_z$group),
                       spec = as.numeric(d_trend_comb_z$species),
                       reg = as.numeric(d_trend_comb_z$zone),
                       
                       T_all = - scaling_parameters$T_all_change_mean$mean / 
                         scaling_parameters$T_all_change_mean$sd,
                       T_seas = d_trend_comb_z$T_seas_change_mean,
                       P = - scaling_parameters$P_warmest_quarter_change_mean$mean / 
                         scaling_parameters$P_warmest_quarter_change_mean$sd,
                       AgrArea_prop = d_trend_comb_z$AgrArea_prop_change_mean,
                       LSUGL =  - scaling_parameters$LSUGL_change_mean$mean / 
                         scaling_parameters$LSUGL_change_mean$sd,
                       IAR =  - scaling_parameters$IAR_change_mean$mean / 
                         scaling_parameters$IAR_change_mean$sd,
                       
                       Tniche = d_trend_comb_z$Tniche_mean,
                       spec_index = d_trend_comb_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_z)[, -1]) %>% 
  bind_rows() %>% 
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "target changes",
         cl_var = "T_seas",
         lu_var = "AgrArea_prop") %>% 
  bind_rows(d_trend_pred, .)


# T_seas & LSUGL ###############################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_z$group),
                       spec = as.numeric(d_trend_comb_z$species),
                       reg = as.numeric(d_trend_comb_z$zone),
                       
                       T_all = - scaling_parameters$T_all_change_mean$mean / 
                         scaling_parameters$T_all_change_mean$sd,
                       T_seas = d_trend_comb_z$T_seas_change_mean,
                       P = - scaling_parameters$P_warmest_quarter_change_mean$mean / 
                         scaling_parameters$P_warmest_quarter_change_mean$sd,
                       AgrArea_prop = - scaling_parameters$AgrArea_prop_change_mean$mean / 
                         scaling_parameters$AgrArea_prop_change_mean$sd,
                       LSUGL =  d_trend_comb_z$LSUGL_change_mean,
                       IAR =  - scaling_parameters$IAR_change_mean$mean / 
                         scaling_parameters$IAR_change_mean$sd,
                       
                       Tniche = d_trend_comb_z$Tniche_mean,
                       spec_index = d_trend_comb_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_z)[, -1]) %>% 
  bind_rows() %>% 
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "target changes",
         cl_var = "T_seas",
         lu_var = "LSUGL") %>% 
  bind_rows(d_trend_pred, .)

# T_seas & IAR ##############################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_z$group),
                       spec = as.numeric(d_trend_comb_z$species),
                       reg = as.numeric(d_trend_comb_z$zone),
                       
                       T_all = - scaling_parameters$T_all_change_mean$mean / 
                         scaling_parameters$T_all_change_mean$sd,
                       T_seas = d_trend_comb_z$T_seas_change_mean,
                       P = - scaling_parameters$P_warmest_quarter_change_mean$mean / 
                         scaling_parameters$P_warmest_quarter_change_mean$sd,
                       AgrArea_prop = - scaling_parameters$AgrArea_prop_change_mean$mean / 
                         scaling_parameters$AgrArea_prop_change_mean$sd,
                       LSUGL =  - scaling_parameters$LSUGL_change_mean$mean / 
                         scaling_parameters$LSUGL_change_mean$sd,
                       IAR =  d_trend_comb_z$IAR_change_mean,
                       
                       Tniche = d_trend_comb_z$Tniche_mean,
                       spec_index = d_trend_comb_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_z)[, -1]) %>% 
  bind_rows() %>% 
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "target changes",
         cl_var = "T_seas",
         lu_var = "IAR") %>% 
  bind_rows(d_trend_pred, .)

# P & AgrArea_prop ##################################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_z$group),
                       spec = as.numeric(d_trend_comb_z$species),
                       reg = as.numeric(d_trend_comb_z$zone),
                       
                       T_all = - scaling_parameters$T_all_change_mean$mean / 
                         scaling_parameters$T_all_change_mean$sd,
                       T_seas = - scaling_parameters$T_seas_change_mean$mean / 
                         scaling_parameters$T_seas_change_mean$sd,
                       P = d_trend_comb_z$P_warmest_quarter_change_mean,
                       AgrArea_prop = d_trend_comb_z$AgrArea_prop_change_mean,
                       LSUGL =  - scaling_parameters$LSUGL_change_mean$mean / 
                         scaling_parameters$LSUGL_change_mean$sd,
                       IAR =  - scaling_parameters$IAR_change_mean$mean / 
                         scaling_parameters$IAR_change_mean$sd,
                       
                       Tniche = d_trend_comb_z$Tniche_mean,
                       spec_index = d_trend_comb_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_z)[, -1]) %>% 
  bind_rows() %>% 
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "target changes",
         cl_var = "P",
         lu_var = "AgrArea_prop") %>% 
  bind_rows(d_trend_pred, .)

# P & LSUGL ####################################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_z$group),
                       spec = as.numeric(d_trend_comb_z$species),
                       reg = as.numeric(d_trend_comb_z$zone),
                       
                       T_all = - scaling_parameters$T_all_change_mean$mean / 
                         scaling_parameters$T_all_change_mean$sd,
                       T_seas = - scaling_parameters$T_seas_change_mean$mean / 
                         scaling_parameters$T_seas_change_mean$sd,
                       P = d_trend_comb_z$P_warmest_quarter_change_mean,
                       AgrArea_prop = - scaling_parameters$AgrArea_prop_change_mean$mean / 
                         scaling_parameters$AgrArea_prop_change_mean$sd,
                       LSUGL =  d_trend_comb_z$LSUGL_change_mean,
                       IAR =  - scaling_parameters$IAR_change_mean$mean / 
                         scaling_parameters$IAR_change_mean$sd,
                       
                       Tniche = d_trend_comb_z$Tniche_mean,
                       spec_index = d_trend_comb_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_z)[, -1]) %>% 
  bind_rows() %>%
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "target changes",
         cl_var = "P",
         lu_var = "LSUGL") %>% 
  bind_rows(d_trend_pred, .)

# P & IAR ###################################################################.

d_trend_pred <- lapply(1:n_iter, 
                       f_pred,
                       group = as.numeric(d_trend_comb_z$group),
                       spec = as.numeric(d_trend_comb_z$species),
                       reg = as.numeric(d_trend_comb_z$zone),
                       
                       T_all = - scaling_parameters$T_all_change_mean$mean / 
                         scaling_parameters$T_all_change_mean$sd,
                       T_seas = - scaling_parameters$T_seas_change_mean$mean / 
                         scaling_parameters$T_seas_change_mean$sd,
                       P = d_trend_comb_z$P_warmest_quarter_change_mean,
                       AgrArea_prop = - scaling_parameters$AgrArea_prop_change_mean$mean / 
                         scaling_parameters$AgrArea_prop_change_mean$sd,
                       LSUGL =  - scaling_parameters$LSUGL_change_mean$mean / 
                         scaling_parameters$LSUGL_change_mean$sd,
                       IAR =  d_trend_comb_z$IAR_change_mean,
                       
                       Tniche = d_trend_comb_z$Tniche_mean,
                       spec_index = d_trend_comb_z$spec_index,
                       
                       elevation = as.numeric(d_trend_comb_z$elevation == "high"),
                       
                       year_start = model.matrix(~ year_start_f, data = d_trend_comb_z)[, -1]) %>% 
  bind_rows() %>% 
  left_join(trend_actual, by = "species",
            suffix = c(".pred", ".obs")) %>% 
  group_by(run) %>% 
  group_split() %>%
  map_dfr(~ data.frame(run = unique(.$run),
                       cor = cor(.$trend_all.obs, .$trend_all.pred))) %>% 
  mutate(type = "target changes",
         cl_var = "P",
         lu_var = "IAR") %>% 
  bind_rows(d_trend_pred, .)

# summarise and plot scenarios -------------------------------------------------.
# ------------------------------------------------------------------------------.

d_pred_actual_data <- d_trend_pred %>% 
  filter(type == "actual data") %>% 
  mutate(R2 = cor^2) %>% 
  summarise(mean = mean(R2),
            lower80 = hdi(R2, ci = .8)$CI_low,
            upper80 = hdi(R2, ci = .8)$CI_high,
            lower95 = hdi(R2, ci = .95)$CI_low,
            upper95 = hdi(R2, ci = .95)$CI_high,
            .groups = "drop") %>% 
  mutate(cl_var = "All together", lu_var = "No change") # "No change" not true here, but needed for accurate plotting


d_trend_pred_agg <- d_trend_pred %>% 
  filter(type == "target changes") %>% 
  mutate(R2 = cor^2) %>% 
  group_by(cl_var, lu_var) %>% 
  summarise(mean = mean(R2),
            lower80 = hdi(R2, ci = .8)$CI_low,
            upper80 = hdi(R2, ci = .8)$CI_high,
            lower95 = hdi(R2, ci = .95)$CI_low,
            upper95 = hdi(R2, ci = .95)$CI_high,
            .groups = "drop") %>% 
  add_row(mean = NA,
          cl_var = "All together",
          lu_var = "none") %>%
  mutate(cl_var = factor(cl_var, levels = c("none", "T_all", "T_seas", "P", "All together"),
                         labels = c("No change", "Annual mean\ntemperature", "Temperature\nseasonality", 
                                    "Summer\nprecipitation", "All together")),
         lu_var = factor(lu_var, levels = c("none", "AgrArea_prop", "LSUGL", "IAR"),
                         labels = c("No change", "Agricultural\narea", "Grassland-use\nintensity", "Crop-use\nintensity")))

d_trend_pred_agg %>%   
  ggplot(aes(x = cl_var, y = mean, fill = lu_var)) +
  geom_point(col = NA) +
  annotate(geom = "rect", xmin = 4.5, xmax = 5.5, ymin = -Inf, ymax = Inf,
           fill = "grey90") +
  annotate(geom = "segment", x = 4.5, xend = 4.5, y = -.042, 
           yend = Inf,
           col = "grey20") +
  geom_col(position = position_dodge(width = .8), width = .8) +
  geom_col(data = d_pred_actual_data,
           fill = "#45773C", width = .8/4) +
  geom_hline(yintercept = 0) +
  geom_linerange(aes(ymin = lower80, ymax = upper80),
                 position = position_dodge(width = .8), size = 1) +
  geom_linerange(aes(ymin = lower95, ymax = upper95),
                 position = position_dodge(width = .8), size = .5) +
  geom_linerange(data = d_pred_actual_data, aes(ymin = lower80, ymax = upper80),
                 size = 1) +
  geom_linerange(data = d_pred_actual_data, aes(ymin = lower95, ymax = upper95),
                 size = .5) +
  scale_fill_manual(values = c("grey30", "#5158BB", "#F26DF9", "#EB4B98"),
                    name = "Land-use change") +
  xlab("Climate change") +
  ylab(expression(Mean~italic(R^2))) +
  coord_cartesian(clip = "off",
                  ylim = c(0, max(c(d_trend_pred_agg$upper95,
                                    d_pred_actual_data$upper95), na.rm = T))) +
  theme(axis.title.x = element_text(hjust = 2.5 / 6.6),
        legend.text = element_text(margin = margin(t = 5, b = 5, unit = "pt")),
        legend.key.height = unit(30, "pt"))



# SENSITIVITY AN.: critical species ############################################
################################################################################.

# define critical species
d_critical <- c(
  # butterflies
  "Cacyreus marshalli" = "(re)introduced",
  "Colias crocea" = "migratory", 
  "Vanessa atalanta" = "migratory", 
  "Lampides boeticus" = "migratory", 
  "Vanessa cardui" = "migratory",
  "Erebia bubastis" = "taxonomic status",
  # grasshoppers
  "Acheta domesticus" = "(re)introduced",
  "Gryllomorpha dalmatina" = "(re)introduced",
  "Epacromius tergestinus" = "(re)introduced",
  # dragonflies
  "Anax ephippiger" = "migratory",
  "Sympetrum fonscolombii" = "migratory"
) %>% 
  data.frame(critical = .) %>% 
  rownames_to_column("species") 

d_critical <- data.frame(species = sort(unique(l_trends$lm_trends_40$species))) %>% 
  filter(grepl("Adscita", species) |
           grepl("Jordanita", species) |
           species %in% c("Colias alfacariensis", "Colias hyale", "Pyrgus warrenensis") |
           grepl("Leptidea", species)) %>% 
  mutate(critical = "difficult identification") %>% 
  bind_rows(d_critical, .)

# define rare species (few records)
sel_rare <-
  d_records %>% 
  filter(species %in% d_trend_comb$species) %>% 
  group_by(group, species) %>% 
  summarise(n = n()) %>% 
  group_by(group) %>% 
  mutate(nspec = length(n),
         threshold = .2 * nspec) %>% # lowest 20%
  arrange(n) %>% 
  mutate(index = seq_along(n)) %>% 
  filter(index < threshold) %>% 
  ungroup() %>%
  select(species) %>% 
  deframe()


# Modelling --------------------------------------------------------------------.
# ------------------------------------------------------------------------------.

d_trend_comb_z_rare <-
  d_trend_comb_z %>% 
  filter(!species %in% sel_rare) %>% 
  arrange(zone, year_start_f, species) %>% 
  droplevels()

# arrange data
l_data_rare <- list(
  N_obs = nrow(d_trend_comb_z_rare),
  N_obs_agrspec = nrow(d_trend_comb_z_rare[d_trend_comb_z_rare$species %in% sel_agricultural_species, ]),
  N_spec = nlevels(d_trend_comb_z_rare$species),
  N_agrspec = n_distinct(d_trend_comb_z_rare$species[d_trend_comb_z_rare$species %in% sel_agricultural_species]),
  N_reg = nlevels(d_trend_comb_z_rare$zone),
  N_int= nlevels(d_trend_comb_z_rare$year_start_f),
  N_group = nlevels(d_trend_comb_z_rare$group),
  
  group = as.numeric(d_trend_comb_z_rare$group),
  spec = as.numeric(d_trend_comb_z_rare$species),
  spec_agrspec = as.numeric(droplevels(d_trend_comb_z_rare$species[d_trend_comb_z_rare$species %in% sel_agricultural_species])),
  reg = as.numeric(d_trend_comb_z_rare$zone),
  
  T_all = d_trend_comb_z_rare$T_all_change_mean,
  T_seas = d_trend_comb_z_rare$T_seas_change_mean,
  P = d_trend_comb_z_rare$P_warmest_quarter_change_mean,
  AgrArea_prop = d_trend_comb_z_rare$AgrArea_prop_change_mean[d_trend_comb_z_rare$species %in% sel_agricultural_species],
  LSUGL = d_trend_comb_z_rare$LSUGL_change_mean[d_trend_comb_z_rare$species %in% sel_agricultural_species],
  IAR = d_trend_comb_z_rare$IAR_change_mean,
  
  Tniche = d_trend_comb_z_rare$Tniche_mean,
  spec_index = d_trend_comb_z_rare$spec_index,
  
  elevation = as.numeric(d_trend_comb_z_rare$elevation == "high"),
  
  year_start = model.matrix(~ year_start_f, data = d_trend_comb_z_rare)[, -1],
  
  index_agrspec = which(d_trend_comb_z_rare$species %in% sel_agricultural_species),
  
  y =  d_trend_comb_z_rare$trend_mean_sqrt,
  
  y_mean = mean(d_trend_comb_z_rare$trend_mean_sqrt)
)



# run model --------------------------------------------------------------------.
set.seed(73)
fit_rare <- stan(file = "Stan_Code/Stan_regression_restricted.stan.stan",
                 data = l_data_rare,
                 chains = 4, iter = 2000)

# Tabular model results --------------------------------------------------------.
# ------------------------------------------------------------------------------.

l_fit <- extract(fit_rare, pars = c("mu_a", "
                                     b_T_all", "b_T_seas", "b_P",  
                                    "b_AgrArea_prop", "b_LSUGL", "b_IAR", 
                                    "b_Tniche", "b_spec_index", 
                                    "b_elevation",
                                    "b_T_all_Tniche", "b_T_seas_Tniche", "b_P_Tniche",
                                    "b_AgrArea_prop_spec_index", "b_LSUGL_spec_index", "b_IAR_spec_index",
                                    "b_T_all_AgrArea_prop", "b_T_seas_AgrArea_prop", "b_P_AgrArea_prop", 
                                    "b_T_all_LSUGL", "b_T_seas_LSUGL", "b_P_LSUGL",
                                    "b_T_all_IAR", "b_T_seas_IAR", "b_P_IAR",
                                    "b_Tniche_elevation", "b_spec_index_elevation", 
                                    "b_T_all_elevation", "b_T_seas_elevation", "b_P_elevation", 
                                    "b_AgrArea_prop_elevation", "b_LSUGL_elevation",
                                    "b_int"))


d_hdi <- data.frame()


# climate main -----------------------------------------------------------------.

for (cl_var_i in c("T_all", "T_seas", "P")) {
  
  hdi_target <- hdi(l_fit[[paste0("b_", cl_var_i)]], ci = c(.8, .95))
  d_hdi <- data.frame(mean = mean(l_fit[[paste0("b_", cl_var_i)]]),
                      lower95 = hdi_target$CI_low[2],
                      lower80 = hdi_target$CI_low[1],
                      upper80 = hdi_target$CI_high[1],
                      upper95 = hdi_target$CI_high[2]) %>% 
    mutate(cl_var = cl_var_i, elevation = "low", category = "climate_main") %>% 
    bind_rows(d_hdi, .)
  
  hdi_target <- hdi(l_fit[[paste0("b_", cl_var_i)]] + l_fit[[paste0("b_", cl_var_i, "_elevation")]], ci = c(.8, .95))
  d_hdi <- data.frame(mean = mean(l_fit[[paste0("b_", cl_var_i)]] + l_fit[[paste0("b_", cl_var_i, "_elevation")]]),
                      lower95 = hdi_target$CI_low[2],
                      lower80 = hdi_target$CI_low[1],
                      upper80 = hdi_target$CI_high[1],
                      upper95 = hdi_target$CI_high[2]) %>% 
    mutate(cl_var = cl_var_i, elevation = "high", category = "climate_main") %>% 
    bind_rows(d_hdi, .)
  
}


# land use main ----------------------------------------------------------------.

for (lu_var_i in c("AgrArea_prop", "LSUGL")) {
  
  hdi_target <- hdi(l_fit[[paste0("b_", lu_var_i)]], ci = c(.8, .95))
  d_hdi <- data.frame(mean = mean(l_fit[[paste0("b_", lu_var_i)]]), 
                      lower95 = hdi_target$CI_low[2],
                      lower80 = hdi_target$CI_low[1],
                      upper80 = hdi_target$CI_high[1],
                      upper95 = hdi_target$CI_high[2]) %>% 
    mutate(lu_var = lu_var_i, elevation = "low", category = "land_use_main") %>% 
    bind_rows(d_hdi, .)
  
  hdi_target <- hdi(l_fit[[paste0("b_", lu_var_i)]] + l_fit[[paste0("b_", lu_var_i, "_elevation")]], ci = c(.8, .95))
  d_hdi <- data.frame(mean = mean(l_fit[[paste0("b_", lu_var_i)]] + l_fit[[paste0("b_", lu_var_i, "_elevation")]]), 
                      lower95 = hdi_target$CI_low[2],
                      lower80 = hdi_target$CI_low[1],
                      upper80 = hdi_target$CI_high[1],
                      upper95 = hdi_target$CI_high[2]) %>% 
    mutate(lu_var = lu_var_i, elevation = "high", category = "land_use_main") %>% 
    bind_rows(d_hdi, .)
  
}

# no altitude interaction for IAR
hdi_target <- hdi(l_fit[["b_IAR"]], ci = c(.8, .95))
d_hdi <- data.frame(mean = mean(l_fit[["b_IAR"]]),
                    lower95 = hdi_target$CI_low[2],
                    lower80 = hdi_target$CI_low[1],
                    upper80 = hdi_target$CI_high[1],
                    upper95 = hdi_target$CI_high[2]) %>% 
  mutate(lu_var = "IAR", elevation = NA, category = "land_use_main") %>% 
  bind_rows(d_hdi, .)


# climate x land use -----------------------------------------------------------.

for (cl_var_i in c("T_all", "T_seas", "P")) {
  for (lu_var_i in c("AgrArea_prop", "LSUGL", "IAR")) {
    
    hdi_target <- hdi(l_fit[[paste0("b_", cl_var_i, "_", lu_var_i)]], ci = c(.8, .95))
    d_hdi <- data.frame(mean = mean(l_fit[[paste0("b_", cl_var_i, "_", lu_var_i)]]), 
                        lower95 = hdi_target$CI_low[2],
                        lower80 = hdi_target$CI_low[1],
                        upper80 = hdi_target$CI_high[1],
                        upper95 = hdi_target$CI_high[2]) %>%
      mutate(cl_var = cl_var_i, lu_var = lu_var_i, category = "climate_land_use") %>% 
      bind_rows(d_hdi, .)
    
  }
}

# T niche main -----------------------------------------------------------------.

hdi_target <- hdi(l_fit$b_Tniche, ci = c(.8, .95))
d_hdi <- data.frame(mean = mean(l_fit$b_Tniche), 
                    lower95 = hdi_target$CI_low[2],
                    lower80 = hdi_target$CI_low[1],
                    upper80 = hdi_target$CI_high[1],
                    upper95 = hdi_target$CI_high[2]) %>%
  mutate(trT_var = "Tniche",  elevation = "low", category = "Tniche_main") %>% 
  bind_rows(d_hdi, .)

hdi_target <- hdi(l_fit$b_Tniche + l_fit$b_Tniche_elevation, ci = c(.8, .95))
d_hdi <- data.frame(mean = mean(l_fit$b_Tniche + l_fit$b_Tniche_elevation), 
                    lower95 = hdi_target$CI_low[2],
                    lower80 = hdi_target$CI_low[1],
                    upper80 = hdi_target$CI_high[1],
                    upper95 = hdi_target$CI_high[2]) %>%
  mutate(trT_var = "Tniche",  elevation = "high", category = "Tniche_main") %>% 
  bind_rows(d_hdi, .)


# T niche x climate ------------------------------------------------------------.

for (cl_var_i in c("T_all", "T_seas", "P")) {
  
  hdi_target <- hdi(l_fit[[paste0("b_", cl_var_i, "_Tniche")]], ci = c(.8, .95))
  d_hdi <- data.frame(mean = mean(l_fit[[paste0("b_", cl_var_i, "_Tniche")]]), 
                      lower95 = hdi_target$CI_low[2],
                      lower80 = hdi_target$CI_low[1],
                      upper80 = hdi_target$CI_high[1],
                      upper95 = hdi_target$CI_high[2]) %>%
    mutate(cl_var = cl_var_i, trT_var = "Tniche", category = "climate_Tniche") %>% 
    bind_rows(d_hdi, .)
  
}

# specialisation ---------------------------------------------------------------.

hdi_target <- hdi(l_fit$b_spec_index, ci = c(.8, .95))
d_hdi <- data.frame(mean = mean(l_fit$b_spec_index), 
                    lower95 = hdi_target$CI_low[2],
                    lower80 = hdi_target$CI_low[1],
                    upper80 = hdi_target$CI_high[1],
                    upper95 = hdi_target$CI_high[2]) %>%
  mutate(spec_var = "spec_index", elevation = "low", category = "spec_main") %>% 
  bind_rows(d_hdi, .)

hdi_target <- hdi(l_fit$b_spec_index + l_fit$b_spec_index_elevation, ci = c(.8, .95))
d_hdi <- data.frame(mean = mean(l_fit$b_spec_index + l_fit$b_spec_index_elevation), 
                    lower95 = hdi_target$CI_low[2],
                    lower80 = hdi_target$CI_low[1],
                    upper80 = hdi_target$CI_high[1],
                    upper95 = hdi_target$CI_high[2]) %>%
  mutate(spec_var = "spec_index", elevation = "high", category = "spec_main") %>% 
  bind_rows(d_hdi, .)

# land use x specialisation ----------------------------------------------------.

for (lu_var_i in c("AgrArea_prop", "LSUGL", "IAR")) {
  
  hdi_target <- hdi(l_fit[[paste0("b_", lu_var_i, "_spec_index")]], ci = c(.8, .95))
  d_hdi <- data.frame(mean = mean(l_fit[[paste0("b_", lu_var_i, "_spec_index")]]), 
                      lower95 = hdi_target$CI_low[2],
                      lower80 = hdi_target$CI_low[1],
                      upper80 = hdi_target$CI_high[1],
                      upper95 = hdi_target$CI_high[2]) %>%
    mutate(lu_var = lu_var_i, spec_var = "spec_index", category = "land_use_spec") %>% 
    bind_rows(d_hdi, .)
}

# elevation -----------------------------------------------------------------------.

hdi_target <- hdi(l_fit$b_elevation, ci = c(.8, .95))
d_hdi <- data.frame(mean = mean(l_fit$b_elevation), 
                    lower95 = hdi_target$CI_low[2],
                    lower80 = hdi_target$CI_low[1],
                    upper80 = hdi_target$CI_high[1],
                    upper95 = hdi_target$CI_high[2]) %>%
  mutate(category = "elevation_main") %>% 
  bind_rows(d_hdi, .)

# group intercepts -------------------------------------------------------------.

for (gr_i in 1:nlevels(d_trend_comb_z_rare2$group)) {
  hdi_target <- hdi(l_fit$mu_a[, gr_i], ci = c(.8, .95))
  d_hdi <- data.frame(mean = mean(l_fit$mu_a[, gr_i]),
                      lower95 = hdi_target$CI_low[2],
                      lower80 = hdi_target$CI_low[1],
                      upper80 = hdi_target$CI_high[1],
                      upper95 = hdi_target$CI_high[2]) %>% 
    mutate(category = "group_intercept", group = paste0("Intercept (", 
                                                        levels(d_trend_comb_z_rare2$group)[gr_i],
                                                        ")")) %>% 
    bind_rows(d_hdi, .)
}

# time interval ----------------------------------------------------------------.

v_intervals <- c("1985-1990", "1990-1995",
                 "1995-2000", "2000-2005", 
                 "2005-2010", "2010-2015",
                 "2015-2020")

for (int_i in 1:7) {
  hdi_target <- hdi(l_fit$b_int[, int_i], ci = c(.8, .95))
  d_hdi <- data.frame(mean = mean(l_fit$b_int[, int_i]),
                      lower95 = hdi_target$CI_low[2],
                      lower80 = hdi_target$CI_low[1],
                      upper80 = hdi_target$CI_high[1],
                      upper95 = hdi_target$CI_high[2]) %>% 
    mutate(category = "time_interval", int = paste("Interval", v_intervals[int_i])) %>% 
    bind_rows(d_hdi, .)
}

# final edits ------------------------------------------------------------------.

d_hdi <-
  d_hdi %>% 
  mutate(cl_var = factor(cl_var,
                         levels = c("T_all",
                                    "T_seas",
                                    "P"),
                         labels = c("T. mean",
                                    "T. seasonality",
                                    "P. summer")),
         lu_var = factor(lu_var,
                         levels = c("AgrArea_prop",
                                    "LSUGL",
                                    "IAR"),
                         labels = c("Agr. area",
                                    "Grassland-use int.",
                                    "Crop-use int.")),
         trT_var = factor(trT_var,
                          levels = c("Tniche"),
                          labels = c("Temp. niche")),
         spec_var = factor(spec_var,
                           levels = c("spec_index"),
                           labels = c("Specialisation")),
         int = factor(int),
         group = factor(group),
         mean = mean * scaling_parameters$trend_mean_sqrt$sd,
         lower95 = lower95 * scaling_parameters$trend_mean_sqrt$sd,
         lower80 = lower80 * scaling_parameters$trend_mean_sqrt$sd,
         upper80 = upper80 * scaling_parameters$trend_mean_sqrt$sd,
         upper95 = upper95 * scaling_parameters$trend_mean_sqrt$sd) %>% # change to original level
  rowwise() %>% 
  mutate(var = ifelse(is.na(elevation), paste(c(cl_var, lu_var, trT_var, spec_var, int, group), collapse = " × "),
                      paste0(paste(c(cl_var, lu_var, trT_var, spec_var, int, group), collapse = " × "), " (", elevation, ")")),
         var = gsub(" × NA", "", var),
         var = gsub("NA × ", "", var),
         var= ifelse(category == "elevation_main", "Elevation (high)", var)) %>% 
  ungroup()

d_hdi %>%
  select(var, lower95, lower80, mean, upper80, upper95)

# SENSITIVITY AN.: Rare species-region combinations ############################
################################################################################.

sel_rare2 <-
  d_records %>% 
  filter(species %in% d_trend_comb$species,
         A > 1979 & A <= 2020) %>% 
  left_join(d_zone, by = "km2_ID") %>% 
  group_by(group, species, zone) %>% 
  summarise(n = n(),
            .groups = "drop") %>% 
  filter(n < 41) %>% # average of 1 observation for every year
  mutate(ID = paste(species, zone))

# Modelling --------------------------------------------------------------------.
# ------------------------------------------------------------------------------.

d_trend_comb_z_rare2 <-
  d_trend_comb_z %>% 
  filter(!paste(species, zone) %in% sel_rare2$ID) %>% 
  arrange(zone, year_start_f, species) %>% 
  droplevels()

# arrange data
l_data_rare2 <- list(
  N_obs = nrow(d_trend_comb_z_rare2),
  N_obs_agrspec = nrow(d_trend_comb_z_rare2[d_trend_comb_z_rare2$species %in% sel_agricultural_species, ]),
  N_spec = nlevels(d_trend_comb_z_rare2$species),
  N_agrspec = n_distinct(d_trend_comb_z_rare2$species[d_trend_comb_z_rare2$species %in% sel_agricultural_species]),
  N_reg = nlevels(d_trend_comb_z_rare2$zone),
  N_int= nlevels(d_trend_comb_z_rare2$year_start_f),
  N_group = nlevels(d_trend_comb_z_rare2$group),
  
  group = as.numeric(d_trend_comb_z_rare2$group),
  spec = as.numeric(d_trend_comb_z_rare2$species),
  spec_agrspec = as.numeric(droplevels(d_trend_comb_z_rare2$species[d_trend_comb_z_rare2$species %in% sel_agricultural_species])),
  reg = as.numeric(d_trend_comb_z_rare2$zone),
  
  T_all = d_trend_comb_z_rare2$T_all_change_mean,
  T_seas = d_trend_comb_z_rare2$T_seas_change_mean,
  P = d_trend_comb_z_rare2$P_warmest_quarter_change_mean,
  Agrare2a_prop = d_trend_comb_z_rare2$Agrare2a_prop_change_mean[d_trend_comb_z_rare2$species %in% sel_agricultural_species],
  LSUGL = d_trend_comb_z_rare2$LSUGL_change_mean[d_trend_comb_z_rare2$species %in% sel_agricultural_species],
  IAR = d_trend_comb_z_rare2$IAR_change_mean,
  
  Tniche = d_trend_comb_z_rare2$Tniche_mean,
  spec_index = d_trend_comb_z_rare2$spec_index,
  
  elevation = as.numeric(d_trend_comb_z_rare2$elevation == "high"),
  
  year_start = model.matrix(~ year_start_f, data = d_trend_comb_z_rare2)[, -1],
  
  index_agrspec = which(d_trend_comb_z_rare2$species %in% sel_agricultural_species),
  
  y =  d_trend_comb_z_rare2$trend_mean_sqrt,
  
  y_mean = mean(d_trend_comb_z_rare2$trend_mean_sqrt)
)



# run model --------------------------------------------------------------------.
set.seed(213)
fit_rare2 <- stan(file = "Stan_Code/Stan_regression_restricted.stan.stan",
                  data = l_data_rare2,
                  chains = 4, iter = 2000)


# Tabular model results --------------------------------------------------------.
# ------------------------------------------------------------------------------.

l_fit <- extract(fit_rare2, pars = c("mu_a", "
                                     b_T_all", "b_T_seas", "b_P",  
                                     "b_AgrArea_prop", "b_LSUGL", "b_IAR", 
                                     "b_Tniche", "b_spec_index", 
                                     "b_elevation",
                                     "b_T_all_Tniche", "b_T_seas_Tniche", "b_P_Tniche",
                                     "b_AgrArea_prop_spec_index", "b_LSUGL_spec_index", "b_IAR_spec_index",
                                     "b_T_all_AgrArea_prop", "b_T_seas_AgrArea_prop", "b_P_AgrArea_prop", 
                                     "b_T_all_LSUGL", "b_T_seas_LSUGL", "b_P_LSUGL",
                                     "b_T_all_IAR", "b_T_seas_IAR", "b_P_IAR",
                                     "b_Tniche_elevation", "b_spec_index_elevation", 
                                     "b_T_all_elevation", "b_T_seas_elevation", "b_P_elevation", 
                                     "b_AgrArea_prop_elevation", "b_LSUGL_elevation",
                                     "b_int"))


d_hdi <- data.frame()

# climate main -----------------------------------------------------------------.

for (cl_var_i in c("T_all", "T_seas", "P")) {
  
  hdi_target <- hdi(l_fit[[paste0("b_", cl_var_i)]], ci = c(.8, .95))
  d_hdi <- data.frame(mean = mean(l_fit[[paste0("b_", cl_var_i)]]),
                      lower95 = hdi_target$CI_low[2],
                      lower80 = hdi_target$CI_low[1],
                      upper80 = hdi_target$CI_high[1],
                      upper95 = hdi_target$CI_high[2]) %>% 
    mutate(cl_var = cl_var_i, elevation = "low", category = "climate_main") %>% 
    bind_rows(d_hdi, .)
  
  hdi_target <- hdi(l_fit[[paste0("b_", cl_var_i)]] + l_fit[[paste0("b_", cl_var_i, "_elevation")]], ci = c(.8, .95))
  d_hdi <- data.frame(mean = mean(l_fit[[paste0("b_", cl_var_i)]] + l_fit[[paste0("b_", cl_var_i, "_elevation")]]),
                      lower95 = hdi_target$CI_low[2],
                      lower80 = hdi_target$CI_low[1],
                      upper80 = hdi_target$CI_high[1],
                      upper95 = hdi_target$CI_high[2]) %>% 
    mutate(cl_var = cl_var_i, elevation = "high", category = "climate_main") %>% 
    bind_rows(d_hdi, .)
  
}


# land use main ----------------------------------------------------------------.

for (lu_var_i in c("AgrArea_prop", "LSUGL")) {
  
  hdi_target <- hdi(l_fit[[paste0("b_", lu_var_i)]], ci = c(.8, .95))
  d_hdi <- data.frame(mean = mean(l_fit[[paste0("b_", lu_var_i)]]), 
                      lower95 = hdi_target$CI_low[2],
                      lower80 = hdi_target$CI_low[1],
                      upper80 = hdi_target$CI_high[1],
                      upper95 = hdi_target$CI_high[2]) %>% 
    mutate(lu_var = lu_var_i, elevation = "low", category = "land_use_main") %>% 
    bind_rows(d_hdi, .)
  
  hdi_target <- hdi(l_fit[[paste0("b_", lu_var_i)]] + l_fit[[paste0("b_", lu_var_i, "_elevation")]], ci = c(.8, .95))
  d_hdi <- data.frame(mean = mean(l_fit[[paste0("b_", lu_var_i)]] + l_fit[[paste0("b_", lu_var_i, "_elevation")]]), 
                      lower95 = hdi_target$CI_low[2],
                      lower80 = hdi_target$CI_low[1],
                      upper80 = hdi_target$CI_high[1],
                      upper95 = hdi_target$CI_high[2]) %>% 
    mutate(lu_var = lu_var_i, elevation = "high", category = "land_use_main") %>% 
    bind_rows(d_hdi, .)
  
}

# no altitude interaction for IAR
hdi_target <- hdi(l_fit[["b_IAR"]], ci = c(.8, .95))
d_hdi <- data.frame(mean = mean(l_fit[["b_IAR"]]),
                    lower95 = hdi_target$CI_low[2],
                    lower80 = hdi_target$CI_low[1],
                    upper80 = hdi_target$CI_high[1],
                    upper95 = hdi_target$CI_high[2]) %>% 
  mutate(lu_var = "IAR", elevation = NA, category = "land_use_main") %>% 
  bind_rows(d_hdi, .)


# climate x land use -----------------------------------------------------------.

for (cl_var_i in c("T_all", "T_seas", "P")) {
  for (lu_var_i in c("AgrArea_prop", "LSUGL", "IAR")) {
    
    hdi_target <- hdi(l_fit[[paste0("b_", cl_var_i, "_", lu_var_i)]], ci = c(.8, .95))
    d_hdi <- data.frame(mean = mean(l_fit[[paste0("b_", cl_var_i, "_", lu_var_i)]]), 
                        lower95 = hdi_target$CI_low[2],
                        lower80 = hdi_target$CI_low[1],
                        upper80 = hdi_target$CI_high[1],
                        upper95 = hdi_target$CI_high[2]) %>%
      mutate(cl_var = cl_var_i, lu_var = lu_var_i, category = "climate_land_use") %>% 
      bind_rows(d_hdi, .)
    
  }
}

# T niche main -----------------------------------------------------------------.

hdi_target <- hdi(l_fit$b_Tniche, ci = c(.8, .95))
d_hdi <- data.frame(mean = mean(l_fit$b_Tniche), 
                    lower95 = hdi_target$CI_low[2],
                    lower80 = hdi_target$CI_low[1],
                    upper80 = hdi_target$CI_high[1],
                    upper95 = hdi_target$CI_high[2]) %>%
  mutate(trT_var = "Tniche",  elevation = "low", category = "Tniche_main") %>% 
  bind_rows(d_hdi, .)

hdi_target <- hdi(l_fit$b_Tniche + l_fit$b_Tniche_elevation, ci = c(.8, .95))
d_hdi <- data.frame(mean = mean(l_fit$b_Tniche + l_fit$b_Tniche_elevation), 
                    lower95 = hdi_target$CI_low[2],
                    lower80 = hdi_target$CI_low[1],
                    upper80 = hdi_target$CI_high[1],
                    upper95 = hdi_target$CI_high[2]) %>%
  mutate(trT_var = "Tniche",  elevation = "high", category = "Tniche_main") %>% 
  bind_rows(d_hdi, .)


# T niche x climate ------------------------------------------------------------.

for (cl_var_i in c("T_all", "T_seas", "P")) {
  
  hdi_target <- hdi(l_fit[[paste0("b_", cl_var_i, "_Tniche")]], ci = c(.8, .95))
  d_hdi <- data.frame(mean = mean(l_fit[[paste0("b_", cl_var_i, "_Tniche")]]), 
                      lower95 = hdi_target$CI_low[2],
                      lower80 = hdi_target$CI_low[1],
                      upper80 = hdi_target$CI_high[1],
                      upper95 = hdi_target$CI_high[2]) %>%
    mutate(cl_var = cl_var_i, trT_var = "Tniche", category = "climate_Tniche") %>% 
    bind_rows(d_hdi, .)
  
}

# specialisation ---------------------------------------------------------------.

hdi_target <- hdi(l_fit$b_spec_index, ci = c(.8, .95))
d_hdi <- data.frame(mean = mean(l_fit$b_spec_index), 
                    lower95 = hdi_target$CI_low[2],
                    lower80 = hdi_target$CI_low[1],
                    upper80 = hdi_target$CI_high[1],
                    upper95 = hdi_target$CI_high[2]) %>%
  mutate(spec_var = "spec_index", elevation = "low", category = "spec_main") %>% 
  bind_rows(d_hdi, .)

hdi_target <- hdi(l_fit$b_spec_index + l_fit$b_spec_index_elevation, ci = c(.8, .95))
d_hdi <- data.frame(mean = mean(l_fit$b_spec_index + l_fit$b_spec_index_elevation), 
                    lower95 = hdi_target$CI_low[2],
                    lower80 = hdi_target$CI_low[1],
                    upper80 = hdi_target$CI_high[1],
                    upper95 = hdi_target$CI_high[2]) %>%
  mutate(spec_var = "spec_index", elevation = "high", category = "spec_main") %>% 
  bind_rows(d_hdi, .)

# land use x specialisation ----------------------------------------------------.

for (lu_var_i in c("AgrArea_prop", "LSUGL", "IAR")) {
  
  hdi_target <- hdi(l_fit[[paste0("b_", lu_var_i, "_spec_index")]], ci = c(.8, .95))
  d_hdi <- data.frame(mean = mean(l_fit[[paste0("b_", lu_var_i, "_spec_index")]]), 
                      lower95 = hdi_target$CI_low[2],
                      lower80 = hdi_target$CI_low[1],
                      upper80 = hdi_target$CI_high[1],
                      upper95 = hdi_target$CI_high[2]) %>%
    mutate(lu_var = lu_var_i, spec_var = "spec_index", category = "land_use_spec") %>% 
    bind_rows(d_hdi, .)
}

# elevation --------------------------------------------------------------------.

hdi_target <- hdi(l_fit$b_elevation, ci = c(.8, .95))
d_hdi <- data.frame(mean = mean(l_fit$b_elevation), 
                    lower95 = hdi_target$CI_low[2],
                    lower80 = hdi_target$CI_low[1],
                    upper80 = hdi_target$CI_high[1],
                    upper95 = hdi_target$CI_high[2]) %>%
  mutate(category = "elevation_main") %>% 
  bind_rows(d_hdi, .)

# group intercepts -------------------------------------------------------------.

for (gr_i in 1:nlevels(d_trend_comb_z_rare2$group)) {
  hdi_target <- hdi(l_fit$mu_a[, gr_i], ci = c(.8, .95))
  d_hdi <- data.frame(mean = mean(l_fit$mu_a[, gr_i]),
                      lower95 = hdi_target$CI_low[2],
                      lower80 = hdi_target$CI_low[1],
                      upper80 = hdi_target$CI_high[1],
                      upper95 = hdi_target$CI_high[2]) %>% 
    mutate(category = "group_intercept", group = paste0("Intercept (", 
                                                        levels(d_trend_comb_z_rare2$group)[gr_i],
                                                        ")")) %>% 
    bind_rows(d_hdi, .)
}

# time interval ----------------------------------------------------------------.

v_intervals <- c("1985-1990", "1990-1995",
                 "1995-2000", "2000-2005", 
                 "2005-2010", "2010-2015",
                 "2015-2020")

for (int_i in 1:7) {
  hdi_target <- hdi(l_fit$b_int[, int_i], ci = c(.8, .95))
  d_hdi <- data.frame(mean = mean(l_fit$b_int[, int_i]),
                      lower95 = hdi_target$CI_low[2],
                      lower80 = hdi_target$CI_low[1],
                      upper80 = hdi_target$CI_high[1],
                      upper95 = hdi_target$CI_high[2]) %>% 
    mutate(category = "time_interval", int = paste("Interval", v_intervals[int_i])) %>% 
    bind_rows(d_hdi, .)
}

# final edits ------------------------------------------------------------------.


d_hdi <-
  d_hdi %>% 
  mutate(cl_var = factor(cl_var,
                         levels = c("T_all",
                                    "T_seas",
                                    "P"),
                         labels = c("T. mean",
                                    "T. seasonality",
                                    "P. summer")),
         lu_var = factor(lu_var,
                         levels = c("AgrArea_prop",
                                    "LSUGL",
                                    "IAR"),
                         labels = c("Agr. area",
                                    "Grassland-use int.",
                                    "Crop-use int.")),
         trT_var = factor(trT_var,
                          levels = c("Tniche"),
                          labels = c("Temp. niche")),
         spec_var = factor(spec_var,
                           levels = c("spec_index"),
                           labels = c("Specialisation")),
         int = factor(int),
         group = factor(group),
         mean = mean * scaling_parameters$trend_mean_sqrt$sd,
         lower95 = lower95 * scaling_parameters$trend_mean_sqrt$sd,
         lower80 = lower80 * scaling_parameters$trend_mean_sqrt$sd,
         upper80 = upper80 * scaling_parameters$trend_mean_sqrt$sd,
         upper95 = upper95 * scaling_parameters$trend_mean_sqrt$sd) %>% # change to original level
  rowwise() %>% 
  mutate(var = ifelse(is.na(elevation), paste(c(cl_var, lu_var, trT_var, spec_var, int, group), collapse = " × "),
                      paste0(paste(c(cl_var, lu_var, trT_var, spec_var, int, group), collapse = " × "), " (", elevation, ")")),
         var = gsub(" × NA", "", var),
         var = gsub("NA × ", "", var),
         var = ifelse(category == "elevation_main", "Elevation (high)", var)) %>% 
  ungroup()

d_hdi %>%
  select(var, lower95, lower80, mean, upper80, upper95) 

