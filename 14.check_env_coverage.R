#--------------------------------------------------------------------------#
# Project: planktyp
# Script purpose: Check environmental coverage
# Date: 18/08/2022
# Author: Thelma Panaïotis
#--------------------------------------------------------------------------#

library(tidyverse)
library(castr)
library(lubridate)
library(tidync)
library(ggpmisc)

load("data/06.all_data.Rdata")
load("data/03.layers.Rdata")

all_data %>% glimpse()

all_data %>% 
  select(psampleid, lat, lon, datetime, temp, sal, oxy_ml, chla, layer)

layers %>% select(psampleid, depth_epi = epi)

## Read WOA data ----
#--------------------------------------------------------------------------#
file <- "data/raw/woa/woa18_A5B7_t00_01.nc"
temp <- tidync(file) %>% 
  hyper_filter(depth = depth <= 500) %>% 
  hyper_tibble() %>% 
  select(lat, lon, depth, temp_woa = t_an) %>% 
  unique()

file <- "data/raw/woa/woa18_A5B7_s00_01.nc"
sal <- tidync(file) %>% 
  hyper_filter(depth = depth <= 500) %>% 
  hyper_tibble() %>% 
  select(lat, lon, depth, sal_woa = s_an) %>% 
  unique()

file <- "data/raw/woa/woa18_all_o00_01.nc"
oxy <- tidync(file) %>% 
  hyper_filter(depth = depth <= 500) %>% 
  hyper_tibble() %>% 
  select(lat, lon, depth, oxy_woa = o_an) %>% 
  unique()

file <- "data/raw/woa/woa18_all_A00_01.nc"
aou <- tidync(file) %>% 
  hyper_filter(depth = depth <= 500) %>% 
  hyper_tibble() %>% 
  select(lat, lon, depth, aou_woa = A_an) %>% 
  unique()


# Store WOA variables together
woa <- left_join(temp, sal) %>% left_join(oxy) %>% left_join(aou)


## Prepare UVP data to match with WOA ----
#--------------------------------------------------------------------------#
# Round station coordinates at 1°
# Keep depth of the epipelagic layer to average WOA data over epi and meso layers
coord <- layers %>% 
  select(lat, lon, epi) %>% 
  mutate(
    lat = roundp(lat, precision = 1, f = floor) + 0.5,
    lon = roundp(lon, precision = 1, f = floor) + 0.5,
  )


## Join WOA with UVP data ----
#--------------------------------------------------------------------------#
woa_uvp <- woa %>% 
  left_join(coord) %>% 
  drop_na(epi) # Drop locations where there are no UVP data

woa_uvp

## Average WOA values over epi and mesosup layers ----
#--------------------------------------------------------------------------#
# Average over epi layer
woa_epi <- woa_uvp %>% 
  filter(depth <= epi) %>% 
  group_by(lat, lon) %>% 
  summarise(
    temp_woa = mean(temp_woa, na.rm = TRUE),
    sal_woa  = mean(sal_woa, na.rm = TRUE),
    oxy_woa  = mean(oxy_woa, na.rm = TRUE),
    aou_woa  = mean(aou_woa, na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  mutate(layer = "epi")

# Average over meso layer
woa_meso <- woa_uvp %>% 
  filter(depth > epi) %>% 
  group_by(lat, lon) %>% 
  summarise(
    temp_woa = mean(temp_woa, na.rm = TRUE),
    sal_woa  = mean(sal_woa, na.rm = TRUE),
    oxy_woa  = mean(oxy_woa, na.rm = TRUE),
    aou_woa  = mean(aou_woa, na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  mutate(layer = "mesosup")

# Plot a map
#woa_epi %>% 
#  ggplot() +
#  geom_point(aes(x = lon, y = lat, color = temp_woa)) +
#  scale_color_viridis_c() +
#  coord_quickmap()

# Store epi and meso data together
woa_epi_meso <- bind_rows(woa_epi, woa_meso)



## Correlation between in situ and WOA at UVP stations ----
#--------------------------------------------------------------------------#
# Join in situ and WOA data at UVP stations
corr_data <- all_data %>% 
  filter(layer %in% c("epi", "mesosup")) %>% 
  select(lat, lon, layer, temp_uvp = temp, sal_uvp = sal, oxy_uvp = oxy_ml, aou_uvp = aou) %>% 
  mutate(
    lat = roundp(lat, precision = 1, f = floor) + 0.5,
    lon = roundp(lon, precision = 1, f = floor) + 0.5,
    station = row_number(), .before = lat
  ) %>% 
  left_join(woa_epi_meso) 

# Nice formatting for variable names
var_names <- tibble(
  variable = c("temp", "sal", "oxy", "aou"),
  nice_name = c("Temperature", "Salinity", "Oxygen", "AOU")
) %>% 
  mutate(nice_name = factor(nice_name, levels = c("Temperature", "Salinity", "Oxygen", "AOU")))

corr_data %>% 
  ggplot() +
  geom_point(aes(x = aou_uvp, y = aou_woa, color = layer))


# Reformat table
corr_data <- corr_data %>% 
  gather(temp_uvp:aou_woa, key = "key", value = "value") %>% 
  separate(key, into = c("variable", "source")) %>% 
  pivot_wider(names_from = source, values_from = value) %>% 
  left_join(var_names)


#corr_data %>% 
#  ggplot(aes(x = uvp, y = woa)) +
#  geom_point( alpha = 0.2) +
#  geom_smooth(method='lm') +
#  facet_grid(layer~nice_name, scales = "free") +
#  theme_minimal() +
#  labs(x = "In situ", y = "WOA") +
#  ggtitle("Correlation between in situ and WOA data at UVP stations in epipelagic and mesopelagic layers")

corr_data %>% 
  ggplot(aes(x = uvp, y = woa, color = layer)) +
  geom_point(alpha = 0.2) +
  #geom_smooth(method='lm', se = FALSE) +
  stat_poly_line(se = FALSE) +
  stat_poly_eq() +
  facet_wrap(~nice_name, scales = "free") +
  scale_color_brewer(palette="Dark2") + 
  theme_classic() +
  #theme_minimal() +
  #ggtitle("Correlation between in situ and WOA data at UVP stations in epipelagic and mesopelagic layers") +
  labs(x = "In situ", y = "WOA", color = "Layer")
ggsave("plots/zoo/14.uvp_woa_corr.png")  
ggsave(file = "plots/paper/14.uvp_woa_corr.png", width = 166, height = 100, unit = "mm", dpi = 300)



## Distribution comparison between all WOA and WOA at UVP stations ----
#--------------------------------------------------------------------------#
woa_layer <- woa %>% 
  mutate(depth_bin = ifelse(depth <= 200, "0-200 m", "200-500 m")) %>% 
  group_by(lon, lat, depth_bin) %>% 
  summarise(
    temp_woa_All = mean(temp_woa, na.rm = TRUE),
    sal_woa_All  = mean(sal_woa, na.rm = TRUE),
    oxy_woa_All  = mean(oxy_woa, na.rm = TRUE),
    aou_woa_All  = mean(aou_woa, na.rm = TRUE)
  ) %>% 
  ungroup()  

woa_uvp_layer <- woa_uvp %>% 
  mutate(depth_bin = ifelse(depth <= 200, "0-200 m", "200-500 m")) %>% 
  group_by(lon, lat, depth_bin) %>% 
  summarise(
    temp_woa_UVP = mean(temp_woa, na.rm = TRUE),
    sal_woa_UVP  = mean(sal_woa, na.rm = TRUE),
    oxy_woa_UVP  = mean(oxy_woa, na.rm = TRUE),
    aou_woa_UVP  = mean(aou_woa, na.rm = TRUE)
  ) %>% 
  ungroup() 

woa_dist <- woa_layer %>% 
  left_join(woa_uvp_layer) %>% 
  gather(temp_woa_All:aou_woa_UVP, key = "key", value = "value") %>% 
  separate(key, into = c("variable", "woa", "location")) %>% 
  select(-woa) %>% 
  left_join(var_names)


woa_dist %>% 
  ggplot() +
  geom_density(aes(x = value, color = location)) +
  scale_color_brewer(palette="Dark2") + 
  facet_grid(depth_bin~nice_name, scales = "free") +
  theme_classic() +
  labs(x = "WOA value", y = "Density", color = "Location") +
  ggtitle("Distribution of WOA environmental data all over the globe and at stations sampled by the UVP")

woa_dist %>% 
  ggplot() +
  geom_density(aes(x = value, color = depth_bin, linetype = location)) +
  scale_color_brewer(palette="Dark2") + 
  facet_wrap(~nice_name, scales = "free") +
  theme_classic() +
  #ggtitle("Distribution of WOA environmental data all over the globe and at stations sampled by the UVP") +
  labs(x = "WOA value", y = "Density", linetype = "Location", color = "Depth")
ggsave("plots/zoo/14.uvp_coverage.png")  
ggsave(file = "plots/paper/14.uvp_coverage.png", width = 166, height = 100, unit = "mm", dpi = 300)



ggplot(df, aes(value)) + 
  stat_ecdf(geom = "step", aes(color = group))

woa_dist %>% 
  ggplot(aes(x = value, color = depth_bin, linetype = location)) +
  #geom_density() +
  stat_ecdf(geom = "step") +
  scale_color_brewer(palette="Dark2") + 
  facet_wrap(~nice_name, scales = "free") +
  theme_classic()

toto <- woa_dist %>% 
  filter(variable == "temp" & depth_bin == "0-200 m") %>% 
  pivot_wider(names_from = location, values_from = value)

titi <- ks.test(na.omit(toto$All), na.omit(toto$UVP))
summary(titi)

woa_dist %>% 
  pivot_wider(names_from = location, values_from = value) %>% 
  group_by(variable, depth_bin) %>% 
  summarise(p_val = ks.test(x = All, y = UVP)$p.val)
