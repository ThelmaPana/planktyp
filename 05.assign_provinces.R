#--------------------------------------------------------------------------#
# Project: planktyo
# Script purpose: Assign profiles to bio-regions (Longhurst provinces, GBP)
# Date: 02/02/2021
# Author: Thelma Panaiotis
#--------------------------------------------------------------------------#

source("lib/set_up.R")
library(sf)
library(castr)

load("data/04.sat.Rdata") # We load satellite data because it only contains profile with complete data

## List profiles ----
#--------------------------------------------------------------------------#
# We read profiles list from satellite data because it contains profiles for which we have all data
my_profiles <- sat %>% select(title, profile, lat, lon, datetime, psampleid)
message(nrow(my_profiles), " profiles with complete data")

# Convert to sf
sf_profiles <- st_as_sf(my_profiles, coords = c("lon", "lat")) %>% st_set_crs(4326)


## Match with Longhurst provinces ----
#--------------------------------------------------------------------------#
longhurst <- read_sf(dsn = "data/raw/regions/longhurst_v4_2010/", layer = "Longhurst_world_v4_2010")
sf_longhurst <- st_as_sf(longhurst)

# Compute intersection of UVP profiles and Longhurst provinces
my_long <- st_join(sf_profiles, sf_longhurst, join = st_intersects) %>% 
  as_tibble() %>% 
  select(title:psampleid, long_code = ProvCode, long = ProvDescr) %>% 
  left_join(my_profiles) %>% # add lat and lon of UVP profiles
  relocate(lat, lon, .after = profile)

message(sum(is.na(my_long$long_code)), " profiles not matched with Longhurst provinces")


## Match with MBGCP (Mesopelagic BioGeoChemical Provinces) ----
#--------------------------------------------------------------------------#
# Read MBGCP data
mbgcp <- read_csv("data/raw/regions/MBGCP.csv") %>% 
  select(nb = MBGCP, lon_m = Lon, lat_m = Lat) %>% 
  filter(!is.nan(nb)) %>% 
  mutate(nb = as.character(nb)) %>% 
  left_join(read_delim("data/raw/regions/MBGCP_desc.csv", 
                       ";", escape_double = FALSE, col_types = cols(nb = col_character()), 
                       trim_ws = TRUE)
            ) %>% 
  select(-nb)

# Compute intersection of UVP profiles and MBGCP
my_mbgcp <- my_profiles %>% 
  mutate( # round lat and lon of UVP profiles to match with MBGCP
    lat_m = roundp(lat, 1, floor) + 0.5,
    lon_m = roundp(lon, 1, floor) + 0.5,
  ) %>% 
  left_join(mbgcp) %>% 
  select(-c(lat_m, lon_m))

message(sum(is.na(my_mbgcp$gbp_code)), " profiles not matched with MBGCP")


## Compute latitude bands ----
#--------------------------------------------------------------------------#

provinces <- my_profiles %>% 
  mutate(b_lat = cut(lat, breaks = seq(-90,90,10))) 


## Join assigned provinces together and save ----
#--------------------------------------------------------------------------#
# Save provinces
provinces <- provinces %>% 
  left_join(my_long) %>% 
  left_join(my_mbgcp)
save(provinces, file = "data/05.provinces.Rdata")


## Plot a few maps ----
#--------------------------------------------------------------------------#
# Latitude bands
provinces %>% 
  ggplot() +
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "gray") +
  geom_point(aes(x = lon, y = lat, colour = b_lat)) +
  coord_quickmap() +
  ggtitle("Latitude bands")
ggsave("plots/provinces/05.map_latitude_bands.png")

# Longhurst provinces
provinces %>% 
  ggplot() +
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "gray") +
  geom_point(aes(x = lon, y = lat, colour = long_code)) +
  coord_quickmap() +
  ggtitle("Longhurst provinces")
ggsave("plots/provinces/05.map_longhurst.png")

# MBGCP 
provinces %>% 
  ggplot() +
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "gray") +
  geom_point(aes(x = lon, y = lat, colour = gbp_code)) +
  coord_quickmap() +
  ggtitle("MBGCP")
ggsave("plots/provinces/05.map_mbgcp.png")


