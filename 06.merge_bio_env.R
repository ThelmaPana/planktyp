#--------------------------------------------------------------------------#
# Project: planktyp
# Script purpose: Merge biological and environmental data
# Date: 02/02/2021
# Author: Thelma Panaiotis
#--------------------------------------------------------------------------#

source("lib/set_up.R")
library(suncalc)

load("data/01.bio_data.Rdata")
load("data/03.ctd_5m.Rdata")
load("data/03.layers.Rdata")
load("data/04.sat.Rdata")
load("data/05.provinces.Rdata")


## Keep only profiles for which we have all data ----
#--------------------------------------------------------------------------#
# sat data contains profiles for which we have all data
my_profiles <- sat %>% select(title, profile, lat, lon, datetime, psampleid)
message(nrow(my_profiles), " profiles with complete data")

# Filter psampleid in zoo_conc
zoo_conc <- zoo_conc %>% filter(psampleid %in% my_profiles$psampleid) %>% select(-sampleid)

# Filter psampleid in detritus
det_conc <- det_conc %>% filter(psampleid %in% my_profiles$psampleid) %>% select(-sampleid)

# Filter psampleid in particles data
part_conc <- part_conc %>% filter(psampleid %in% my_profiles$psampleid) %>% select(-c(sampleid, watervolume))

# Filter psampleid in CTD data
ctd_5m <- ctd_5m %>% filter(psampleid %in% my_profiles$psampleid)

# Filter psampleid in layers data
layers <- layers %>% filter(psampleid %in% my_profiles$psampleid) %>% select(-max_depth)

# Filter psampleid in provinces data
provinces <- provinces %>% filter(psampleid %in% my_profiles$psampleid)


## Assemble env data together ----
#--------------------------------------------------------------------------#
# Pieces to assemble: CTD data, layers, det and parts
# Some profiles might be deeper in det/parts data while others are deeper in ctd data 
# --> perform an inner-join to avoid any problem

env <- ctd_5m %>% # start with ctd data
  inner_join(det_conc %>% rename(snow_conc = detritus)) %>% # add detritus as marine snow
  inner_join(part_conc %>% rename(bulk_conc = conc, bulk_biov = biovol)) %>% # add particles as bulk
  left_join(layers) %>%  # add layers
  left_join(sat) # and satellite data


## Join env data and zooplankton data ----
#--------------------------------------------------------------------------#
# Use an inner-join for same reason as previously
all <- inner_join(zoo_conc, env) %>% 
  relocate(datetime, .after = lon)


## Cut layers ----
#--------------------------------------------------------------------------#
# Cut data into 4 layers:
# - epipelagic layer: surface to epipelagic layer depth (epi)
# - mesopelagic superior (mesosup) layer: epipelagic layer depth (epi) to 500 m
# - mesopelagic inferior (mesoinf) layer: 500 m to 1000 m
# - bathypelagic layer: below 1000 m

# Make sure that profiles cover at least > 80% of theoretical layers.
# This concerns only mesosup and mesoinf layers:
# - if epipelagic layer could not be defined, profile was ignored
# - for bathypelagic, there is no bottom limit

# Epipelagic layer
epi <- all %>% 
  filter(depth < epi) %>% # keep bins shallower than epipelagic layer depth
  group_by(title, profile, lat, lon, datetime, psampleid) %>% # group by profile
  summarise_at(vars(Acantharea:kd490), mean, na.rm = TRUE) %>% # compute mean concentration for all plankton groups and env data
  ungroup() %>% 
  mutate(layer = "epi", .after = psampleid)
message(nrow(epi), " profiles for epipelagic layer")


# Mesopelagic superior layer
mesosup <- all %>% 
  filter(depth > epi & depth < 500) %>% # keep bins between epipelagic layer depth and 500 m
  group_by(title, profile, lat, lon, datetime, psampleid) %>% # group by profile
  mutate( # for each profile compute:
    depth_min = min(depth), # minimum depth for mesosup layer
    depth_max = max(depth), # maximum depth for mesosup layer
    range = depth_max - depth_min, # depth range of mesosup layer
    range_max = 497.5 - depth_min, # maximum range of mesosup layer (from epipelagic layer to last bin before 500 m)
    .after = depth
    ) %>% 
  filter(range >= 0.8*range_max) %>% # keep only profiles which cover >80% of maximum range
  summarise_at(vars(Acantharea:kd490), mean, na.rm = TRUE) %>% # compute mean concentration for all plankton groups and env data
  ungroup()%>% 
  mutate(layer = "mesosup", .after = psampleid)
message(nrow(mesosup), " profiles for mesosup layer")


# Mesopelagic superior layer
mesoinf <- all %>% 
  filter(depth > 500 & depth < 1000) %>% # keep bins between 500 and 1000 m
  group_by(title, profile, lat, lon, datetime, psampleid) %>% # group by profile
  mutate( # for each profile compute:
    depth_min = min(depth), # minimum depth for mesoinf layer
    depth_max = max(depth), # maximum depth for mesoinf layer
    range = depth_max - depth_min, # depth range of mesoinf layer
    range_max = 997.5 - 502.5, # maximum range of mesoinf layer (from first bin after 500 m to last bin before 1000 m)
    .after = depth
  ) %>% 
  filter(range >= 0.8*range_max) %>% # keep only profiles which cover >80% of maximum range
  summarise_at(vars(Acantharea:kd490), mean, na.rm = TRUE) %>% # compute mean concentration for all plankton groups and env data
  ungroup()%>% 
  mutate(layer = "mesoinf", .after = psampleid)
message(nrow(mesoinf), " profiles for mesoinf layer")


# Bathypelagic layer
bathy <- all %>% 
  filter(depth > 1000) %>% # keep bins deeper than 1000 m
  group_by(title, profile, lat, lon, datetime, psampleid) %>% # group by profile
  summarise_at(vars(Acantharea:kd490), mean, na.rm = TRUE) %>% # compute mean concentration for all plankton groups and env data
  ungroup() %>% 
  mutate(layer = "bathy", .after = psampleid)
message(nrow(bathy), " profiles for bathy layer")


# Put everything together
all_in_layers <- bind_rows(epi, mesosup, mesoinf, bathy)


## Compute day/night and season ----
#--------------------------------------------------------------------------#
# Get sun altitude from lat, lon and datetime. 
# Set day_night to "day" if sun is above -6° below horizon (after nautical dawn and before nautical dusk)
# Compute season from latitude and month

# Read seasons definition
seasons <- read_csv("data/raw/seasons.csv") %>% 
  #select(-zone) %>% 
  gather(`1`:`12`, key = "month", value = "prod") %>% # prod is for productivity
  rename(lat_zone = zone) %>% 
  mutate(month = as.numeric(month))

# Latitude zones are:
# - "North frigid": North pole to 66.5°N (Arctic circle)
# - "North temperate": 66.5°N to 23.5°N
# - "Torrid": 23.5°N to 23.5°S
# - "South temperate": 23.5°S to 66.5°S (Antarctic circle)
# - "South frigid": 66.5°S to South pole

# Compute sun altitude
sun_alt <- getSunlightPosition(data = select(my_profiles, date=datetime, lat=lat, lon=lon), keep = "altitude") %>% 
  as_tibble() %>% 
  rename(datetime = date)

# Compute day/night and productivity
group_profiles <- my_profiles %>% 
  left_join(sun_alt) %>% 
  mutate(day_night = ifelse(altitude > -deg2rad(6), "day", "night")) %>% # compute day_night
  select(-altitude) %>% 
  mutate(
    month = month(datetime), # compute month
    lat_zone = ifelse( # compute latitude zone
      lat > 66.5,
      "North frigid", # if above 66.5°N, "North frigid"
      ifelse(
        lat > 23.5, # if above 23.5°N, "North temperate"
        "North temperate",
        ifelse(
          lat > -23.5, # if above 23.5°S, "Torrid"
          "Torrid",
          ifelse(
            lat > -66.5, # if above 66.5°S, "South temperate"
            "South temperate",
            "South frigid" # else (below 66.5°S), "South frigid"
          )
        )
      )
    )
  ) %>% 
  left_join(seasons) %>% 
  select(-c(month, lat_zone))

# Add to table with all data
all_in_layers <- all_in_layers %>% 
  left_join(group_profiles) %>% 
  relocate(day_night, prod, .after = psampleid)


## Join with bio-regions ----
#--------------------------------------------------------------------------#
all_data <- all_in_layers %>% 
  left_join(provinces) %>% 
  relocate(b_lat, long_code, long, gbp, gbp_code, .after = prod)


## Save data ----
#--------------------------------------------------------------------------#
save(all_data, file = "data/06.all_data.Rdata")
