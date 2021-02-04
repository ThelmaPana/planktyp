#--------------------------------------------------------------------------#
# Project: plantyp
# Script purpose: Read and process satellite data
# Date: 02/02/2021
# Author: Thelma Panaiotis
#--------------------------------------------------------------------------#

source("lib/set_up.R")

load("data/03.layers.Rdata")


## Read satellite data ----
#--------------------------------------------------------------------------#
# Read, rename, replace missing data by NA and ignore them
sat3 <- read_csv("data/raw/satellite/matchup_sat_uvp_config_3.csv") %>% 
  rename(chla_surf = Chla, bbp = Bbp, poc = POC, pic = PIC, par = daily_PAR, kd490 = Kd490) %>% 
  replace(., is.na(.), NA) %>% 
  filter(!is.na(chla_surf))


## Match with profile selection ----
#--------------------------------------------------------------------------#
# This is made difficult because profile names have changed since satellite data extraction.

profiles <- layers %>% select(title, profile, lat_ctd = lat, lon_ctd = lon, datetime, psampleid)
message(nrow(profiles), " profiles of CTD data")

sat3 <- sat3 %>% select(-c(datetime)) %>% rename(lat_sat = lat, lon_sat = lon)
message(nrow(sat3), " profiles with satellite data")

joined <- inner_join(profiles, sat3)
message(nrow(joined), " profiles joined between CTD and satellite data with profile names")

# Look for profiles in one dataset and not in the other
# in ctd but not in sat
ctd_not_sat <- anti_join(profiles, sat3)
ctd_not_sat %>% pull(title) %>% unique()
unique(ctd_not_sat$title)
# CTD data could not be matched with satellite data for 5 projects.
# Let's inspect individually to see if more profiles can be matched. 

# in sat but not in ctd
sat_not_ctd <- anti_join(sat3, profiles)
unique(sat3$project)

ggplot() +
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "gray") +
  geom_point(data = ctd_not_sat,  aes(x = lon_ctd, y = lat_ctd, color = title)) +
  coord_quickmap() +
  theme(legend.position = "bottom",  legend.direction = "vertical") +
  ggtitle("Projects with CTD data not joined with satellite data")


## Specific matches ----
#--------------------------------------------------------------------------#

## Match for "UVP5 Geomar 2014 m105"  
ctd_m105 <- ctd_not_sat %>% filter(title == "UVP5 Geomar 2014 m105")
sat_m105 <- sat_not_ctd %>% filter(project == "M 105")

joined_m105 <- ctd_m105 %>% 
  left_join(sat_m105 %>% mutate(profile = str_split_fixed(profile, "_", n=2)[,2])) %>% 
  mutate(
    lat_diff = abs(lat_ctd - lat_sat),
    lon_diff = abs(lon_ctd - lon_sat),
    ) %>% 
  select(title, profile, lat = lat_ctd, lon = lon_ctd, psampleid, datetime, chla_surf:lon_diff) %>% 
  filter(!is.na(lat_diff))


## Match for "UVP5 Geomar 2015 m116"  
# Nothing to join for m116: some profiles are missing from satellite data
ctd_m116 <- ctd_not_sat %>% filter(title == "UVP5 Geomar 2015 m116")
# no profiles of M116 in left satellite data

## Match for "UVP5 Sargasso 2014 (sargasso_a sargasso_b)"
ctd_sarg <- ctd_not_sat %>% filter(title == "UVP5 Sargasso 2014 (sargasso_a sargasso_b)")
sat_sarg <- sat_not_ctd %>% filter(project == "SARGASSO A")

joined_sarg <- ctd_sarg %>% 
  left_join(sat_sarg %>% mutate(profile = str_split_fixed(profile, "_", n=2)[,2])) %>% 
  mutate(
    lat_diff = abs(lat_ctd - lat_sat),
    lon_diff = abs(lon_ctd - lon_sat),
  ) %>% 
  select(title, profile, lat = lat_ctd, lon = lon_ctd, psampleid, datetime, chla_surf:lon_diff) %>% 
  filter(!is.na(lat_diff))


## Match for "UVP5 Tara Oceans (tara2009, tara2010, tara2011, tara2012, tara2013)"
# Nothing to join for Tara
ctd_tara <- ctd_not_sat %>% filter(title == "UVP5 Tara Oceans (tara2009, tara2010, tara2011, tara2012, tara2013)")
sat_tara <- sat_not_ctd %>% filter(project %in% c("TARA 2009", "TARA 2010", "TARA 2012 A", "TARA 2012 B", "Tara 2011", "TARA 2013"))
# no profile in common between left ctd and satellite data


## Match for "UVP5hd GreenEdge 2016"
# Nothing to join for Greenedge 2016
ctd_ge <- ctd_not_sat %>% filter(title == "UVP5hd GreenEdge 2016")
sat_ge <- sat_not_ctd %>% filter(project %in% c("GREENEDGE 2016 A", "GREENEDGE 2016 C"))
# no profile in common between left ctd and satellite data


## Match for "UVP5zd KEOPS"
ctd_keops <- ctd_not_sat %>% filter(title == "UVP5zd KEOPS")
sat_keops <- sat_not_ctd %>% filter(project == "KEOPS 2")

joined_keops <- ctd_keops %>% 
  left_join(sat_keops %>% mutate(profile = str_split_fixed(profile, "_", n=2)[,2])) %>% 
  mutate(
    lat_diff = abs(lat_ctd - lat_sat),
    lon_diff = abs(lon_ctd - lon_sat),
  ) %>% 
  select(title, profile, lat = lat_ctd, lon = lon_ctd, psampleid, datetime, chla_surf:lon_diff) %>% 
  filter(!is.na(lat_diff))
# no profiles of naames in left satellite data


## Match for "UVP5hd NAAMES02"  
# Nothing to join for naames
ctd_naames <- ctd_not_sat %>% filter(title == "UVP5hd NAAMES02")
sat_naames <- sat_not_ctd %>% filter(project == "NAAMES 02")


## Put back all matches together ----
#--------------------------------------------------------------------------#
sat <- bind_rows(
  joined %>% select(title, profile, lat = lat_ctd, lon = lon_ctd, datetime, psampleid, chla_surf, bbp, poc, pic, par, kd490),
  joined_m105 %>% select(title, profile, lat, lon, datetime, psampleid, chla_surf, bbp, poc, pic, par, kd490),
  joined_sarg %>% select(title, profile, lat, lon, datetime, psampleid, chla_surf, bbp, poc, pic, par, kd490),
  joined_keops %>% select(title, profile, lat, lon, datetime, psampleid, chla_surf, bbp, poc, pic, par, kd490),
  )

message(nrow(sat), " profiles joined between CTD and satellite data after specific matches")
# Some data are missing for pic or par but keep them still


## Save data ----
#--------------------------------------------------------------------------#
save(sat, file = "data/04.sat.Rdata")

