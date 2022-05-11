#--------------------------------------------------------------------------#
# Project: planktyp
# Script purpose: Clean and reformat bio data: split between zoo and detritus; extract parts data; compute concentrations
# Date: 29/01/2021
# Author: Thelma Panaïotis
#--------------------------------------------------------------------------#

source("lib/set_up.R")
library(feather)
library(ecotaxar)


## Read extracted objects and add project and sample info ----
#--------------------------------------------------------------------------#
load("data/00.extraction.Rdata")


## Keep only relevant taxa for zoo and det tables ----
#--------------------------------------------------------------------------#
# zoo: only living objects in relevant groups, biological data to analyze
zoo <- o %>%
  filter(!taxon %in% c("?", "artefact", "bubble", "misc?", "detritus", NA)) %>% 
  select(title, profile, sampleid, lat, lon, depth, taxon) %>% 
  group_by(title, profile, sampleid, lat, lon, depth, taxon) %>% 
  summarise(abund = n()) %>% 
  ungroup() %>% 
  arrange(title, profile, depth, taxon)

# det: only detritus, used as environmental descriptor of richness in marine snow
det <- o %>%
  filter(taxon == "detritus") %>% 
  select(title, profile, sampleid, lat, lon, depth, taxon) %>% 
  group_by(title, profile, sampleid, lat, lon, depth, taxon) %>% 
  summarise(abund = n()) %>% 
  ungroup() %>% 
  arrange(title, profile, depth, taxon)


## Delete Trichodesmium below 500 m ----
#--------------------------------------------------------------------------#
# After inspections, these objects are nor puffs neither tuffs
deep_tricho <- zoo %>% filter(taxon == "Trichodesmium" & depth > 500) 
deep_tricho
zoo <- zoo %>% filter(!(taxon == "Trichodesmium" & depth > 500))
message("Removed ", sum(deep_tricho$abund), " bins with deep Trichodesmiums (below 500 m)")


## Read particles concentration in sampled bins ----
#--------------------------------------------------------------------------#
part_conc <- read_csv("data/raw/part_conc.csv")


## Get all sampled bins and compute concentrations for zoo and detritus ----
#--------------------------------------------------------------------------#
# Nice thing is that particles table contains all sampled bins with their watervolume
bins <- part_conc %>% select(title:watervolume)

# For zoo data, generate all combinations of bins by taxa
zoo_taxon <- zoo %>% pull(taxon) %>% unique()
zoo_conc <- crossing(bins, taxon = zoo_taxon) %>% 
  left_join(zoo, by = c("title", "profile", "sampleid", "lat", "lon", "depth", "taxon")) %>% 
  replace_na(list(abund = 0)) %>% 
  mutate(conc = abund / watervolume) %>% 
  select(-c(watervolume, abund)) %>% 
  spread(key = taxon, value = conc)


# For det and parts, proceed to a join on bins data
det_conc <- bins %>% 
  left_join(det, by = c("title", "profile", "sampleid", "lat", "lon", "depth")) %>% 
  replace_na(list(abund = 0, taxon = "detritus")) %>% 
  mutate(conc = abund / watervolume) %>% 
  select(-c(watervolume, abund)) %>% 
  spread(key = taxon, value = conc)


## Save data ----
#--------------------------------------------------------------------------#
# Save concentrations
save(zoo, det, zoo_conc, det_conc, part_conc, file="data/01.bio_data.Rdata")

