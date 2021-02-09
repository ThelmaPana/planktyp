#--------------------------------------------------------------------------#
# Project: planktyp
# Script purpose: Clean and reformat bio data: split between zoo and detritus; extract parts data; compute concentrations
# Date: 29/01/2021
# Author: Thelma Panaïotis
#--------------------------------------------------------------------------#

# TODO: check beginning depth of profiles (Tara)

source("lib/set_up.R")
library(feather)
library(ecotaxar)


## Read extracted objects and add project and sample info ----
#--------------------------------------------------------------------------#
o <- read_feather("data/00.all_zoo.feather")

# Keep only useful columns
o <- o %>% select(title, profile, sampleid, lat, lon, depth=mid_depth_bin, taxon, lineage, group, group_lineage)
# List taxonomy
o_taxon <- o %>% select(taxon, group) %>% unique() %>% arrange(taxon)

# Plot a map of all available UVP profiles
uvp_profiles <- o %>% 
  select(title, sampleid, profile, lat, lon) %>% 
  unique() 

uvp_profiles %>%   
  ggplot() +
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "gray") +
  geom_point(aes(x = lon, y = lat)) +
  coord_quickmap() +
  ggtitle("All UVP profiles") +
  ggsave("plots/zoo/00.all_uvp_profiles.png")


## Keep only relevant taxa for zoo and det tables ----
#--------------------------------------------------------------------------#
# zoo: only living objects in relevant groups, biological data to analyze
zoo <- o %>%
  filter(!group %in% c("?", "artefact", "bubble", "misc?", "detritus", NA)) %>% 
  select(title, profile, lat, lon, depth, taxon=group) %>% 
  group_by(title, profile, lat, lon, depth, taxon) %>% 
  summarise(abund = n()) %>% 
  ungroup() %>% 
  arrange(title, profile, depth, taxon)

# det: only detritus, used as environmental descriptor of richness in marine snow
det <- o %>%
  filter(group == "detritus") %>% 
  select(title, profile, lat, lon, depth, taxon=group) %>% 
  group_by(title, profile, lat, lon, depth, taxon) %>% 
  summarise(abund = n()) %>% 
  ungroup() %>% 
  arrange(title, profile, depth, taxon)


## Delete Trichodesmium below 500 m ----
#--------------------------------------------------------------------------#
# After inspections, these objects are nor puffs neither tuffs
deep_tricho <- zoo %>% filter(taxon == "Trichodesmium" & depth > 500) 
zoo <- zoo %>% filter(!(taxon == "Trichodesmium" & depth > 500))
message("Removed ", sum(deep_tricho$abund), " deep Trichodesmiums (below 500 m)")


## Extract particles data along with sampled bins ----
#--------------------------------------------------------------------------#
# Objects are extracted from ecotaxa but watervolumes are associated with ecopart data. 
# Use sampleid and psampleid to establish link between these. 
db <- db_connect_ecotaxa()

# Link ecotaxa sampleid to ecopart psampleid
psamples <- tbl(db, "part_samples") %>% select(psampleid, sampleid) %>% filter(sampleid %in% !!unique(o$sampleid)) %>% collect()
# 1 sampleid = 1 psampleid, no duplicates

# Add psampleid to objects table to make life easier
o <- left_join(o, psamples) %>% relocate(psampleid, .after = sampleid)

# Particles data are stored in part_histopart_XXX. Use small table as we’ll sum all classes. 
parts <- tbl(db, "part_histopart_reduit") %>% 
  select(-c(lineno, datetime)) %>% 
  filter(psampleid %in% !!psamples$psampleid) %>% # keep only relevant psamples
  collect() %>% 
  left_join(o %>% select(title, profile, sampleid, psampleid, lat, lon) %>% unique()) %>%  # add title, profile name and coordinates
  select(title, profile, lat, lon, sampleid, everything()) # reorder columns

# Sum abundances and biovolumes over size classes and compute concentrations
part_conc <- parts %>% 
  mutate(
    sum_abund = reduce(select(., c(class01:class15)), `+`), # sum abundances for all size classes
    conc = sum_abund / watervolume,
    sum_biovol = reduce(select(., c(biovol01:biovol15)), `+`),
    biovol = sum_biovol / watervolume,
    ) %>% 
  select(-c(class01:class15, biovol01:biovol15, sum_abund, sum_biovol))


## Get all sampled bins and compute concentrations for zoo and det ----
#--------------------------------------------------------------------------#
# Nice thing is that particles table contains all sampled bins with their watervolume
bins <- parts %>% select(title:watervolume)

# For zoo data, generate all combinations of bins by taxa
zoo_taxon <- zoo %>% pull(taxon) %>% unique()
zoo_conc <- crossing(bins, taxon = zoo_taxon) %>% 
  left_join(zoo, by = c("title", "profile", "lat", "lon", "depth", "taxon")) %>% 
  replace_na(list(abund = 0)) %>% 
  mutate(conc = abund / watervolume) %>% 
  select(-c(watervolume, abund)) %>% 
  spread(key = taxon, value = conc)


# For det and parts, proceed to a join on bins data
det_conc <- bins %>% 
  left_join(det, by = c("title", "profile", "lat", "lon", "depth")) %>% 
  replace_na(list(abund = 0, taxon = "detritus")) %>% 
  mutate(conc = abund / watervolume) %>% 
  select(-c(watervolume, abund)) %>% 
  spread(key = taxon, value = conc)


## Save data ----
#--------------------------------------------------------------------------#
# Save concentrations
save(zoo_conc, det_conc, part_conc, file="data/01.bio_data.Rdata")

# and disconnect from database
db_disconnect_ecotaxa(db)
