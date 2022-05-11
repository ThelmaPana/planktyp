#--------------------------------------------------------------------------#
# Project: planktyp
# Script purpose: Prepare data from Laetitia’s extraction
# Date: 01/02/2021
# Author: Thelma Panaiotis
#--------------------------------------------------------------------------#

source("lib/set_up.R")
library(tidyverse)
library(googlesheets4)
library(data.tree)
gs4_auth(use_oob = TRUE)
1

# Grouping sheet
ss <- "https://docs.google.com/spreadsheets/d/1o7D5-Q-UBZw3JgaNpsUpH3-5ineLdWR1kBreYu5cSGs/edit?usp=sharing"
update_sheet <- FALSE


## Laetita data ----
#--------------------------------------------------------------------------#
# Load objects
load("data/raw/uvp_objects.Rdata")

# Select relevant columns
obj <- obj %>% 
  select(title, projid, sampleid, profile, objid, depth, datetime, lon, lat, taxon, lineage) %>% 
  # mutate depth bin
  mutate(depth = floor(depth/5)*5 + 2.5)


## Generate grouping sheet ----
#--------------------------------------------------------------------------#
if (update_sheet){
  # Build the tree 
  tc <- count(obj, taxon, lineage) %>% 
    # convert it into a tree
    rename(pathString=lineage) %>%
    arrange(pathString) %>%
    as.Node()
  
  print(tc, "taxon","n", limit = 50)
  # Convert to dataframe
  tcd <- ToDataFrameTree(tc, "taxon", "n")%>% 
    as_tibble() %>% 
    rename(level0=taxon, nb_level0=n)
  
  # Write new tree
  range_write(ss, data = tcd, sheet = "counts") 
  # Open it in browser tab to make sure everything is ok
  gs4_browse(ss)
}


## Regroup organisms ----
#--------------------------------------------------------------------------#
# Read google sheet with groups
groups <- read_sheet(ss, sheet = "counts") %>% 
  select(taxon, group) %>% 
  drop_na(taxon)

# Join objects with taxo groups
obj <- obj %>% 
  left_join(groups, by = "taxon") %>% 
  drop_na(group) %>% # ignore objects not matched with a group
  select(-taxon) %>% # use group as new taxon
  rename(taxon = group)


## Select only my samples ----
#--------------------------------------------------------------------------#
# Read list of my samples
samples <- read_tsv("data/raw/UVP5_samples_selected.tsv")
# Keep only objects in my samples
o <- obj %>% filter(sampleid %in% samples$sampleid)


## Save extracted data ----
#--------------------------------------------------------------------------#
save(o, file = "data/00.extraction.Rdata")
