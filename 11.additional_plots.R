#--------------------------------------------------------------------------#
# Project: planktyp
# Script purpose: Generate additional plots for the paper
# Date: 04/02/2021
# Author: Thelma Panaiotis
#--------------------------------------------------------------------------#

# TODO check tricho below 500 m

source("lib/set_up.R")
library(feather)

load("data/01.bio_data.Rdata")
load("data/06.all_data.Rdata")


## Get plankton organisms in selected profiles and groups ----
#--------------------------------------------------------------------------#
# List samples and psamples used in analysis
my_samples <- all_data %>% 
  select(psampleid) %>% 
  unique() %>% 
  left_join(zoo_conc %>% select(psampleid, sampleid) %>% unique())

# List taxa used in analysis
my_taxa <- zoo_conc %>% select(Acantharea:Trichodesmium) %>% colnames()
# add detritus
my_taxa <- c(my_taxa, "detritus")

# Read table with all objects
o <- read_feather("data/00.all_zoo.feather")

# Keep only objects in selected profiles
o <- o %>% filter(sampleid %in% my_samples$sampleid)

# Keep only objects in selected taxa
o <- o %>% filter(group %in% my_taxa)


## Bar chart of plankton organisms count ----
#--------------------------------------------------------------------------#
# Plot bar chart
o %>% 
  count(group) %>% 
  ungroup() %>% 
  filter(group != "detritus") %>% # ignore detritus 
  arrange(desc(n)) %>% # sort by descending number of objects
  mutate(group = factor(group, group)) %>% # group as factor for a nice plot
  ggplot() +
  geom_col(aes(x = group, y = n)) +
  theme_classic() +
  # rotate x axis text
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # remove white space on y axis
  scale_y_continuous(expand = c(0,0)) +
  # Axes labels
  xlab("Taxonomic group") +
  ylab("n organisms") +
  ggtitle("Plankton dataset composition")
ggsave("plots/zoo/11.plankton_dataset_comp.png")  

# Plot bar chart with log y axis
o %>% 
  count(group) %>% 
  ungroup() %>% 
  filter(group != "detritus") %>% # ignore detritus 
  arrange(desc(n)) %>% # sort by descending number of objects
  mutate(group = factor(group, group)) %>% # group as factor for a nice plot
  ggplot() +
  geom_col(aes(x = group, y = n)) +
  theme_classic() +
  # rotate x axis text
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # remove white space on y axis
  scale_y_continuous(trans = "log1p", expand = c(0,0)) +
  # Axes labels
  xlab("Taxonomic group") +
  ylab("n organisms") +
  ggtitle("Plankton dataset composition") +
  ggsave("plots/zoo/11.plankton_dataset_comp_log.png") 


## Maps of plankton distribution ----
#--------------------------------------------------------------------------#


## Look for deep Trichodesmiums ----
#--------------------------------------------------------------------------#
zoo_conc %>% 
  select(c(title:depth, Trichodesmium))

tricho <- o %>% 
  select(c(title:datetime, depth, taxon = group)) %>% 
  filter(taxon == "Trichodesmium")
summary(tricho)
# Most Trichodesmium are < 50 m but a few are very deep
tricho %>% filter(depth > 500) %>% select(title, profile, sampleid) %>% unique()
# 169 profiles have Tricho > 500 m


tricho %>%   
  ggplot() +
  geom_histogram(aes(x = depth), binwidth = 10) +
  scale_y_continuous(trans = "log1p") +
  ggtitle("Histogram of Trichodesmium depth occurences")
ggsave("plots/zoo/11.histo_depth_tricho.png") 


## Map of projects ----
#--------------------------------------------------------------------------#


## Concentrations at day and night for CCE ----
#--------------------------------------------------------------------------#
