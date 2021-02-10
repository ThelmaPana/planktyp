#--------------------------------------------------------------------------#
# Project: planktyp
# Script purpose: Generate additional plots for the paper
# Date: 04/02/2021
# Author: Thelma Panaiotis
#--------------------------------------------------------------------------#

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


## Inspect dataset composition ----
#--------------------------------------------------------------------------#
# Dates
all_data %>% select(datetime) %>% summary()

# Number of profiles
length(unique(all_data$psampleid))
all_data %>% select(psampleid, layer) %>% count(layer)

# Number of objects
nrow(o) # objects
nrow(o %>% filter(group != "detritus")) # living objects

# Number and list of taxonomic groups
length(unique(o$group)) # 29 groups including detritus --> 28 taxonomic groups
sort(unique(o$group))


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
ggsave("plots/paper/11.plankton_dataset.svg")

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


## Plankton clusters composition ----
#--------------------------------------------------------------------------#
load("data/07.epi_plankton_cluster_comp.Rdata")
epi_comp <- comp
load("data/08.mesosup_plankton_cluster_comp.Rdata")
mesosup_comp <- comp

comp <- bind_rows(epi_comp, mesosup_comp) %>% 
  mutate(layer = ifelse(layer == "epi", "Epipelagic", "Mesopelagic sup"))
comp %>% 
  ggplot() +
  geom_tile(aes(x = clust_zoo, y = taxon, fill = prop)) +
  # Choose a gray color scale
  scale_fill_distiller(palette = "Greys", direction = 1, limits = c(0,1)) +
  # Minimal theme
  theme_minimal() +
  # Remove background grid
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  # Reorder taxon
  scale_y_discrete(limits = rev(levels(comp$taxon))) +
  # Rename axes and legend
  labs(x = "Plankton clusters", y = "Taxon", fill = "Proportion") +
  facet_wrap(~layer)
ggsave("plots/zoo/11.plankton_clusters_composition.png")  
ggsave("plots/paper/11.plankton_clusters_composition.svg")  


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


## Check watervolumes ----
#--------------------------------------------------------------------------#
summary(part_conc$watervolume)
# max watervolume is > 25000. This seems wrong
part_conc %>% filter(watervolume > 500) # 212 bins with watervolume > 500
part_conc %>% filter(watervolume > 1000) # 93 bins with watervolume > 500

part_conc %>% 
  ggplot() +
  geom_histogram(aes(x = watervolume), binwidth = 100) +
  scale_y_continuous(trans = "log1p")


## Map of stations ----
#--------------------------------------------------------------------------#
all_data %>% 
  select(title, profile, lat, lon,  psampleid) %>% 
  unique() %>% 
  ggplot() +
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "gray") +
  geom_point(aes(x = lon, y = lat), size = 0.5) +
  coord_quickmap() +
  theme_minimal() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x= "Longitude", y = "Latitude")
ggsave("plots/zoo/11.map_profiles.png")  
ggsave("plots/paper/11.map_profiles.svg")


# Try to color according to projects
all_data %>% 
  select(title, profile, lat, lon,  psampleid) %>% 
  unique() %>% 
  ggplot() +
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "gray") +
  geom_point(aes(x = lon, y = lat, color = title), size = 0.5) +
  coord_quickmap() +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x= "Longitude", y = "Latitude") 


## Concentrations at day and night for CCE ----
#--------------------------------------------------------------------------#
# Get concentrations at CCE in epipelagic layer
cce_epi <- all_data %>% 
  filter(layer == "epi") %>% # concentrations in epipelagic layer
  filter(str_detect(title, "CCELTER")) %>% # at CCELTER
  select(psampleid, datetime, day_night, Acantharea:Trichodesmium)

# Look for the 5 most abundant taxa
cce_taxa <- cce_epi %>% 
  gather(Acantharea:Trichodesmium, key = "taxon", value = "conc") %>% 
  group_by(taxon) %>% 
  summarise(conc = mean(conc)) %>% 
  ungroup() %>% 
  arrange(desc(conc)) %>% 
  slice(1:5) 
cce_taxa

cce_epi %>% 
  gather(Acantharea:Trichodesmium, key = "taxon", value = "conc") %>% 
  filter(taxon %in% cce_taxa$taxon) %>% 
  ggplot() +
  geom_boxplot(aes(x = taxon, y = conc, color = day_night)) +
  # y scale as log
  scale_y_continuous(trans = "log1p") +
  theme_classic() +
  labs(x = "Taxon", y = expression(paste("Concentration (ind.", L^{-1}, ")")), color = "")
ggsave("plots/zoo/11.cce_conc_day_night.png")  
ggsave("plots/paper/11.cce_conc_day_night.svg")


## Look at anosim results ----
#--------------------------------------------------------------------------#
# Load results
load("data/07.epi_anosim.Rdata")
anosim_epi <- anosim_res
load("data/08.mesosup_anosim.Rdata")
anosim_meso <- anosim_res

anosim_epi
anosim_meso

# Number of selected profiles for anosim
load("data/07.epi_partitionings.Rdata")
groups_epi <- groups
nrow(groups_epi) # number of included profiles
groups_epi %>% count(day_night) # number of day and night profiles
groups_epi %>% count(prod) # number of productive and non-productive profiles

load("data/08.mesosup_partitionings.Rdata")
groups_meso <- groups 
nrow(groups_meso) # number of included profiles
groups_meso %>% count(day_night) # number of day and night profiles
