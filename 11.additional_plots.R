#--------------------------------------------------------------------------#
# Project: planktyp
# Script purpose: Generate additional plots for the paper
# Date: 04/02/2021
# Author: Thelma Panaiotis
#--------------------------------------------------------------------------#

source("lib/set_up.R")
library(feather)
library(scales)

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
  mutate(layer = ifelse(layer == "epi", "Epipelagic", "Mesopelagic"))
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
  labs(x = "Plankton clusters", y = "Taxon", fill = "Proportion \n") +
  # Legend position at the bottom
  theme(legend.position = "bottom") +
  facet_wrap(~layer)
ggsave("plots/zoo/11.plankton_clusters_composition.png")  
ggsave("plots/paper/11.plankton_clusters_composition.svg")  

# Other version with difference beween actual and equal proportions
comp %>% 
  mutate(
    equal_prop = 1/length(unique(taxon)), # compute equal proportions of all taxa
    diff_prop = prop - equal_prop, # compute difference between prop and equal prop
    ) %>% 
  ggplot() +
  geom_tile(aes(x = clust_zoo, y = taxon, fill = diff_prop)) +
  scale_fill_gradient2(high = muted("blue"), low = muted("red")) +
  # Choose a gray color scale
  #scale_fill_distiller(palette = "Greys", direction = 1, limits = c(0,1)) +
  # Minimal theme
  theme_minimal() +
  # Remove background grid
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  # Reorder taxon
  scale_y_discrete(limits = rev(levels(comp$taxon))) +
  # Rename axes and legend
  labs(x = "Plankton clusters", y = "Taxon", fill = "Difference to \nequal proportions") +
  # Legend position at the bottom
  theme(legend.position = "bottom") +
  facet_wrap(~layer)



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

# Keep only data from these taxa
cce_epi <- cce_epi %>% 
  gather(Acantharea:Trichodesmium, key = "taxon", value = "conc") %>% 
  filter(taxon %in% cce_taxa$taxon) 

# Test if there is a difference in concentrations
# Wilcoxon–Mann–Whitney test because concentrations are far from a normal distribution
p_val <- c() # initiate empty vector to store p-values
for (taxa in cce_taxa$taxon){ # Loop over taxa
  cce_day <- cce_epi %>% filter(taxon == taxa) %>% filter(day_night == "day") # extract day concentrations
  cce_night <- cce_epi %>% filter(taxon == taxa) %>% filter(day_night == "night") # extract night concentrations
  
  test_res <- wilcox.test(x = cce_day$conc, y = cce_night$conc) # perform test
  
  p_val <- c(p_val, test_res$p.value) # store p_value
}

# Store test p-values and taxa in a table
cce_p_vals <- tibble(
  taxon = cce_taxa$taxon,
  p_val = p_val
) %>% 
  mutate( # compute significance for various thesholds
    sig_05 = p_val < 0.05, # *
    sig_01 = p_val < 0.01, # **
    sig_001 = p_val < 0.001, # ***
    )
cce_p_vals
# Concentrations are different between day and night for Copepoda and Eumalacostraca, and maybe for Phaeodaria. 


# Boxplot of concentrations
cce_epi %>% 
  ggplot() +
  geom_boxplot(aes(x = taxon, y = conc, color = day_night)) +
  # y scale as log
  scale_y_continuous(trans = "log1p") +
  theme_classic() +
  labs(x = "Taxon", y = expression(paste("Concentration (ind.", L^{-1}, ")")), color = "")
ggsave("plots/zoo/11.cce_conc_day_night.png")  
ggsave("plots/paper/11.cce_conc_day_night.svg")


# Check high value of copepoda at day
cce_psample <- cce_epi %>% 
  filter(day_night == "day" & taxon == "Copepoda") %>% 
  arrange(desc(conc)) %>% 
  slice(1) %>% 
  pull(psampleid)

cce_prof <- all_data %>% filter(psampleid == cce_psample) %>% select(title:datetime) %>% pull(profile)
o %>% 
  select(title:profile, lat, lon, datetime, depth, taxon, lineage, group, group_lineage) %>% 
  filter(profile == cce_prof & group == "Copepoda") %>% 
  filter(depth < 200) %>% 
  arrange(depth)
# 130 copepods between surface and 200 m  

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
