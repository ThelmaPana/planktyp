#--------------------------------------------------------------------------#
# Project: planktyp
# Script purpose: Generate additional plots for the paper
# Date: 04/02/2021
# Author: Thelma Panaiotis
#--------------------------------------------------------------------------#

source("lib/set_up.R")
library(feather)
library(scales)
library(vegan)
#library(ggsankey)
library(cmocean)
library(castr)
library(fuzzyjoin)

load("data/01.bio_data.Rdata")
load("data/06.all_data.Rdata")
load("data/03.layers.Rdata")
#load("data/00.extraction.Rdata")


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
#o <- read_feather("data/00.all_zoo.feather")

# Keep only objects in selected profiles
zoo <- zoo %>% filter(sampleid %in% my_samples$sampleid)

# Keep only objects in selected taxa
zoo <- zoo %>% filter(taxon %in% my_taxa)


## Inspect dataset composition ----
#--------------------------------------------------------------------------#
# Dates
all_data %>% select(datetime) %>% summary()

# Number of profiles
length(unique(all_data$psampleid))
all_data %>% select(psampleid, layer) %>% count(layer)

# Number of objects
sum(zoo$abund) # living objects
sum(det$abund) # detritus
sum(det$abund) + sum(zoo$abund) # all objects

# Number and list of taxonomic groups
length(unique(zoo$taxon)) # 28 groups including detritus --> 27 taxonomic groups
sort(unique(zoo$taxon))


## Bar chart of plankton organisms count ----
#--------------------------------------------------------------------------#
# Plot bar chart
zoo %>% 
  group_by(taxon) %>% 
  summarise(n = sum(abund)) %>% 
  ungroup() %>% 
  filter(taxon != "detritus") %>% # ignore detritus 
  arrange(desc(n)) %>% # sort by descending number of objects
  mutate(taxon = factor(taxon, taxon)) %>% # taxon as factor for a nice plot
  ggplot() +
  geom_col(aes(x = taxon, y = n)) +
  theme_classic() +
  # rotate x axis text
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # remove white space on y axis
  scale_y_continuous(expand = c(0,0)) +
  # Axes labels
  xlab("Taxonomic group") +
  ylab("Number of organisms") +
  ggtitle("Plankton dataset composition")
ggsave("plots/zoo/11.plankton_dataset_comp.png")  
ggsave(file = "plots/paper/11.plankton_dataset.pdf", width = 166, height = 100, unit = "mm", dpi = 300)

# Plot bar chart with log y axis
zoo %>% 
  group_by(taxon) %>% 
  summarise(n = sum(abund)) %>% 
  ungroup() %>% 
  filter(taxon != "detritus") %>% # ignore detritus 
  arrange(desc(n)) %>% # sort by descending number of objects
  mutate(taxon = factor(taxon, taxon)) %>% # taxon as factor for a nice plot
  ggplot() +
  geom_col(aes(x = taxon, y = n)) +
  theme_classic() +
  # rotate x axis text
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # remove white space on y axis
  scale_y_continuous(trans = "log1p", expand = c(0,0)) +
  # Axes labels
  xlab("Taxonomic group") +
  ylab("n organisms") +
  ggtitle("Plankton dataset composition")
#ggsave("plots/zoo/11.plankton_dataset_comp_log.png") 


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

# Other version with difference between actual and equal proportions
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

comp %>% 
  mutate(clust_zoo = ifelse(layer == "Epipelagic", paste0(clust_zoo, "e"), paste0(clust_zoo, "m"))) %>% 
  ggplot() +
  geom_col(aes(x = prop, y = taxon, fill = clust_zoo), position = "dodge") +
  scale_y_discrete(limits = rev(levels(comp$taxon))) +
  scale_x_continuous(limits = c(0,1), expand = c(0,0)) +
  # Rename axes and legend
  labs(x = "Proportion", y = "Taxon", fill = "Epipelagic \nPlankton \ncluster") +
  theme_minimal() +
  theme(panel.spacing = unit(2, "lines")) +
  facet_wrap(~layer)


## Maps of plankton distribution ----
#--------------------------------------------------------------------------#


## Look for deep Trichodesmiums ----
#--------------------------------------------------------------------------#
zoo_conc %>% 
  select(c(title:depth, Trichodesmium))

tricho <- zoo %>% 
  #select(c(title:datetime, depth, taxon)) %>% 
  filter(taxon == "Trichodesmium")
summary(tricho)

# Most Trichodesmium are < 50 m but a few are very deep
tricho %>% filter(depth > 500) %>% select(title, profile, sampleid) %>% unique()

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
  geom_point(aes(x = lon, y = lat), size = 1, shape = 1, alpha = 0.5) +
  coord_quickmap() +
  theme_minimal() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x= "Longitude", y = "Latitude", title = "Stations map")
ggsave("plots/zoo/11.map_profiles.png")  
ggsave(file = "plots/paper/11.map_profiles.pdf", width = 166, height = 123, unit = "mm", dpi = 300)

# Try to color according to projects
all_data %>% 
  mutate(year = factor(year(datetime))) %>% 
  select(title, profile, lat, lon,  psampleid, year) %>% 
  unique() %>% 
  ggplot() +
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "gray") +
  geom_point(aes(x = lon, y = lat, color = year), size = 0.5) +
  coord_quickmap() +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x= "Longitude", y = "Latitude") 

t_missions <- all_data %>% 
  select(title:psampleid) %>% 
  mutate(year = year(datetime)) %>% 
  unique() %>% 
  count(title, year)


## Hovmöller diagram of sampling ----
#--------------------------------------------------------------------------#
stations <- all_data %>% 
  select(title, profile, lat, lon,  psampleid, datetime) %>% 
  unique()


break.vec <- seq(from = as.Date("1970-01-01"), to = as.Date("1970-12-01"), by = "month")

stations %>% 
  mutate(yday = as.Date(yday(datetime), "1970-01-01")) %>% 
  ggplot(aes(x = yday, y = lat)) +
  geom_density_2d_filled(aes(x = yday, y = lat), bins = 12) +
  geom_point(aes(x = yday, y = lat), size = 0.5, alpha = 1) +
  #scale_x_date(breaks = break.vec, , date_labels="%b", expand = c(0,0)) +
  scale_x_date(breaks = break.vec, limits = c(min(break.vec), as.Date("1970-12-31")), date_labels="%b", expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(-80, 80)) +
  scale_fill_viridis_d(option = "G", direction = -1) +
  #scale_fill_brewer(palette = "YlGnBu") +
  labs(x = "Month", y = "Latitude", fill = "Density")
ggsave("plots/zoo/11.hovmoller.png")  
ggsave(file = "plots/paper/11.hovmoller.pdf", width = 166, height = 80, unit = "mm", dpi = 300)





## Generate pairs of day and night stations ----
#--------------------------------------------------------------------------#
# Separate day and night profiles
all_data_d <- all_data %>% filter(layer == "epi" & day_night == "day") %>% rename(psampleid_d = psampleid)
all_data_n <- all_data %>% filter(layer == "epi" & day_night == "night") %>% rename(psampleid_n = psampleid)


# Geographic distance
dist_geo <- geo_left_join(
  all_data_n %>% select(title, psampleid_n, lat, lon), 
  all_data_d %>% select(psampleid_d, lat, lon), 
  by = c("lat", "lon"), max_dist = 2, unit = "km", distance_col = "dist_geo"
  ) %>% 
  drop_na("dist_geo") %>% 
  mutate(pair = paste(psampleid_n, psampleid_d)) %>% 
  select(title, pair, dist_geo)
dist_geo

# Time distance
dist_dt <- difference_left_join(
  all_data_n %>% select(title, psampleid_n, datetime), 
  all_data_d %>% select(psampleid_d, datetime),
  by = "datetime", max_dist = 1440, distance_col = "dist_dt" #720 mins for 12h
) %>% 
  drop_na("dist_dt") %>% 
  mutate(pair = paste(psampleid_n, psampleid_d)) %>% 
  select(title, pair, dist_dt)

to_comp <- dist_geo %>% 
  inner_join(dist_dt) %>% 
  separate(pair, into = c("psampleid_n", "psampleid_d")) %>% 
  select(-dist_dt)

# Some night stations are joined with several day stations, select the closest day station
to_comp <- to_comp %>% 
  group_by(title, psampleid_n) %>% 
  filter(dist_geo == min(dist_geo)) %>% 
  ungroup() %>% 
  mutate(pair = row_number())
to_comp # 172 pairs of stations


## Concentrations at day and night for paired stations ----
#--------------------------------------------------------------------------#
# Most five abundant taxa
all_data %>% 
  filter(layer == "epi") %>% 
  pivot_longer(Acantharea:Trichodesmium, names_to = "taxon", values_to = "conc") %>% 
  group_by(taxon) %>% 
  summarise(conc = mean(conc)) %>% 
  ungroup() %>% 
  arrange(desc(conc))

taxa_comp <- c("Trichodesmium", "Copepoda", "Phaeodaria", "Acantharea")


# Night data
to_comp_n <- to_comp %>% 
  mutate(psampleid_n = as.numeric(psampleid_n)) %>% 
  select(title, psampleid = psampleid_n, pair) %>% 
  left_join(all_data %>% filter(layer == "epi")) %>% 
  select(title, psampleid, pair, all_of(taxa_comp))

# Day data
to_comp_d <- to_comp %>% 
  mutate(psampleid_d = as.numeric(psampleid_d)) %>% 
  select(title, psampleid = psampleid_d, pair) %>% 
  left_join(all_data %>% filter(layer == "epi")) %>% 
  select(title, psampleid, pair, all_of(taxa_comp))


# Test if there is a difference in concentrations
# Wilcoxon–Mann–Whitney test because concentrations are far from a normal distribution
p_val <- c() # initiate empty vector to store p-values
for (taxa in taxa_comp){ # Loop over taxa
  
  test_res <- wilcox.test(x = to_comp_n %>% pull(taxa), y = to_comp_d %>% pull(taxa), paired = TRUE) # perform test
  
  p_val <- c(p_val, test_res$p.value) # store p_value
}

# Store test p-values and taxa in a table
p_vals <- tibble(
  taxon = taxa_comp,
  p_val = p_val
) %>% 
  mutate( # compute significance for various thresholds
    sig_05 = p_val < 0.05, # *
    sig_01 = p_val < 0.01, # **
    sig_001 = p_val < 0.001, # ***
  )

p_vals

# Retirer les zéros
# GLM par paire. family quasi-poisson après avoir transformé les données en x10000. La plus petite valeur qui n’est pas zéro doit être à 1. 


# Boxplot of concentrations

to_comp_long <- to_comp_n %>% 
  mutate(day_night = "night") %>% 
  #pivot_longer(Trichodesmium:Acantharea, names_to = "taxon", values_to = "conc") 
  bind_rows(to_comp_d %>% mutate(day_night = "day")) %>% 
  pivot_longer(Trichodesmium:Acantharea, names_to = "taxon", values_to = "conc")


to_comp_long %>% 
  mutate(conc = conc * 1000) %>% 
  ggplot() +
  geom_boxplot(aes(x = taxon, y = conc, color = day_night)) +
  scale_color_manual(values = c("gray60", "black")) +
  # y scale as log
  scale_y_continuous(trans = "log1p", limits = quantile(to_comp_n$conc, c(0.1, 0.9))) +
  theme_classic() +
  labs(x = "Taxonomic group", y = expression(paste("Concentration (ind.", L^{-1}, ")")), color = "")



## Remove double zeros
to_comp_long <- to_comp_long %>% 
  select(-psampleid) %>% 
  pivot_wider(names_from = day_night, values_from = conc) %>% 
  filter(night + day > 0) 

# Stat test
p_vals <- to_comp_long %>% 
  arrange(taxon) %>% 
  group_by(taxon) %>% 
  summarise(p_val = wilcox.test(x = night, y = day, paired = TRUE)$p.value)


# Plot boxplots
to_comp_long %>% 
  pivot_longer(night:day, names_to = "day_night", values_to = "conc") %>% 
  mutate(conc = conc * 1000) %>% 
  ggplot() +
  geom_boxplot(aes(x = taxon, y = conc, color = day_night)) +
  scale_color_manual(values = c("gray60", "black")) +
  # y scale as log
  scale_y_continuous(trans = "log1p") +
  theme_classic() +
  labs(x = "Taxonomic group", y = expression(paste("Concentration (ind.", m^{-3}, ")")), color = "")
ggsave("plots/zoo/11.ww_conc_day_night.png")  
ggsave(file = "plots/paper/11.ww_conc_day_night.pdf", width = 166, height = 100, unit = "mm", dpi = 300)




## Concentrations at day and night for CCE ----
#--------------------------------------------------------------------------#
# Get concentrations at CCE in epipelagic layer
cce_epi <- all_data %>% 
  filter(layer == "epi") %>% # concentrations in epipelagic layer
  filter(str_detect(title, "CCELTER")) %>% # at CCELTER
  select(psampleid, datetime, day_night, Acantharea:Trichodesmium)

# Keep one version for Hellinger transformation
cce_epi_meta <- cce_epi %>% select(psampleid:day_night)
cce_epi_hel <- tibble(decostand(cce_epi %>% select(Acantharea:Trichodesmium), "hellinger"))
cce_epi_hel <- bind_cols(cce_epi_meta, cce_epi_hel)

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

cce_epi_hel <- cce_epi_hel %>% 
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
  mutate( # compute significance for various thresholds
    sig_05 = p_val < 0.05, # *
    sig_01 = p_val < 0.01, # **
    sig_001 = p_val < 0.001, # ***
    )
cce_p_vals
# Concentrations are different between day and night for Copepoda and Eumalacostraca, and maybe for Phaeodaria. 


# Boxplot of concentrations
cce_epi %>% 
  ggplot() +
  geom_boxplot(aes(x = taxon, y = conc, color = day_night), size = 0.3) +
  scale_color_manual(values = c("gray60", "black")) +
  # y scale as log
  scale_y_continuous(trans = "log1p") +
  theme_classic() +
  labs(x = "Taxonomic group", y = expression(paste("Concentration (ind.", L^{-1}, ")")), color = "")
ggsave("plots/zoo/11.cce_conc_day_night.png")  
ggsave(file = "plots/paper/11.cce_conc_day_night.pdf", width = 166, height = 100, unit = "mm", dpi = 300)


# Do the same with Hellinger-transformed data
# Test if there is a difference in concentrations
# Wilcoxon–Mann–Whitney test because concentrations are far from a normal distribution
p_val <- c() # initiate empty vector to store p-values
for (taxa in cce_taxa$taxon){ # Loop over taxa
  cce_day <- cce_epi_hel %>% filter(taxon == taxa) %>% filter(day_night == "day") # extract day concentrations
  cce_night <- cce_epi_hel %>% filter(taxon == taxa) %>% filter(day_night == "night") # extract night concentrations
  
  test_res <- wilcox.test(x = cce_day$conc, y = cce_night$conc) # perform test
  
  p_val <- c(p_val, test_res$p.value) # store p_value
}

# Store test p-values and taxa in a table
cce_hel_p_vals <- tibble(
  taxon = cce_taxa$taxon,
  p_val = p_val
) %>% 
  mutate( # compute significance for various thresholds
    sig_05 = p_val < 0.05, # *
    sig_01 = p_val < 0.01, # **
    sig_001 = p_val < 0.001, # ***
  )
cce_hel_p_vals


# Check high value of copepoda at day
cce_psample <- cce_epi %>% 
  filter(day_night == "day" & taxon == "Copepoda") %>% 
  arrange(desc(conc)) %>% 
  slice(1) %>% 
  pull(psampleid)

cce_prof <- all_data %>% filter(psampleid == cce_psample) %>% select(title:datetime) %>% pull(profile)
zoo %>% 
  select(title:profile, lat, lon, depth, taxon, abund) %>% 
  filter(profile == cce_prof & taxon == "Copepoda") %>% 
  filter(depth < 200) %>% 
  arrange(depth) %>% 
  summarise(sum(abund))
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

cols_epi <- c("#FFBF1A", "#53D3D3", "#10698D", "#123249", "#BF6322")

anosim_epi %>% 
  mutate(
    full = str_replace_all(full, " ", "\n"),
    full = factor(full, levels = full)
    ) %>% 
  filter(!short %in% c("day_night", "prod")) %>% 
  ggplot() +
  geom_col(aes(x = full, y = prop_exp, fill = full), show.legend = FALSE) +
  scale_fill_manual(values = cols_epi) +
  labs(x = "Regionalisation", y = "Explainable variance") +
  scale_y_continuous(expand = c(0,0), labels = scales::percent) +
  theme_classic() +
  theme(text = element_text(size = 18))


cols_meso <- c("#FFBF1A", "#53D3D3", "#50571F", "#10698D", "#123249", "#BF6322")

anosim_meso %>% 
  mutate(
    full = str_replace_all(full, " ", "\n"),
    full = factor(full, levels = full)
  ) %>% 
  mutate(full = factor(full, levels = full)) %>% 
  filter(!short %in% c("day_night", "prod")) %>% 
  ggplot() +
  geom_col(aes(x = full, y = prop_exp, fill = full), show.legend = FALSE) +
  scale_fill_manual(values = cols_meso) +
  labs(x = "Regionalisation", y = "Explainable variance") +
  scale_y_continuous(expand = c(0,0), labels = scales::percent) +
  theme_classic() +
  theme(text = element_text(size = 18))

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


## Eigenvalues and layers ----
#--------------------------------------------------------------------------#
# Investigate structure level of layers by comparing eigenvalues
# Do this with the same number of profiles in each layer
# Count profiles per layer
all_data %>% count(layer) %>% arrange(n)
# 827 profiles in bathypelagic layer
# Run a PCA on a subsample of plankton data
n_prof <- all_data %>% count(layer) %>% filter(n == min(n)) %>% pull(n)

layer_names <- c("epi", "mesosup", "mesoinf", "bathy") # layers to run PCA on
eig <- tibble() # initiate empty tibble to store results

for (my_layer in layer_names){ # loop over layers
  # Extract plankton data
  zoo <- all_data %>% 
    filter(layer == my_layer) %>% # relevant layer
    select(Acantharea:Trichodesmium) %>% # select plankton concentrations
    sample_n(n_prof) # sample <n_prof> profiles 
  
  # Apply Hellinger transformation
  zoo_hel <- tibble(decostand(zoo, "hellinger"))
  
  # Run PCA
  pca <- rda(zoo_hel)
  
  # Extract eigenvalues for this layer
  eig_layer <- as.data.frame(t(summary(eigenvals(pca)))) %>% 
    rownames_to_column(var = "PC") %>% 
    as_tibble() %>% 
    rename(eigenvalue = Eigenvalue, prop_exp = `Proportion Explained`, cum_prop = `Cumulative Proportion`) %>% 
    mutate(
      PC = str_remove(PC, "PC") %>% as.numeric() %>% as.factor(),
      eig_cor = eigenvalue - mean(eigenvalue), # compute centered eigenvalue
      layer = my_layer, # add a column with layer
      n = row_number(),
    )
  
  # Bind rows with eigvals from other layers
  eig <- bind_rows(eig, eig_layer)
}

eig %>% 
  ggplot() +
  geom_path(aes(x = n, y = eig_cor, colour = layer)) +
  # dotted black line at y=0
  geom_hline(aes(yintercept = 0), colour = "black", alpha = 0.5, lty = 3) +
  theme_classic() +
  theme(legend.position = c(.8,.7), legend.justification=c(0.5,0.5)) +
  scale_x_continuous(breaks = c(1:28)) +
  labs(colour = "Layer", x = "Principal componant", y = "Centered eigenvalue") +
  scale_color_manual(
    values = c("#a1dab4", "#41b6c4", "#2c7fb8", "#253494"), 
    name = "Layer", 
    labels = c("Epipelagic", "Mesopelagic sup", "Mesopelagic inf", 'Bathypelagic')
  )
ggsave("plots/zoo/11.eigenvalues_layers.png")  


## Depth of epi-meso pelagic boundary ----
#--------------------------------------------------------------------------#

summary(layers$epi)
layers %>% 
  ggplot() +
  geom_histogram(aes(x = epi), binwidth = 5, fill = "gray50") +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  labs(x = "Depth of epipelagic layer (m)", y = "Count") +
  theme_classic()

toto <- layers %>% 
  mutate(lat_bin = cut(abs(lat), breaks = c(0, 30, 60, 90))) %>% 
  filter(is.na(lat_bin))

layers %>% 
  mutate(lat_bin = cut(abs(lat), breaks = c(0, 30, 60, 90), right = FALSE)) %>% 
  ggplot() +
  geom_histogram(aes(x = epi, fill = lat_bin), binwidth = 5) +
  #scale_fill_viridis_d() +
  scale_fill_brewer(palette = "Blues", direction = -1) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  labs(x = "Depth of epipelagic layer (m)", y = "Count", fill = "Absolute \nlatitude") +
  theme_classic()
ggsave("plots/zoo/11.hist_epi.png")  
ggsave(file = "plots/paper/11.hist_epi.pdf", width = 166, height = 80, unit = "mm", dpi = 300)


layers %>% 
  ggplot() +
  geom_polygon(aes(x = lon, y = lat, group = group), fill = "gray", data = world) +
  geom_point(aes(x = lon, y = lat, color = epi), size= 0.8) +
  #scale_colour_viridis_c(option = "G") +
  scale_colour_cmocean(name = "deep") +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  labs(x = "Longitude", y = "Latitude", color ="Depth of \nepipelagic \nlayer (m)") +
  guides(colour = guide_colourbar(reverse = TRUE)) +
  theme_minimal() +
  coord_quickmap()
ggsave("plots/zoo/11.map_epi.png")  
ggsave(file = "plots/paper/11.map_epi.pdf", width = 166, height = 123, unit = "mm", dpi = 300)

# epipelagic layer deeper than 150 m in tropics
layers %>% 
  filter(between(lat, -30, 30)) %>% 
  mutate(deep = epi > 150) %>% 
  summarise(n = sum(deep) / length(deep))


## Nostocales ----
#--------------------------------------------------------------------------#
all_data %>% 
  filter(layer == "epi") %>% 
  ggplot() +
  geom_point(aes(x = lon, y = lat, color = Nostocales > 0)) +
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "gray") +
  #scale_color_viridis_c() +
  coord_quickmap()

all_data %>% 
  filter(layer == "epi") %>% 
  filter(Nostocales > 0)


## Epipelagic and mesopelagic plankton clusters ----
#--------------------------------------------------------------------------#
load("data/07.epi_profiles_cluster.Rdata")
ind_clust_epi <- ind_clust %>% mutate(clust_zoo = as.character(clust_zoo))

load("data/08.mesosup_profiles_cluster.Rdata")
ind_clust_mesosup <- ind_clust %>% mutate(clust_zoo = as.character(clust_zoo))

ind_clust <- ind_clust_epi %>% 
  select(title:psampleid, epi_clust = clust_zoo) %>% 
  left_join(ind_clust_mesosup %>% select(title:psampleid, meso_clust = clust_zoo))

ind_clust %>% 
  ggplot() +
  geom_point(aes(x = epi_clust, y = meso_clust)) +
  geom_jitter(aes(x = epi_clust, y = meso_clust))

#
#df <- mtcars %>%
#  make_long(cyl, vs, am, gear, carb)
#
#ggplot(df, aes(x = x, 
#               next_x = next_x, 
#               node = node, 
#               next_node = next_node,
#               fill = factor(node))) +
#  geom_sankey()
#
#ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
#  geom_sankey(flow.alpha = .6,
#              node.color = "gray30") +
#  geom_sankey_label(size = 3, color = "white", fill = "gray40") +
#  scale_fill_viridis_d() +
#  theme_sankey(base_size = 18) +
#  labs(x = NULL) +
#  theme(legend.position = "none",
#        plot.title = element_text(hjust = .5)) +
#  ggtitle("Car features")

ind_clust %>% 
  replace_na(list(meso_clust = "NA")) %>% 
  rename(epipelagic = epi_clust, mesopelagic = meso_clust) %>% 
  make_long(epipelagic, mesopelagic) %>% 
  ggplot(aes(x = x, next_x = next_x, node = node, next_node = next_node, label = node, fill = node)) +
  geom_sankey(flow.fill = "gray", flow.alpha = .6, node.color = "gray30") +
  geom_sankey_label(size = 3.5, color = "white", fill = "gray40") +
  scale_fill_viridis_d() +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none")

ggsave("plots/zoo/11.connec_epi_meso.png")  
ggsave(file = "plots/paper/11.connec_epi_meso.pdf", width = 166, height = 200, unit = "mm", dpi = 300)

ind_clust %>% count(epi_clust)
ind_clust %>% count(meso_clust)

ind_clust %>% 
  replace_na(list(meso_clust = "NA")) %>% 
  rename(epipelagic = epi_clust, mesopelagic = meso_clust) %>% 
  count(epipelagic, mesopelagic) %>% 
  group_by(epipelagic) %>% 
  mutate(n = n / sum(n))

ind_clust %>% 
  replace_na(list(meso_clust = "NA")) %>% 
  rename(epipelagic = epi_clust, mesopelagic = meso_clust) %>% 
  count(mesopelagic, epipelagic) %>% 
  group_by(mesopelagic) %>% 
  mutate(n = n / sum(n))

ind_clust %>% 
  #replace_na(list(meso_clust = "NA")) %>% 
  drop_na(meso_clust) %>% 
  rename(epipelagic = epi_clust, mesopelagic = meso_clust) %>% 
  count(epipelagic, mesopelagic) %>% 
  group_by(epipelagic) %>% 
  mutate(n = n / sum(n))

# Find mixed-type in epi and copepod in meso
ind_clust %>% 
  filter(epi_clust == "1" & meso_clust == "2") %>% 
  ggplot() +
  geom_polygon(aes(x = lon, y = lat, group = group), fill = "gray", data = world) +
  geom_point(aes(x = lon, y = lat)) +
  coord_quickmap()


ind_clust %>% 
  drop_na(meso_clust) %>% 
  mutate(
    epi_clust = ifelse(epi_clust == "1", "rhizaria", ifelse(epi_clust == "2", "copepod", "trichodesmium")),
    meso_clust = ifelse(meso_clust == "1", "copepod", ifelse(meso_clust == "2", "mixed", "phaeodaria")),
    combination = paste(epi_clust, meso_clust, sep = " // ")
    ) %>% 
  ggplot() +
  geom_polygon(aes(x = lon, y = lat, group = group), fill = "gray", data = world) +
  geom_point(aes(x = lon, y = lat), alpha = 0.2) +
  coord_quickmap() +
  facet_wrap(~combination)


## Uncommon taxa detected with chi-square ----
#--------------------------------------------------------------------------#

prop_zero <- function(x){sum(x == 0)/length(x)}

all_data %>% 
  filter(layer == "epi") %>% 
  select(taxon = Collodaria_colonial) %>% 
  transmute(presence = taxon > 0) %>% 
  sum()


## High Tricho concentrations ----
#--------------------------------------------------------------------------#
all_data %>% 
  filter(layer == "epi") %>% 
  filter(Trichodesmium > 0) %>% 
  ggplot() +
  geom_density(aes(x = Trichodesmium)) +
  scale_y_continuous(trans = "log1p", expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x = "Trichodesmium (L⁻¹)") +
  theme_minimal()
