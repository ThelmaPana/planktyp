#--------------------------------------------------------------------------#
# Project: planktyp
# Script purpose: Run analyses for epipelagic layer
# Date: 03/02/2021
# Author: Thelma Panaiotis
#--------------------------------------------------------------------------#

# TODO explain variables selection on PCA plots

source("lib/set_up.R")
library(vegan)
library(ggrepel)
library(grid)
library(gridExtra)

load("data/06.all_data.Rdata")

# Set study layer
study_layer <- "epi"
subset <- all_data %>% filter(layer == study_layer)
message(nrow(subset), " profiles in epipelagic layer")

# Minimum number of profiles per regionalisation modality
n_min <- 50
n_target_max <- 13
n_target_env <- 15
n_target_null <- 12

#n_min <- 25
#n_target_max <- 15
#n_target_env <- 15
#n_target_null <- 15

## Create a nice color scale for the 3 groups
# - group 1 will be collodaria --> light green #a6d854 (we'll use dark green for phaeodaria in meso)
# - group 2 will be copepods --> red/orange #fc8d62
# - group 3 will be tricho --> purple #756bb1 
my_colors <- c("#a6d854", "#fc8d62", "#756bb1")


## Prepare data ----
#--------------------------------------------------------------------------#
# Extract metadata
meta <- select(subset, title:gbp_code)

# Extract env data
env <- select(subset, temp:kd490) %>% 
  # replace missing data by variable average
  replace_na(as.list(colMeans(.,na.rm=T))) %>% 
  # log transform snow_conc and bulk_conc to approach normality
  mutate_at(vars(snow_conc, bulk_conc), list(log1p)) %>% 
  # standardize data to mean = 0 and sd = 1
  transmute_all(list(scale2))


# Extract zooplankton
zoo <- select(subset, Acantharea:Trichodesmium)
# Checj classes with only zeros
cat("Groups with no organisms for this layer:", zoo[colSums(zoo) == 0] %>% names)

# Hellinger transformation
zoo_hel <- tibble(decostand(zoo, "hellinger"))


## Run PCA ----
#--------------------------------------------------------------------------#
# Run PCA
pca <- rda(zoo_hel) # don't scale data because we use Hellinger transformation

# Fit standardized env data on PCA axes
ef <- envfit (pca ~ ., data = env, perm = 999, na.rm = T, choices=c(1:5))

# Extract coordinates of env data projection
env_proj <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>% 
  rownames_to_column(var = "variable") %>% 
  as_tibble() %>% 
  mutate(orig = 0)


## Clustering of zoo profiles ----
#--------------------------------------------------------------------------#
# Extract PCA scores for clustering using scaling 1 
ind_clust <- vegan::scores(pca, display="sites", choices=c(1:5), scaling=1) %>% 
  as_tibble() %>% 
  bind_cols(meta, .)

# Compute euclidean distance
euc_dist <- dist(select(ind_clust, c(PC1:PC5)), method = "euclidian")
# Compute hierarchical clustering with Ward method
clust_zoo <- hclust(euc_dist, method = "ward.D2")
# Plot dendrogram
png(file = "plots/analysis/epi/07.plankton_dendrogram.png", width = 9.79, height = 7.96, units = "in", res = 300)
plot(clust_zoo, main = "Cluster dendrogram on plankton data in epipelagic layer")
# Choose number of clusters
nclust <- 3
# Plot clusters
rect.hclust(clust_zoo, k=nclust, border=my_colors)
dev.off()
# Add clusters to table of individuals
ind_clust$clust_zoo <- as.factor(cutree(clust_zoo, k = nclust))

# Save plot for paper
svg(file = "plots/paper/07.plankton_dendrogram.svg", width = 9.79, height = 7.96)
plot(clust_zoo, main = "Cluster dendrogram on plankton data in epipelagic layer")
nclust <- 3
rect.hclust(clust_zoo, k=nclust, border=my_colors)
dev.off()

# Check number of profile per cluster
ind_clust %>% 
  group_by(clust_zoo) %>% 
  summarise(nb = n()) %>% 
  ungroup() %>% 
  ggplot() +
  geom_col(aes(x = clust_zoo, y = nb, fill = nb>50)) +
  ggtitle("Number of profiles per plankton cluster in epi layer") +
  theme_classic() +
  # y scale as a percentage and remove white space
  scale_y_continuous(expand = c(0,0))
ggsave(file = "plots/analysis/epi/07.plankton_clusters_profile_nb.png")


## Prepare PCA data for plots ----
#--------------------------------------------------------------------------#
# Extract eigenvalues
eig <- as.data.frame(t(summary(eigenvals(pca)))) %>% 
  rownames_to_column(var = "PC") %>% 
  as_tibble() %>% 
  rename(eigenvalue = Eigenvalue, prop_exp = `Proportion Explained`, cum_prop = `Cumulative Proportion`) %>% 
  mutate(PC = str_remove(PC, "PC") %>% as.numeric() %>% as.factor())

# Variance explained by first 2 PC
eig %>% slice(1:2) %>% summarise_if(is.numeric, sum)
# and by first 4 PC
eig %>% slice(1:4) %>% summarise_if(is.numeric, sum)

# Extract PCA scores of profiles for plots (scaling 2)
ind_plot <- vegan::scores(pca, display="sites", choices=c(1:5), scaling=2) %>% 
  as_tibble() %>% 
  bind_cols(meta, .) %>% 
  left_join(ind_clust %>% select(-c(PC1:PC5))) # ignore coordinates in scaling 1 

# Select variables to plot with scaling 1
# Only variables whose contribution is higher than equivalent contribution will be plotted
var_sel <- vegan::scores(pca, display="species", choices=c(1:5), scaling=1) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "variable") %>% 
  as_tibble() %>% 
  mutate(
    p12 = (PC1^2+PC2^2)^0.5,
    p23 = (PC2^2+PC3^2)^0.5,
    p34 = (PC3^2+PC4^2)^0.5,
    const = ((nrow(ind_plot)-1)*sum(eig$eigenvalue))^0.25,
    r = (2/nrow(eig))^0.5*const,
    disp12 = p12 > r,
    disp23 = p23 > r,
    disp34 = p34 > r,
  )  %>% 
  select(variable, disp12, disp23, disp34)

# Extract variables scores in scaling 2
var_plot <- vegan::scores(pca, display="species", choices=c(1:5), scaling=2) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "variable") %>% 
  as_tibble() %>% 
  mutate(orig = 0) %>% 
  left_join(var_sel)


## Plot PCA eigenvalues ----
#--------------------------------------------------------------------------#
eig %>% 
  ggplot() +
  # cols for variance of each PC
  geom_col(aes(x = PC, y = prop_exp)) +
  # line for cumulative variance
  geom_line(aes(x = PC, y = cum_prop, group=1)) +
  # rename axes
  xlab("Principal componant") +
  ylab("Explained variance") +
  ggtitle("PC and explained variance of plankton PCA in epi layer") +
  theme_classic() +
  # y scale as a percentage and remove white space
  scale_y_continuous(labels=scales::percent, expand = c(0,0))
ggsave(file = "plots/analysis/epi/07.plankton_pca_variance_explained.png")


## Biplot PCA axes 1-2 ----
#--------------------------------------------------------------------------#
# Scaling factor for plots
k <- 10

# circle of equivalent contributions
circ <- circle_fun(c(0,0), 2*sqrt(2/nrow(eig)), npoints = 500)

#var_plot %>% filter(disp12)

ggplot() +
  # Vertical and horizontal lines
  geom_vline(xintercept = 0, color = "gray") +
  geom_hline(yintercept = 0, color = "gray") +
  # Circle of equivalent contributions
  #geom_path(aes(k*x, k*y), data = circ, color = "gray") +
  # Points for individuals (i.e. profiles) with color according to plankton cluster
  geom_point(aes(x = k*PC1, y = k*PC2, colour = clust_zoo), data = ind_plot, alpha = 0.5, size = 0.5) +
  # Points for species with high contribution
  geom_point(aes(x = PC1, y = PC2), data = filter(var_plot, disp12)) +
  # Labels for species with high contribution
  geom_label_repel(aes(x = PC1, y = PC2, label = variable), data = filter(var_plot, disp12), fill="white") +
  # Segments for projected environmental data
  geom_segment(aes(x = orig, y = orig, xend = k*PC1, yend = k*PC2), data = env_proj, colour = "grey40", arrow = arrow(length = unit(0.2, "cm"))) +
  # Labels for projected environmental data
  geom_text(aes(x = k*PC1, y = k*PC2, label = variable), data = env_proj, colour = "grey40", hjust=ifelse(env_proj$PC1>0, -0.05, 1.05)) +
  # Fixed ratio between axes
  coord_fixed() +
  # Add a little padding
  scale_x_continuous(expand = c(.1, 0)) +
  # Nice color palette
  #scale_color_brewer(palette = "Set2") +
  scale_color_manual(values = my_colors) +
  # Put explained variance in axes names
  xlab(paste0("PC 1", " (", format(round(100*eig$prop_exp[1], 1), nsmall = 1), "%)")) +
  ylab(paste0("PC 2", " (", format(round(100*eig$prop_exp[2], 1), nsmall = 1), "%)")) +
  # Legend title
  labs(color = "Plankton clusters") +
  # Plot title
  ggtitle("PCA biplot of epipelagic profiles, PC1 and PC2") +
  # Nice theme
  theme_classic() +
  # Legend position at the bottom
  #theme(legend.position = "bottom") +
  # Legend inside plot
  theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.direction = "horizontal",
  ) +
  # Bigger points in the legend, legend title on top
  guides(colour = guide_legend(override.aes = list(size=2), title.position="top")) 
ggsave(file = "plots/analysis/epi/07.plankton_pca_biplot_1_2.png")
ggsave(file = "plots/paper/07.plankton_pca_biplot_1_2.svg")


## Biplot PCA axes 2-3 ----
#--------------------------------------------------------------------------#
# Scaling factor for plots
k <- 10

# circle of equivalent contributions
circ <- circle_fun(c(0,0), 2*sqrt(2/nrow(eig)), npoints = 500)

#var_plot %>% filter(disp12)

ggplot() +
  # Vertical and horizontal lines
  geom_vline(xintercept = 0, color = "gray") +
  geom_hline(yintercept = 0, color = "gray") +
  # Circle of equivalent contributions
  #geom_path(aes(k*x, k*y), data = circ, color = "gray") +
  # Points for individuals (i.e. profiles) with color according to plankton cluster
  geom_point(aes(x = k*PC2, y = k*PC3, colour = clust_zoo), data = ind_plot, alpha = 0.5, size = 0.5) +
  # Points for species with high contribution
  geom_point(aes(x = PC2, y = PC3), data = filter(var_plot, disp23)) +
  # Labels for species with high contribution
  geom_label_repel(aes(x = PC2, y = PC3, label = variable), data = filter(var_plot, disp23), fill="white") +
  # Segments for projected environmental data
  geom_segment(aes(x = orig, y = orig, xend = k*PC2, yend = k*PC3), data = env_proj, colour = "grey40", arrow = arrow(length = unit(0.2, "cm"))) +
  # Labels for projected environmental data
  geom_text(aes(x = k*PC2, y = k*PC3, label = variable), data = env_proj, colour = "grey40", hjust=ifelse(env_proj$PC2>0, -0.05, 1.05)) +
  # Fixed ratio between axes
  coord_fixed() +
  # Add a little padding
  scale_x_continuous(expand = c(.1, 0)) +
  # Nice color palette
  #scale_color_brewer(palette = "Set2") +
  scale_color_manual(values = my_colors) +
  # Put explained variance in axes names
  xlab(paste0("PC 2", " (", format(round(100*eig$prop_exp[2], 1), nsmall = 1), "%)")) +
  ylab(paste0("PC 3", " (", format(round(100*eig$prop_exp[3], 1), nsmall = 1), "%)")) +
  # Legend title
  labs(color = "Plankton clusters") +
  # Plot title
  ggtitle("PCA biplot of epipelagic profiles, PC2 and PC3") +
  # Nice theme
  theme_classic() +
  # Legend position at the bottom
  theme(legend.position = "bottom") +
  # Bigger points in the legend
  guides(colour = guide_legend(override.aes = list(size=3)))
ggsave(file = "plots/analysis/epi/07.plankton_pca_biplot_2_3.png")


## Plot map of plankton clusters ----
#--------------------------------------------------------------------------#
ind_plot %>% 
  ggplot() +
  # Map background
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "gray") +
  # Map coordinates
  coord_quickmap() +
  # Points for profiles coordinates with color according to plankton cluster
  geom_point(aes(x = lon, y = lat, color = clust_zoo), size = 1) +
  # Nice theme
  theme_minimal() +
  # Use same color palette as for PCA biplot
  #scale_color_brewer(palette = "Set2") +
  scale_color_manual(values = my_colors) +
  # Rename axes
  labs(x = "Longitude", y = "Latitude", color = "Plankton clusters") +
  # Legend at the bottom
  #theme(legend.position="bottom") +
  # Legend inside plot
  theme(
    legend.position = c(0.025, 0.05),
    legend.justification = c("left", "bottom"),
    legend.background = element_rect(fill="white", size=0, linetype="solid"),
    legend.direction = "horizontal",
  ) +
  # Remove white spaces on the border
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  # Bigger points in the legend, legend title on top
  guides(colour = guide_legend(override.aes = list(size=2), title.position="top")) +
  ggtitle("Map of epipelagic plankton clusters")
ggsave(file = "plots/analysis/epi/07.plankton_clusters_map.png")
ggsave(file = "plots/paper/07.plankton_clusters_map.svg")


## Plot clusters composition ----
#--------------------------------------------------------------------------#
# In each cluster, compute mean concentration for each taxa
# Then compute proportion of each taxon per cluster

comp <- bind_cols(zoo, ind_clust) %>% 
  select(clust_zoo, Acantharea:Trichodesmium) %>% 
  gather(Acantharea:Trichodesmium, key = "taxon", value = "conc") %>% 
  # Group by plankton cluster and taxon and compute mean concentration
  group_by(clust_zoo, taxon) %>% 
  summarise(conc = mean(conc)) %>% 
  ungroup() %>% 
  # Group by cluster
  group_by(clust_zoo) %>% 
  # Compute the relative proportion of each taxa in the cluster
  mutate(prop = conc / sum(conc)) %>% 
  ungroup() %>% 
  # Make taxon a factor for the plot
  mutate(taxon = factor(taxon, levels = unique(taxon)))

comp %>% 
  ggplot() +
  # Plankton clusters on x axis, taxon on y axis, fill with proportion
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
  ggtitle("Composition of epipelagic plankton clusters")
ggsave(file = "plots/analysis/epi/07.plankton_clusters_composition.png")


## Compute number of modalities to aim for in partitionings ----
#--------------------------------------------------------------------------#
# Three groupings to generate:
# - maximal model on plankton data themselve
# - immediate environment from environmental data
# - null model with random clusters
# All these groupings should have ~ the same number of modalities as the ones that already exist.
# We already have 3 partitionings: Longhurst provinces, latitude bands and mbgcp (but this lats one is only for mesopelagic layer)
# Let's count modalities in these partitionings, but consider only modalities with more than 50 profiles.

# for latitude bands
count_b_lat <- meta %>% 
  select(psampleid, b_lat, long_code) %>% 
  count(b_lat) 
message(nrow(count_b_lat %>% filter(n > n_min)), " modalities for latitude bands")
count_b_lat %>% 
  ggplot() +
  geom_col(aes(x = b_lat, y = n, fill = n > n_min)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Number of profiles per latitude band")

# for longhurst provinces
count_long <- meta %>% 
  select(psampleid, b_lat, long_code) %>% 
  count(long_code)
message(nrow(count_long %>% filter(n > n_min)), " modalities for longhurst provinces")
count_long %>% 
  ggplot() +
  geom_col(aes(x = long_code, y = n, fill = n > n_min)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Number of profiles per longhurst province")


# These counts account for all profiles.
# Now exclude profiles which do not belong in selected modalities.
counts <- meta %>% 
  select(psampleid, b_lat, long_code) %>% 
  # ignore profiles not assigned to latitude band or Longhurst province
  drop_na(b_lat, long_code) %>% 
  # count profiles per latitude band
  add_count(b_lat, name = "n_b_lat") %>% 
  # count profiles per longhurst province
  add_count(long_code, name = "n_long") %>% 
  filter(n_b_lat > n_min & n_long > n_min)

counts %>% count(b_lat)
counts %>% count(long_code)
message("Only ", nrow(counts), " profiles left in selected modalities")
message("We should aim for ~ 12 modalities in other partitionings.")


## Generate clustering for maximal model ----
#--------------------------------------------------------------------------#
# Clustering for maximal model in generated with plankton data themselves.
# It shows the maximum part of variance that can be explained by a partitioning with n modalities.
# Use same clustering as before but change number of clusters

png(file = "plots/analysis/epi/07.plankton_dendrogram_max_model.png", width = 9.79, height = 7.96, units = "in", res = 300)
plot(clust_zoo, main = "Cluster dendrogram on plankton data in epipelagic layer for max model")
# Choose number of clusters
nclust <- n_target_max
# Plot clusters
rect.hclust(clust_zoo, k=nclust)
dev.off()
# Add clusters to table of individuals
ind_clust$mod_max <- as.factor(cutree(clust_zoo, k = nclust))

# And join to table of selected profiles
counts <- counts %>% left_join(ind_clust %>% select(psampleid, mod_max))
message(nrow(counts %>% count(mod_max) %>% filter(n > n_min)), " modalities for maximal model")

# Remove modalities with not enough profiles
counts <- counts %>% 
  add_count(mod_max) %>% 
  filter(n > n_min) %>% 
  select(-n)
message(nrow(counts), " profiles left")

## Generate clustering on immediate environment ----
#--------------------------------------------------------------------------#
# PCA on env data
pca_env <- rda(env) # no need to scale because we already did it

# Extract eigenvalues for env PCA
eig_env <- as.data.frame(t(summary(eigenvals(pca_env)))) %>% 
  rownames_to_column(var = "PC") %>% 
  as_tibble() %>% 
  rename(eigenvalue = Eigenvalue, prop_exp = `Proportion Explained`, cum_prop = `Cumulative Proportion`) %>% 
  mutate(PC = str_remove(PC, "PC") %>% as.numeric() %>% as.factor())

# Plot variance of env PCA
eig_env %>% 
  ggplot() +
  # cols for variance of each PC
  geom_col(aes(x = PC, y = prop_exp)) +
  # line for cumulative variance
  geom_line(aes(x = PC, y = cum_prop, group=1)) +
  # rename axes
  xlab("Principal componant") +
  ylab("Explained variance") +
  ggtitle("PC and explained variance of env data PCA in epi layer") +
  theme_classic() +
  # y scale as a percentage and remove white space
  scale_y_continuous(labels=scales::percent, expand = c(0,0))
ggsave(file = "plots/analysis/epi/07.env_pca_variance_explained.png")

# Biplot of env PCA
png(file = "plots/analysis/epi/07.env_pca_biplot_1_2.png", width = 9.79, height = 7.96, units = "in", res = 300)
biplot(pca_env, xlab = "", ylab = "",)
orditorp(pca_env, display="sp")
title(
  xlab = paste0("PC 1", " (", format(round(100*eig_env$prop_exp[1], 2), nsmall = 2), "%)"), 
  ylab = paste0("PC 2", " (", format(round(100*eig_env$prop_exp[2], 2), nsmall = 2), "%)"), 
  main = "PCA biplot of env data in epi layer"
  ) 
dev.off()

# Extract scores in scaling 1 for clustering on env data
# Use 5 firts PCs
env_ind_clust <- vegan::scores(pca_env, display="sites", choices=c(1:5), scaling=1) %>% 
  as_tibble() %>% 
  bind_cols(meta, .)

# Compute euclidean distance
env_euc_dist <- dist(select(env_ind_clust, c(PC1:PC5)), method = "euclidian")
# Compute hierarchical clustering with Ward method
clust_env <- hclust(env_euc_dist, method = "ward.D2")
# Plot dendrogram
png(file = "plots/analysis/epi/07.env_dendrogram.png", width = 9.79, height = 7.96, units = "in", res = 300)
plot(clust_env, main = "Cluster dendrogram env data epipelagic layer")
# Choose number of clusters
nclust <- n_target_env
# Plot clusters
rect.hclust(clust_env, k=nclust)
dev.off()
# Add clusters to table of individuals
ind_clust$clust_env <- as.factor(cutree(clust_env, k = nclust))

# And join to table of selected profiles
counts <- counts %>% left_join(ind_clust %>% select(psampleid, clust_env))
message(nrow(counts %>% count(clust_env) %>% filter(n > n_min)), " modalities for env clustering")

# Remove modalities with not enough profiles
counts <- counts %>% 
  add_count(clust_env) %>% 
  filter(n > n_min) %>% 
  select(-n)
message(nrow(counts), " profiles left")


## Generate random groups for a null model ----
#--------------------------------------------------------------------------#
# Compute a random group for null model (12 groups)
counts <- counts %>% 
  mutate(mod_null = factor(sample(c(1:n_target_null), nrow(.), replace = TRUE)))


## Plot number of profiles per partitioning modalities ----
#--------------------------------------------------------------------------#
# Get groups for selected profiles
groups <- counts %>% 
  select(psampleid, mod_max, b_lat, long_code, clust_env, mod_null) %>% 
  left_join(meta %>% select(psampleid, day_night, prod))

# Look at number of modalities for each partitioning
modalities <- groups %>% 
  select(mod_max:prod) %>% 
  summarise_all(list(count_mods))
modalities

p1 <- as.data.frame(t(modalities)) %>% 
  rownames_to_column(var = "partitioning") %>% 
  rename(n_mod = V1) %>% 
  ggplot() +
  geom_col(aes(x = partitioning, y = n_mod)) +
  coord_flip() +
  xlab("Number of modalities") +
  theme_light()

# Plot number of profiles for each modality of partitionings
p2 <- ggplot(groups) +
  geom_bar(aes(x = mod_max)) +
  ggtitle("Profile number in maximal model groups") +
  theme_light()

p3 <- ggplot(groups) +
  geom_bar(aes(x = b_lat)) +
  ggtitle("Profile number in latitude band groups") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p4 <- ggplot(groups) +
  geom_bar(aes(x = long_code)) +
  ggtitle("Profile number in longhurst provinces groups") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p5 <- ggplot(groups) +
  geom_bar(aes(x = clust_env)) +
  ggtitle("Profile number in env clustering groups") +
  theme_light()

p6 <- ggplot(groups) +
  geom_bar(aes(x = mod_null)) +
  ggtitle("Profile number in null model groups") +
  theme_light()

p7 <- ggplot(groups) +
  geom_bar(aes(x = day_night)) +
  ggtitle("Profile number at day and night") +
  theme_light()

p8 <- ggplot(groups) +
  geom_bar(aes(x = prod)) +
  ggtitle("Profile number at productive and non productive seasons") +
  theme_light()
g <- arrangeGrob(p1, p2, p3, p4, p5, p6, p7, p8, ncol=3)
ggsave(g, file = "plots/analysis/epi/07.partitioning_modalities.png")


## Perform anosim to compute variance explained by each partitioning ----
#--------------------------------------------------------------------------#
# Join groups and hellinger transformed zooplankton data
data_anosim <- groups %>% 
  left_join(bind_cols(select(meta, psampleid), zoo_hel)) # need to bind meta and zoo data because psampleid not in zoo data

# Make a separate table with plankton concentrations to make vegan happy
zoo_anosim <- data_anosim %>% select(Acantharea:Trichodesmium)

# List of partitionings to test
short_names <- c(
  "mod_max",
  "long_code",
  "b_lat",
  "clust_env",
  "mod_null",
  "day_night",
  "prod"
)
full_names <- c(
  "Maximal model",
  "Longhurst provinces",
  "Latitude bands",
  "Local environment",
  "Null model",
  "Circadian cycle",
  "Seasonal cycle"
)
# Put this in a table
models <- tibble(short = short_names, full = full_names)


# Initiate empty list to store results
list_mod <- c()
list_Radj <- c()
list_p <- c()

# Loop over partitionings
for (mod in models$short){
  # Run RDA
  rda_mod <- rda(formula(paste0('zoo_anosim~', mod)), data=data_anosim, scale = F)
  # Run anova
  anov <- anova(rda_mod, permutations=how(nperm=999))
  # Compute R square
  r_square <- RsquareAdj(rda_mod)
  
  # Store results in lists
  list_mod <- c(list_mod, anov$Df[1]+1)
  list_Radj <- c(list_Radj, r_square$adj.r.squared)
  list_p <- c(list_p, anov$`Pr(>F)`[1])
}

# Convert lists to table and compute part of explainable variance explained by each model
anosim_res <- 
  tibble(
  modalities = list_mod,
  Radj = list_Radj,
  p = list_p
  ) %>% 
  bind_cols(models, .) %>% 
  mutate(prop_exp = Radj / max(Radj))
anosim_res

# Plot results
anosim_res %>% 
  # ignore models of day/night and seasons
  filter(!(short %in% c("day_night", "prod"))) %>% 
  # sort data and reset factor for a nice plot
  arrange(desc(prop_exp)) %>% 
  mutate(full = factor(full, full)) %>% 
  ggplot() +
  geom_col(aes(x = full, y = prop_exp, fill = full), show.legend = FALSE) +
  # Nice theme
  theme_classic() +
  # y scale as a percentage and remove white space
  scale_y_continuous(labels=scales::percent, expand = c(0,0)) +
  # Axes labels
  xlab("Regionalisation model") +
  ylab("Explainable variance") +
  ggtitle("Part of plankton data variance explained by partitionings in epi layer")
ggsave(file = "plots/analysis/epi/07.partitioning_anosim.png")


## Save useful data ----
#--------------------------------------------------------------------------#
# For each saved table, generate a column with layer name

# Anosim results
anosim_res <- anosim_res %>% mutate(layer = study_layer) %>% 
save(anosim_res, file="data/07.epi_anosim.Rdata")

# Plankton composition of plankton clusters
comp <- comp %>% mutate(layer = study_layer)
save(comp, file="data/07.epi_plankton_cluster_comp.Rdata")

# Partitionings for anosim
groups <- groups %>% mutate(layer = study_layer)
save(groups, file="data/07.epi_partitionings.Rdata")

# Plankton PCA eigenvalues
eig <- eig %>% mutate(layer = study_layer)
save(eig, file="data/07.epi_plankton_eigvals.Rdata")
