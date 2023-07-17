#--------------------------------------------------------------------------#
# Project: planktyp
# Script purpose: Run analyses on a subsample of epipelagic profiles to limit oversampling effect of some areas
# Date: Thelma Panaïotis
# Author: 11/02/2020
#--------------------------------------------------------------------------#  

source("lib/set_up.R")
library(vegan)
library(castr)
library(ggrepel)
library(sp)
library(gstat)
library(pastecs)

load("data/06.all_data.Rdata")



# Set study layer
study_layer <- "epi"
subset <- all_data %>% filter(layer == study_layer)
message(nrow(subset), " profiles in epipelagic layer")

# Extract metadata
meta <- select(subset, title:gbp_code)

prec <- 10 # precision at which to round coordinates
thres <- 20 # max number of profile per <prec>*<prec> square


## Scale of autocorrelation ----
#--------------------------------------------------------------------------#
# Variogram for Copepoda / Trichodesmium / Collodaria
subset

hsp <- SpatialPointsDataFrame(
  # define the coordinates      
  coords=select(subset, lon, lat),
  # provide the associated data
  # (we need only SiO3 but we can include several other variables here)
  data=select(subset, Copepoda, Trichodesmium, Collodaria),
  # define the type of projection (required if coordinates are lon,lat)
  # here we specify that coordinates are longitudes and latitudes taken in the WGS84 datum
  # (this is true for most lon,lat data you will encounter)
  proj4string=CRS("+proj=longlat +datum=WGS84")
)
plot(hsp)

v1 <- variogram(Copepoda ~ 1, data=hsp, alpha = c(0,90)) # same variation in both direction
v1 <- variogram(Copepoda ~ 1, data=hsp, width = 50)
plot(v1)
# Difficult to fit a variogram

show.vgms(max = 5000, nugget=0.004, sill=0.01, range=50, models=c("Wav", "Sph", "Exc"))

v2 <- variogram(Trichodesmium ~ 1, data=hsp, width = 200)
plot(v2)

show.vgms(max = 5000, nugget=0.005, sill=0.04, range=1500, models=c("Wav", "Sph", "Exc"))
v_fitted <- fit.variogram(v2, vgm(model="Wav", psill=0.04, nugget=0.005, range=900))
plot(v2, v_fitted)


v3 <- variogram(Collodaria ~ 1, data=hsp, width = 200)
plot(v3)
v_fitted <- fit.variogram(v3, vgm(model="Sph", psill=0.00001, nugget=0.0001, range=1000))
plot(v3, v_fitted)


# Multivariate: self-Mahalanobis distance
?disto
# disto is for a timeseries




## Round lat and lon to compute profile density ----
#--------------------------------------------------------------------------#
prof_dens <- meta %>% 
  mutate( # round lon and lat
    lon_r = roundp(lon, prec, floor) + prec/2,
    lat_r = roundp(lat, prec, floor) + prec/2
    ) %>% 
  group_by(lon_r, lat_r) %>% # group by square
  summarise(n = n()) %>% # compute number of profiles
  ungroup()

message(nrow(prof_dens %>% filter(n>thres)), " squares with more profiles than threshold")

# Plot map of profile density with colorbar
prof_dens %>% 
  ggplot() +
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "gray") +
  geom_tile(aes(x = lon_r, y = lat_r, fill = n)) +
  coord_quickmap() +
  scale_fill_viridis_c(trans="log10") + # scale color as log
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(x = "Longitude", y = "Latitude", fill = "Profile nb") 

# Plot map of profile density with threshold limit
prof_dens %>% 
  ggplot() +
  geom_polygon(data = world, aes(x = lon, y = lat, group = group), fill = "gray") +
  geom_tile(aes(x = lon_r, y = lat_r, fill = n > thres)) +
  coord_quickmap() +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(x = "Longitude", y = "Latitude", fill = "Profile nb > thres") 


## Subsample profiles in high density squares ----
#--------------------------------------------------------------------------#
sub_meta <- meta %>% 
  mutate( # round lon and lat
    lon_r = roundp(lon, prec, floor) + prec/2, 
    lat_r = roundp(lat, prec, floor) + prec/2
    ) %>% 
  group_by(lon_r, lat_r) %>% # group by square
  mutate(n=n()) %>% # compute number of profiles in square
  sample_n(min(n, thres)) %>% # sample at most <thres> profiles in each square
  ungroup()

message("Subset of ", nrow(sub_meta), " profiles from ", nrow(meta), " in epipelagic layer")


# Keep only selected profiles
subset <- subset %>% filter(psampleid %in% sub_meta$psampleid)


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
biplot(pca, display = c("sites", "species"), type = c("text", "points"))


## Predict discarded stations ----
#--------------------------------------------------------------------------#
meta_all <- all_data %>% filter(layer == study_layer) %>% select(title:datetime)
zoo_all <- all_data %>% filter(layer == study_layer) %>% select(Acantharea:Trichodesmium)
# Hellinger transformation
zoo_all_hel <- tibble(decostand(zoo_all, "hellinger"))


## S3 method for class 'rda'
sub_proj <- predict(pca, zoo_all_hel, type = "wa", scaling = 1) %>% 
  as_tibble() %>% 
  select(sub_PC1 = PC1, sub_PC2 = PC2)
sub_proj <- bind_cols(meta_all, sub_proj)
summary(sub_proj)


## PCA with all stations ----
#--------------------------------------------------------------------------#
pca_all <- rda(zoo_all_hel) # don't scale data because we use Hellinger transformation

proj <- vegan::scores(pca_all, display="sites", choices=c(1:2), scaling=1) %>% 
  as_tibble() %>% 
  bind_cols(meta_all, .)
summary(proj)


## Correlation between PCs ----
#--------------------------------------------------------------------------#
proj <- proj %>% left_join(sub_proj)

#hist(proj$PC1)
#hist(proj$PC2)
#hist(proj$sub_PC1)
#hist(proj$sub_PC2)

cor.test(proj$PC1, proj$sub_PC1)
cor.test(proj$PC2, proj$sub_PC2)

cor(proj$PC1, proj$sub_PC1) ^ 2
cor(proj$PC2, proj$sub_PC2) ^ 2


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
png(file = "plots/subsample/epi/12.plankton_dendrogram.png", width = 9.79, height = 7.96, units = "in", res = 300)
plot(clust_zoo, main = "Cluster dendrogram on plankton data in epipelagic layer")
# Choose number of clusters
nclust <- 3
# Plot clusters
rect.hclust(clust_zoo, k=nclust)
dev.off()
# Add clusters to table of individuals
ind_clust$clust_zoo <- as.factor(cutree(clust_zoo, k = nclust))

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
ggsave(file = "plots/subsample/epi/12.plankton_clusters_profile_nb.png")


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


## Plot PCA eigenvalues plot ----
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
ggsave(file = "plots/subsample/epi/12.plankton_pca_variance_explained.png")


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
  scale_color_brewer(palette = "Set2") +
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
  theme(legend.position = "bottom") +
  # Legend inside plot
  #theme(
  #  legend.position = c(.95, .95),
  #  legend.justification = c("right", "top"),
  #  legend.direction = "horizontal",
  #) +
  # Bigger points in the legend, legend title on top
  guides(colour = guide_legend(override.aes = list(size=2), title.position="top")) 
ggsave(file = "plots/subsample/epi/12.plankton_pca_biplot_1_2.png")


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
  scale_color_brewer(palette = "Set2") +
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
ggsave(file = "plots/subsample/epi/12.plankton_pca_biplot_2_3.png")


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
  scale_color_brewer(palette = "Set2") +
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
ggsave(file = "plots/subsample/epi/12.plankton_clusters_map.png")


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
ggsave(file = "plots/subsample/epi/12.plankton_clusters_composition.png")

message("Same pattern after subsampling than before")
