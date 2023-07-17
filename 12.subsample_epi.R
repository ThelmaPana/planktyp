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


cor.test(proj$PC1, proj$sub_PC1)
cor.test(proj$PC2, proj$sub_PC2)

cor(proj$PC1, proj$sub_PC1) ^ 2
cor(proj$PC2, proj$sub_PC2) ^ 2

