#--------------------------------------------------------------------------#
# Project: plantyp
# Script purpose: Clean and process CTD data
# Date: 01/02/2021
# Author: Thelma Panaiotis
#--------------------------------------------------------------------------#

source("lib/set_up.R")
library(multidplyr)
library(parallel)
library(castr)
library(gsw)

load("data/02.ctd_raw.Rdata")


## Cure CTD data depths (disordered, up/down cast...) ----
#--------------------------------------------------------------------------#

## Remove problematic profiles
pb_profiles <- c(
  "p1604_10", # abnormal depth range: from -1 to -26982
  "tara_205_03" # abnormal depth range: from 57707 to 86380
)
ctd <- ctd_raw %>% filter(!profile %in% pb_profiles)
rm(ctd_raw)

## Correct depth
# remove all negative depth
ctd <- ctd %>% filter(depth > 0)

## Look for down/up cast and proceed to corrections to keep only dowcasts
ctd %>% 
  group_by(title, profile) %>% 
  mutate(row = row_number()) %>% 
  ungroup() %>% 
  ggplot() +
  geom_path(aes(x = row, y = -depth, group = profile, color = profile), show.legend = FALSE) +
  facet_wrap(~title, scales = "free") +
  ggsave("plots/env/03.depth_dist.png")

# There are three main cases of depth distribution for env data:
# - downcast and upcast -> keep downcast only
# - downcast only but disordered -> reorder depths
# - downcast only and ordered -> ok 
# Plus there are also a few profiles to check and fix by hand
# Split data into three groups to cure them separately


# Projects with downcast and upcast
down_up <- c(
  "UVP5 AN Arctique (an1304 an1405 an1407) part A to validate",
  "UVP5 BOUM 2008",
  "UVP5 CASSIOPEE 2015",
  "UVP5 CCELTER 2016",
  "UVP5 Dyfamed 2013 2014 2015 2016 2017 2018",
  "UVP5 GREEN EDGE Ice Camp 2015",
  "UVP5 LOHAFEX 2009",
  "UVP5 MooseGE 2012 2013 2014 2015 2016",
  "UVP5hd NAAMES02"
)
# Fix using slice for beginning to max depth
ctd_du <- ctd %>% 
  filter(title %in% down_up) %>% 
  group_by(title, profile, lat, lon, datetime, psampleid) %>% 
  slice(1:max(depth)) %>% 
  arrange(depth) %>% # fix for profile "lohafex_ps73_204", see below
  ungroup()


# Projects with disordered downcast
down_desord <- c(
  "UVP5 CCELTER (Iter2008 ccelter_2012 ccelter_2014)",
  "UVP5zd CCELTER 2011",
  "UVP5zd KEOPS"
)
# Fix using arrange
ctd_dd <- ctd %>% 
  filter(title %in% down_desord) %>% 
  group_by(title, profile, lat, lon, datetime, psampleid) %>% 
  arrange(depth)
  

# Projects which are ok, nor in down_cast neither in down_desord
ctd_ok <- ctd %>% 
  filter(!(title %in% down_up | title %in% down_desord))

# Bind rows from three groups after fix
ctd_fix <- bind_rows(ctd_du, ctd_dd, ctd_ok)
rm(ctd_du, ctd_dd, ctd_ok)

# Projects with profiles to check
to_check <- c(
  "UVP5 CCELTER 2016", # see profile "p1604_22", only upcast, not sure if downcast. Depths match with zoo data. 
  "UVP5 CCELTER (Iter2008 ccelter_2012 ccelter_2014)", # see profile "cycle_04_cast03", lower resolution below 172 m -> ignore data below 172 m
  "UVP5 Geomar 2012 msm22", # see profile "c_msm22_040", data gap between 1853 and 2094 m -> keep only data below 1853m
  "UVP5 LOHAFEX 2009", # see profile "lohafex_ps73_204", downcast data is shuffled -> keep only downcast and them reorder by depth
                       # see profile "lohafex_ps73_162", gap between 591 m and 65 1m -> keep only data below 591m
                       # see profile "lohafex_ps73_195", gap between 519 m and 579 m -> keep only data below 519m
  "UVP5 Tara Oceans (tara2009, tara2010, tara2011, tara2012, tara2013)" # see profiles "tara_037_00_b" "tara_124_00_e" "tara_146_00_a" "tara_163_03" "tara_163_04" "tara_175_10" "tara_r13" 
                                                                        # which start deeper than 10 meters while zooplankton data is available shallower
)

# Fix problems relative to specific profiles
ctd_fix <- ctd_fix %>% 
  filter(!(profile == "c_msm22_040" & depth > 1853)) %>% # fix for profile "c_msm22_040" in "UVP5 Geomar 2012 msm22": ignore data below 1853 m
  filter(!(profile == "cycle_04_cast03" & depth > 172)) %>% # fix for profile "cycle_04_cast03" in "UVP5 CCELTER (Iter2008 ccelter_2012 ccelter_2014)": ignore data below 172 m
  filter(!(profile == "lohafex_ps73_162" & depth > 591)) %>% # fix for profile "lohafex_ps73_162" in "UVP5 LOHAFEX 2009": ignore data below 591 m
  filter(!(profile == "lohafex_ps73_195" & depth > 519)) %>% # fix for profile "lohafex_ps73_195" in "UVP5 LOHAFEX 2009": ignore data below 519 m
  filter(profile != "p1604_22") %>%  # ignore profile "p1604_22" in "UVP5 CCELTER 2016" for which we don't know if we have downcast or upcast data
  filter(!(profile %in% c(
    "tara_037_00_b", 
    "tara_124_00_e", 
    "tara_146_00_a", 
    "tara_163_03", 
    "tara_163_04", 
    "tara_175_10", 
    "tara_r13")
    )) # ignore 7 profiles from "UVP5 Tara Oceans (tara2009, tara2010, tara2011, tara2012, tara2013)" starting too deep

# Check depth distribution after fix
ctd_fix %>% 
  group_by(title, profile) %>% 
  mutate(row = row_number()) %>% 
  ungroup() %>% 
  ggplot() +
  geom_path(aes(x = row, y = -depth, group = profile, color = profile), show.legend = FALSE) +
  facet_wrap(~title, scales = "free") +
  ggsave("plots/env/03.depth_dist_fix.png")

message("Depth distribution now seems ok!")
rm(ctd)


## Add fluo data from GREENEDGE 2016 ----
#--------------------------------------------------------------------------#
# Read raw fluo data
fluo_ge <- read_csv("data/raw/raw_greenedge_fluo.csv") %>% select(profile = station, depth = Pres, chla = Fluo) # ignore lon and lat because we already have them

# Extract ctd data for "UVP5hd GreenEdge 2016"
ctd_ge <- ctd_fix %>% filter(title %in% c("UVP5hd GreenEdge 2016")) %>% select(-chla)
ctd_others <- ctd_fix %>% filter(!title %in% c("UVP5hd GreenEdge 2016"))

# Join GREENEDGE fluo data on ctd data
ctd_ge <- left_join(ctd_ge, fluo_ge, by = c("profile", "depth"))

# Put back all profiles together
ctd_fix <- bind_rows(ctd_ge, ctd_others)
rm(ctd_ge, ctd_others)
message("Added fluo data for UVP5hd GreenEdge 2016")


## Clean CTD data values (abnormal values) ----
#--------------------------------------------------------------------------#
# Some anomalies to correct:
# - set temp < -9 to NA
# - set salinity < 0 or > 50 to NA
# - set chla > 100 and values of -9.99 and -999 (codes for missing values) to NA
# - set chla < 0 to 0
# - set oxy_mass < 0 to 0
# - set oxy_vol < 0 to 0
# - look at chla mean and variance and set profiles with very low variance (constant value) or high mean (>10) to NA

ctd_clean <- ctd_fix %>%
  # set aberrant values to NA
  mutate(temp = ifelse(temp < -9, NA, temp), # remove temp < 9
         sal = ifelse(sal <= 0 | sal > 50, NA, sal), # remove sal < 0 or > 50
         chla = ifelse(chla == -9.99 | chla == -999 | chla > 100, NA, chla), # remove chla > 100 and values of -9.99 and -999 (codes for missing values)
         chla = ifelse(chla <0, 0, chla), # set chla < 0 to 0
         oxy_mass = ifelse(oxy_mass < 0, 0, oxy_mass), # set oxy_mass < 0 to 0
         oxy_vol = ifelse(oxy_vol < 0, 0, oxy_vol)) %>% # set oxy_vol < 0 to 0
  # look at chla mean and variance per profile and set aberrant profiles to NA
  group_by(title, profile, lat, lon, datetime, psampleid) %>% 
  mutate(
    var_chla = var(chla, na.rm = T), # compute chla variance for each profile
    mean_chla = mean(chla, na.rm = T), # compute chla mean for each profile
    chla = ifelse(var_chla < 1e-6 | mean_chla > 10, NA, chla) # remove chla if variance is < 1e-6 (i.e. constant value of chla) or if mean is > 10 (too high values)
    ) %>% 
  ungroup() %>% 
  select(-var_chla, -mean_chla)
rm(ctd_fix)

## Look for still problematic profiles for chla
ctd_clean %>% 
  ggplot() +
  geom_path(aes(x = chla, y = -depth, group = profile, colour = profile), show.legend = F) +
  facet_wrap(~title, scales = "free_x")

# chla is missing for all profiles in two projects:
# - "UVP5 Dyfamed 2013 2014 2015 2016 2017 2018"
# - "UVP5 Geomar 2014 m108"
# Ignore all CTD data for these projects

# Abnormal values for 
# - profile "p16n_006" in "uvp5_sn009_2015_p16n": large spike around 1200 m with sensor drift -> ignore chla below 1200 m
# - profile "p16n_024" in "uvp5_sn009_2015_p16n": large spikes below 850 m with sensor drift -> ignore chla below 850 m
# - profile "gn2015_l2_002" in "UVP5 GREEN EDGE Ice Camp 2015": large spikes + max depth is 7 m -> ignore profile
# - profile "moose2013_ge_leg6_011" in "UVP5 MooseGE 2012 2013 2014 2015 2016": abnormal shape -> ignore profile
# - profile "moose2013_ge_leg6_017" in "UVP5 MooseGE 2012 2013 2014 2015 2016": abnormal values below 650 m -> ignore chla below 650 m
# - profile "moose2015_ge_leg1_003" in "UVP5 MooseGE 2012 2013 2014 2015 2016": abnormal values below 300 m -> ignore chla below 300 m
# - profile "outpace_001" in "UVP5 OUTPACE 2015": one outlier -> ignore chla values above 0.05
# - profile "outpace_002" in "UVP5 OUTPACE 2015": abnormal shape -> ignore profile
# - profile "dewex_leg1_20" in "UVP5 DEWEX 2013 (winter)": abnormal shape -> ignore profile

# Clean this!
ctd_clean <- ctd_clean %>% 
  filter(title != "UVP5 Dyfamed 2013 2014 2015 2016 2017 2018") %>% # ignore all CTD for "UVP5 Dyfamed 2013 2014 2015 2016 2017 2018"
  filter(title != "UVP5 Geomar 2014 m108") %>% # ignore all CTD for "UVP5 Geomar 2014 m108"
  # ignore problematic profiles
  mutate(
    chla = ifelse(profile == "gn2015_l2_002", NA, chla), # ignore chla in profile "gn2015_l2_002" in "UVP5 GREEN EDGE Ice Camp 2015"
    chla = ifelse(profile == "moose2013_ge_leg6_011", NA, chla), # ignore chla in profile "moose2013_ge_leg6_011" in "UVP5 MooseGE 2012 2013 2014 2015 2016"
    chla = ifelse(profile == "outpace_002", NA, chla), # ignore chla in profile "outpace_002" in "UVP5 OUTPACE 2015"
    chla = ifelse(profile == "dewex_leg1_20", NA, chla), # ignore chla in profile "dewex_leg1_20" in "UVP5 DEWEX 2013 (winter)"
  ) %>% 
  # ignore specific points within problematic profiles
  mutate(
    chla = ifelse(profile == "p16n_006" & depth > 1200, NA, chla), # ignore chla below 1200 m in profile "p16n_006" in "uvp5_sn009_2015_p16n"
    chla = ifelse(profile == "p16n_024" & depth > 850, NA, chla), # ignore chla below 850 m in profile "p16n_024" in "uvp5_sn009_2015_p16n"
    chla = ifelse(profile == "moose2013_ge_leg6_017" & depth > 650, NA, chla), # ignore chla below 650 m in profile "moose2013_ge_leg6_017" in "UVP5 MooseGE 2012 2013 2014 2015 2016"
    chla = ifelse(profile == "moose2015_ge_leg1_003" & depth > 300, NA, chla), # ignore chla below 300 m in profile "moose2015_ge_leg1_003" in "UVP5 MooseGE 2012 2013 2014 2015 2016"
    chla = ifelse(profile == "outpace_001" & chla > 0.05, NA, chla), # ignore chla above 0.05 in profile "outpace_001" in "UVP5 OUTPACE 2015"
  )

ctd_clean %>% 
  ggplot() +
  geom_path(aes(x = chla, y = -depth, group = profile), show.legend = F) +
  facet_wrap(~title, scales = "free_x")
# Some abnormalities are left, they should disappear with despiking and smoothing


## Bin at 1 meter to avoid multiple values per depth bin ----
#--------------------------------------------------------------------------#
ctd_bin <- ctd_clean %>% 
  mutate(depth = roundp(depth, precision = 1)) %>% 
  group_by(title, profile, lat, lon, datetime, psampleid, depth) %>% 
  summarise_all(list(mean), na.rm = T) %>% 
  ungroup()
rm(ctd_clean)


## Ignore profiles where one variable is entirely missing ----
#--------------------------------------------------------------------------#
ctd_profiles <- ctd_bin %>% 
  group_by(title, profile, lat, lon, datetime, psampleid) %>% # group by profile
  summarise_all(list(na_prop)) %>% # compute proportion of missing data for each variable
  ungroup() %>% 
  mutate(ignore = (temp == 1 | sal == 1 | chla == 1 | (oxy_mass == 1 & oxy_vol == 1))) %>% # ignore profile if one of the variable is entirely missing
  filter(!ignore) # keep only non ignored profile

# Ignore these profiles in ctd data
ctd_bin <- ctd_bin %>% filter(psampleid %in% ctd_profiles$psampleid)


## Despike CTD data to remove outliers ----
#--------------------------------------------------------------------------#
# Parameters for despiking:
# - k = 7 --> window = 15 (we have 1 point per m)
# - mult = 3.92
# - n.max = 5

cluster <- new_cluster(36) # use 36 cores
cluster_library(cluster, "castr") # load "castr" functions onto clusters

ctd_despike <- ctd_bin %>%
  group_by(title, profile, lat, lon, datetime, psampleid) %>%
  partition(cluster) %>% 
  mutate(
    temp     = despike(temp,     k = 7, mult = 3.92 , n.max = 5),
    sal      = despike(sal,      k = 7, mult = 3.92 , n.max = 5),
    oxy_mass = despike(oxy_mass, k = 7, mult = 3.92 , n.max = 5),
    oxy_vol  = despike(oxy_vol,  k = 7, mult = 3.92 , n.max = 5),
    chla     = despike(chla,     k = 7, mult = 3.92 , n.max = 5)
    ) %>% 
  collect() %>% 
  ungroup() %>% 
  replace(., is.na(.), NA)
rm(ctd_bin)


## Compute new variables: density and AOU ----
#--------------------------------------------------------------------------#
# Compute density from temperature and salinity
# Convert oxy_mass (µmol/kg) to oxy_vol (mL/L) (see https://ocean.ices.dk/tools/unitconversion.aspx)
# Compute AOU (apparent oxygen utilization) from oxy_umol

ctd_dens <- ctd_despike %>% 
  mutate(
    press = gsw_p_from_z(-depth, lat), # compute pressure
    SA = gsw_SA_from_SP(sal, press, lon, lat), # compute absolute salinity
    CT = gsw_CT_from_t(SA, temp, press), # compute conservative temperature
    sigma = gsw_sigma0(SA, CT), # compute sigma
    oxy_ml = ifelse(is.na(oxy_vol), oxy_mass * (sigma + 1000)/1000 * 0.022391, oxy_vol), # if oxy_vol is missing, compute oxy_ml from oxy_umol
    oxy_ml = ifelse(oxy_ml < 0, 0, oxy_ml), #  set oxy_ml < 0 to 0
    oxy_umol = oxy_ml*(1/0.022391)*1000/(sigma + 1000), # compute oxy_umol (to compute AOU) from oxy_ml 
    aou = aou(sal, temp, oxy_umol) # compute AOU
    ) %>% 
  select(title, profile, lat, lon, datetime, psampleid, depth, temp, sal, sigma, chla, oxy_ml, aou) 
rm(ctd_despike)


## Despike newly computed variables ----
#--------------------------------------------------------------------------#
ctd_dens <- ctd_dens %>%
  group_by(title, profile, lat, lon, datetime, psampleid) %>%
  partition(cluster) %>% 
  mutate(
    oxy_ml = despike(oxy_ml, k = 7, mult = 3.92 , n.max = 5),
    aou    = despike(aou,    k = 7, mult = 3.92 , n.max = 5)
    ) %>% 
  collect() %>% 
  ungroup() %>% 
  replace(., is.na(.), NA) 


## Smooth variables ----
#--------------------------------------------------------------------------#
# Parameters for smooting:
# - k = 3 --> window = 7 (we have 1 point per m)
# - n = 5 --> smooth data 5 times

ctd_smooth <- ctd_dens %>%
  group_by(title, profile, lat, lon, datetime, psampleid) %>%
  partition(cluster) %>% 
  mutate(
    temp   = smooth(temp,   k = 3, n = 5),
    sal    = smooth(sal,    k = 3, n = 5),
    sigma  = smooth(sigma,  k = 3, n = 5),
    chla   = smooth(chla,   k = 3, n = 5),
    oxy_ml = smooth(oxy_ml, k = 3, n = 5),
    aou    = smooth(aou,    k = 3, n = 5)
    ) %>% 
  collect() %>% 
  ungroup() %>% 
  replace(., is.na(.), NA)
rm(ctd_dens)


## Interpolate data over 1 meter bins ----
#--------------------------------------------------------------------------#
# Interpolation is done only if proportion of missing data for given variable and profile is lower than 20%

# Create all 1 meter bins
ctd_1m <- ctd_smooth %>% 
  select(title, profile, lat, lon, datetime, psampleid, depth) %>% 
  group_by(title, profile, lat, lon, datetime, psampleid) %>% 
  complete(depth = seq(1, max(depth), by = 1)) %>% 
  ungroup() 


ctd_use <- ctd_smooth %>% 
  group_by(title, profile, psampleid) %>% 
  summarise_at(vars(temp, sal, sigma, chla, oxy_ml, aou), list(na_prop)) %>%
  mutate(use = (temp < 0.2) & (sal < 0.2) & (sigma < 0.2) & (chla < 0.2) & (oxy_ml < 0.2) & (aou < 0.2))
sum(!ctd_use$use)
# 11 profiles have too much missing data in one of their variable
# remove these profiles from profiles to interpolate
ctd_use <- ctd_use %>% filter(use)

# List variables to interpolate
variables <- c("temp", "sal", "sigma", "chla", "oxy_ml", "aou")


# Run parallel interpolation
intl <- mclapply(unique(ctd_use$psampleid), function(id) {
  
  # Initiate empty tibble for this profile
  profile_int <- tibble()
  
  # Bins to interpolate
  ci <- ctd_1m %>% select(title, profile, lat, lon, datetime, psampleid, depth) %>% filter(psampleid == id)
  
  # Loop over variables to interpolate
  for (my_var in variables){
    # Compute interpolation
    # Available data
    wi <- ctd_smooth %>% filter(psampleid == id) %>% select(depth, all_of(my_var))
    
    # Run interpolation
    cint <- ci %>% 
      mutate(
        value = interpolate(x=wi$depth, y=pull(wi[my_var]) , xo=ci$depth, extrapolate = TRUE),
        variable = my_var
        ) %>% 
      spread(variable, value)
    
    # Join to table with other variables
    if (length(profile_int) == 0) { # If first interpolated variable on this profile
      profile_int <- cint # Replace transect table by newly computed interpolation
    } else { # Else perform a left join with previously interpolated variables
      profile_int <- left_join(profile_int, cint, by = c("title", "profile", "psampleid", "depth", "lat", "lon", "datetime"))
    } 
  }
  
  return(profile_int)
}, mc.cores=24) # on 24 cores
# this returns a list, recombine it into a tibble
ctd_int <- do.call(bind_rows, intl)
rm(ctd_smooth, ctd_1m)


## Bin CTD data at 5 m to match with plankton data ----
#--------------------------------------------------------------------------#
# Run in parallel
ctd_5m <- ctd_int %>%
  mutate(depth = roundp(depth + 1.2, 5, f=round) + 2.5) %>% # add 1.2 m to depth to account for distance between CTD and UVP and compute closest 5 m bin
  #mutate(depth = roundp(depth, 5, f=round) + 2.5) %>%
  group_by(title, profile, lat, lon, datetime, psampleid, depth) %>%
  partition(cluster) %>% 
  summarise_all(list(mean), na.rm = T) %>%
  collect() %>% 
  ungroup() %>% 
  replace(., is.na(.), NA) %>% 
  arrange(title, profile, depth)


## Compute clines ----
#--------------------------------------------------------------------------#
# Compute thermocline, halocline, mixed layer depth, pycnocline, deep chlorophyll maximum, euphotic depth and stratification index
# Compute epipelagic layer depth as the maximum depth between euphotic zone and pycnocline or MLD

layers <- ctd_int %>% 
  group_by(title, profile, lat, lon, datetime, psampleid) %>% # group by profile
  partition(cluster) %>% 
  summarise( # compute:
    max_depth = max(depth),                                                # maximum depth of profile
    thermo = clined(temp, depth, n.smooth=2, k=2),                         # thermocline from temperature
    halo = clined(sal, depth, n.smooth=2, k=2),                            # halocline from salinity
    mld = mld(sigma, depth),                                               # mixed layer depth from density
    pyc = depth[which.max(slide(sigma, k = 2, sd))],                       # pycnocline from density
    dcm = maxd(chla),                                                      # deep chlorophyll maximum from chla
    ze = ze_from_chla(chla, depth),                                        # euphotic depth from chla
    strat = stratif(sigma, depth, min.depths = 1:5, max.depths = 195:200)  # stratification index from density between surface (1:5 m) and (195:200 m)
    ) %>% 
  collect() %>% 
  ungroup() %>% 
  mutate( # compute epipelagic layer depth
    epi_pyc = ifelse(ze > pyc, ze, pyc), # from ze and pycnocline
    epi_mld = ifelse(ze > mld, ze, mld), # from ze and mld
    ) 

# Inspect depth of epipelagic layer
layers %>% 
  select(psampleid, epi_pyc, epi_mld) %>% 
  gather(epi_pyc:epi_mld, key = "method", value = "depth") %>% 
  ggplot() +
  geom_histogram(aes(x = depth, fill = method), alpha = 0.5, binwidth = 5, position = "identity") +
  ggtitle("Epipelagic layer depth by two computation methods")
ggsave("plots/env/03.epi_histo.png")

# Plot sigma profile, pyc and MLD for these profiles
# When using pycnocline, 8 profiles have an epipelagic layer deeper than 500 m. 
ctd_int %>% 
  filter(psampleid %in% (layers %>% filter(pyc > 500) %>% pull(psampleid))) %>% 
  left_join(layers) %>% 
  select(profile, psampleid, depth, sigma, pyc, mld) %>% 
  gather(pyc:mld, key = "cline", value = "cline_depth") %>% 
  ggplot() +
  geom_point(aes(x = sigma, y = -depth, group = profile), size = 0.1) +
  geom_hline(aes(yintercept = -cline_depth, color = cline)) +
  facet_wrap(~profile) +
  ggtitle("Sigma profiles of deep pycnocline cases") +
  ggsave("plots/env/03.sigma_profiles_deep_pyc.png")

# This does not happen when using MLD, while the overall distribution is similar --> use epi computed with MLD
layers <- layers %>% 
  rename(epi = epi_mld) %>% 
  select(-epi_pyc)
layers %>% filter(is.na(epi)) %>% summary()
# Epipelagic layer depth could not be computed for 84 profiles which are too shallow (max depth < 150 m) --> ignore these profiles
layers <- layers %>% filter(!is.na(epi))

# Inspect chla profiles where Ze is 180 m (63 profiles)
ctd_int %>% 
  filter(psampleid %in% (layers %>% filter(ze == 180) %>% pull(psampleid))) %>% 
  left_join(layers) %>% 
  select(profile, psampleid, depth, chla, ze, mld) %>% 
  gather(ze:mld, key = "cline", value = "cline_depth") %>% 
  ggplot() +
  geom_point(aes(x = chla, y = -depth, group = profile), size = 0.1) +
  geom_hline(aes(yintercept = -cline_depth, color = cline)) +
  ylim(-300,0) +
  facet_wrap(~profile, scales = "free") +
  ggtitle("Chla profiles of cases where Ze = 180 m")
ggsave("plots/env/03.chla_profiles_deep_ze.png")
# There does not seem to be something wrong with these profiles.
# They probably just have low chla values, leading to deep ze values.


## Remove profiles with incomplete data  ----
#--------------------------------------------------------------------------#
# Keep only profiles for which we were able to compute epipelagic layer depth
ctd_5m <- ctd_5m %>% filter(psampleid %in% layers$psampleid)
ctd_int <- ctd_int %>% filter(psampleid %in% layers$psampleid)


## Save data ----
#--------------------------------------------------------------------------#
# Save CTD data on 5 meter bins
save(ctd_5m, file = "data/03.ctd_5m.Rdata")

# Save interpolated CTD data on 1 meter bins
save(ctd_int, file = "data/03.ctd_int.Rdata")

# Save layers data
save(layers, file = "data/03.layers.Rdata")

