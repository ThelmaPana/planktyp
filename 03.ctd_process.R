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
library(cmocean)

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
  facet_wrap(~title, scales = "free")
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


## Convert oxy_mass (µmol/kg) to oxy_vol (mL/L) ----
#--------------------------------------------------------------------------#
# (see https://www.ices.dk/data/tools/Pages/Unit-conversions.aspx)

# For UVP5hd GreenEdge 2016 and UVP5 MALINA 2009, values in oxy_mass (µmol/kg) are actually oxy_vol (mL/L)
ctd_oxy <- ctd_bin %>% 
  # Correct oxy_mass to oxy_vol for GreenEdge 2016 and UVP5 MALINA 2009
  mutate(
    oxy_vol = ifelse(title %in% c("UVP5hd GreenEdge 2016", "UVP5 MALINA 2009"), oxy_mass, oxy_vol),
    oxy_mass = ifelse(title %in% c("UVP5hd GreenEdge 2016", "UVP5 MALINA 2009"), NA, oxy_mass)
  ) %>% 
  mutate(
    press = gsw_p_from_z(-depth, lat), # compute pressure
    SA = gsw_SA_from_SP(sal, press, lon, lat), # compute absolute salinity
    CT = gsw_CT_from_t(SA, temp, press), # compute conservative temperature
    sigma = gsw_sigma0(SA, CT), # compute sigma
    oxy_ml = ifelse(is.na(oxy_vol), oxy_mass * (sigma + 1000)/1000 * 0.022391, oxy_vol), # if oxy_vol is missing, compute oxy_ml from oxy_umol
    oxy_ml = ifelse(oxy_ml < 0, 0, oxy_ml), #  set oxy_ml < 0 to 0
    #oxy_umol = oxy_ml*(1/0.022391)*1000/(sigma + 1000), # compute oxy_umol (to compute AOU) from oxy_ml 
  ) %>% 
  select(title, profile, lat, lon, datetime, psampleid, depth, temp, sal, chla, oxy_ml) 
rm(ctd_bin)


## Interpolate data over 1 meter bins ----
#--------------------------------------------------------------------------#
# Interpolation is done only if proportion of missing data for given variable and profile is lower than 20%

# Create all 1 meter bins
ctd_1m <- ctd_oxy %>% 
  select(title, profile, lat, lon, datetime, psampleid, depth) %>% 
  group_by(title, profile, lat, lon, datetime, psampleid) %>% 
  complete(depth = seq(1, max(depth), by = 1)) %>% 
  ungroup() 

ctd_use <- ctd_oxy %>% 
  group_by(title, profile, psampleid) %>% 
  summarise_at(vars(temp, sal, chla, oxy_ml), list(na_prop)) %>%
  mutate(use = (temp < 0.2) & (sal < 0.2) & (chla < 0.2) & (oxy_ml < 0.2))
sum(!ctd_use$use)
# 18 profiles have too much missing data in one of their variable
# remove these profiles from profiles to interpolate
ctd_use <- ctd_use %>% filter(use)

# List variables to interpolate
variables <- c("temp", "sal", "chla", "oxy_ml")


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
    wi <- ctd_oxy %>% filter(psampleid == id) %>% select(depth, all_of(my_var))
    
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
rm(ctd_oxy, ctd_1m)


## Despike CTD data to remove outliers ----
#--------------------------------------------------------------------------#
# Parameters for despiking:
# - k = 7 --> window = 15 (we have 1 point per m)
# - mult = 3.92
# - n.max = 5

cluster <- new_cluster(36) # use 36 cores
cluster_library(cluster, "castr") # load "castr" functions onto clusters

ctd_despike <- ctd_int %>%
  group_by(title, profile, lat, lon, datetime, psampleid) %>%
  partition(cluster) %>% 
  mutate(
    temp     = despike(temp,   k = 7, mult = 3.92 , n.max = 5),
    sal      = despike(sal,    k = 7, mult = 3.92 , n.max = 5),
    chla     = despike(chla,   k = 7, mult = 3.92 , n.max = 5),
    oxy_ml   = despike(oxy_ml, k = 7, mult = 3.92 , n.max = 5)
    ) %>% 
  collect() %>% 
  ungroup() %>% 
  replace(., is.na(.), NA)
rm(ctd_int)


## Smooth variables ----
#--------------------------------------------------------------------------#
# Parameters for smoothing:
# - k = 3 --> window = 7 (we have 1 point per m)
# - n = 5 --> smooth data 5 times

ctd_smooth <- ctd_despike %>%
  group_by(title, profile, lat, lon, datetime, psampleid) %>%
  partition(cluster) %>% 
  mutate(
    temp   = smooth(temp,   k = 3, n = 5),
    sal    = smooth(sal,    k = 3, n = 5),
    chla   = smooth(chla,   k = 3, n = 5),
    oxy_ml = smooth(oxy_ml, k = 3, n = 5)
    ) %>% 
  collect() %>% 
  ungroup() %>% 
  replace(., is.na(.), NA)
rm(ctd_despike)


## Compute new variables: density and AOU ----
#--------------------------------------------------------------------------#
# Compute density from temperature and salinity
# Compute AOU (apparent oxygen utilization) from oxy_umol





#gsw_o2sol_sp_pt <- function(sp, pt){
#  # Calculates the oxygen concentration expected at equilibrium with air at
#  # an Absolute Pressure of 101325 Pa (sea pressure of 0 dbar) including 
#  # saturated water vapor.  This function uses the solubility coefficients 
#  # derived from the data of Benson and Krause (1984), as fitted by Garcia 
#  # and Gordon (1992, 1993).
#  #
#  # Note that this algorithm has not been approved by IOC and is not work 
#  # from SCOR/IAPSO Working Group 127. It is included in the GSW
#  # Oceanographic Toolbox as it seems to be oceanographic best practice.
#  #
#  # SP  :  Practical Salinity  (PSS-78)                         [ unitless ]
#  # pt  :  potential temperature (ITS-90) referenced               [ dbar ]
#  #         to one standard atmosphere (0 dbar).
#  #
#  # gsw_o2sol_sp_pt : solubility of oxygen in micro-moles per kg     [umol/kg]
#  
#  x = sp;
#  
#  pt68 = pt*1.00024;
#  
#  y = log((298.15 - pt68)/(273.15 + pt68));
#  
#  a0 =  5.80871; 
#  a1 =  3.20291;
#  a2 =  4.17887;
#  a3 =  5.10006;
#  a4 = -9.86643e-2;
#  a5 =  3.80369;
#  b0 = -7.01577e-3;
#  b1 = -7.70028e-3;
#  b2 = -1.13864e-2;
#  b3 = -9.51519e-3;
#  c0 = -2.75915e-7;
#  
#  o2sol = exp(a0 + y*(a1 + y*(a2 + y*(a3 + y*(a4 + a5*y))))
#              + x*(b0 + y*(b1 + y*(b2 + b3*y)) + c0*x));
#  
#  return (o2sol)
#  
#}
#
#gsw_o2sol_sp_pt(sp = toto$sal, pt = toto$temp)
#toto$oxy_ml
#
#
#library(tidync)
#file <- "data/raw/woa/woa18_all_A00_01.nc"
#aou <- tidync(file) %>% 
#  hyper_filter(depth = depth <= 500) %>% 
#  hyper_tibble() %>% 
#  select(lat, lon, depth, aou_woa = A_an) %>% 
#  unique()
#
#O2sol = gsw_O2sol_SP_pt(SP,pt)

ctd_dens <- ctd_smooth %>% 
  mutate(
    press = gsw_p_from_z(-depth, lat), # compute pressure
    SA = gsw_SA_from_SP(sal, press, lon, lat), # compute absolute salinity
    CT = gsw_CT_from_t(SA, temp, press), # compute conservative temperature
    sigma = gsw_sigma0(SA, CT), # compute sigma
    oxy_umol = oxy_ml*(1/0.022391)*1000/(sigma + 1000), # compute oxy_umol (to compute AOU) from oxy_ml 
    aou = aou(sal, temp, oxy_umol) # compute AOU
  ) %>% 
  select(title, profile, lat, lon, datetime, psampleid, depth, temp, sal, sigma, chla, oxy_ml, aou) 
rm(ctd_smooth)


## Bin CTD data at 5 m to match with plankton data ----
#--------------------------------------------------------------------------#
# Run in parallel
ctd_5m <- ctd_dens %>%
  #mutate(depth = roundp(depth + 1.2, 5, f=round) + 2.5) %>% # add 1.2 m to depth to account for distance between CTD and UVP and compute closest 5 m bin
  mutate(depth = roundp(depth, 5, f=round) + 2.5) %>% # compute closest 5 m bin
  group_by(title, profile, lat, lon, datetime, psampleid, depth) %>%
  partition(cluster) %>% 
  summarise_all(list(mean), na.rm = T) %>%
  collect() %>% 
  ungroup() %>% 
  replace(., is.na(.), NA) %>% 
  arrange(title, profile, depth)

summary(ctd_5m)


## Compute clines ----
#--------------------------------------------------------------------------#
# Compute thermocline, halocline, mixed layer depth, pycnocline, deep chlorophyll maximum, euphotic depth and stratification index
# Compute epipelagic layer depth as the maximum depth between euphotic zone and pycnocline or MLD

layers <- ctd_dens %>% 
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

# Inspect depth of pycnocline
layers %>% 
  ggplot() +
  geom_histogram(aes(x = pyc), binwidth = 5)
sum(layers$pyc > 250)
# 32 profiles with pycnocline deeper than 250 m

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
ctd_dens %>% 
  filter(psampleid %in% (layers %>% filter(pyc > 500) %>% pull(psampleid))) %>% 
  left_join(layers) %>% 
  select(profile, psampleid, depth, sigma, pyc, mld) %>% 
  gather(pyc:mld, key = "cline", value = "cline_depth") %>% 
  ggplot() +
  geom_point(aes(x = sigma, y = -depth, group = profile), size = 0.1) +
  geom_hline(aes(yintercept = -cline_depth, color = cline)) +
  facet_wrap(~profile) +
  ggtitle("Sigma profiles of deep pycnocline cases")
ggsave("plots/env/03.sigma_profiles_deep_pyc.png")

# This does not happen when using MLD, while the overall distribution is similar --> use epi computed with MLD
layers <- layers %>% 
  rename(epi = epi_mld) %>% 
  select(-epi_pyc)
layers %>% filter(is.na(epi)) %>% summary()
# Epipelagic layer depth could not be computed for 84 profiles which are too shallow (max depth < 150 m) --> ignore these profiles
layers <- layers %>% filter(!is.na(epi))

# Inspect chla profiles where Ze is 180 m (63 profiles)
ctd_dens %>% 
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
ctd_1m <- ctd_dens %>% filter(psampleid %in% layers$psampleid)







## Plot maps of surface data ----
#--------------------------------------------------------------------------#

ctd_5m %>% summary()
ctd_5m %>% 
  ggplot() +
  geom_density(aes(x = aou))

# Temp
ctd_5m %>% 
  filter(depth == 2.5) %>% 
  ggplot() +
  geom_polygon(aes(x = lon, y = lat, group = group), fill = "gray", data = world) +
  geom_point(aes(x = lon, y = lat, color = temp)) +
  scale_color_cmocean(name = "thermal") +
  theme_minimal() +
  coord_quickmap()

# Sal
ctd_5m %>% 
  filter(depth == 2.5) %>% 
  ggplot() +
  geom_polygon(aes(x = lon, y = lat, group = group), fill = "gray", data = world) +
  geom_point(aes(x = lon, y = lat, color = sal)) +
  scale_color_cmocean(name = "haline") +
  theme_minimal() +
  coord_quickmap()

# Sal
ctd_5m %>% 
  filter(depth == 2.5) %>% 
  ggplot() +
  geom_polygon(aes(x = lon, y = lat, group = group), fill = "gray", data = world) +
  geom_point(aes(x = lon, y = lat, color = sigma)) +
  scale_color_cmocean(name = "dense") +
  theme_minimal() +
  coord_quickmap()

# Oxy
ctd_5m %>% 
  filter(depth == 2.5) %>% 
  ggplot() +
  geom_polygon(aes(x = lon, y = lat, group = group), fill = "gray", data = world) +
  geom_point(aes(x = lon, y = lat, color = oxy_ml)) +
  scale_color_distiller(palette = "Blues") +
  theme_minimal() +
  coord_quickmap()

# Chla
ctd_5m %>% 
  filter(depth == 2.5) %>% 
  ggplot() +
  geom_polygon(aes(x = lon, y = lat, group = group), fill = "gray", data = world) +
  geom_point(aes(x = lon, y = lat, color = chla)) +
  scale_color_cmocean(name = "algae") +
  theme_minimal() +
  coord_quickmap()

# AOU
ctd_5m %>% 
  filter(depth == 2.5) %>% 
  ggplot() +
  geom_polygon(aes(x = lon, y = lat, group = group), fill = "gray", data = world) +
  geom_point(aes(x = lon, y = lat, color = aou)) +
  scale_color_cmocean(name = "matter") +
  theme_minimal() +
  coord_quickmap()

ctd_5m %>% 
  group_by(title) %>% 
  summarise(aou = mean(aou)) %>% 
  arrange(desc(aou))


## Save data ----
#--------------------------------------------------------------------------#
# Save CTD data on 5 meter bins
save(ctd_5m, file = "data/03.ctd_5m.Rdata")

# Save interpolated CTD data on 1 meter bins
save(ctd_1m, file = "data/03.ctd_1m.Rdata")

# Save layers data
save(layers, file = "data/03.layers.Rdata")

