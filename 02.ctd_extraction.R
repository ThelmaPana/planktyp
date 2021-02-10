#--------------------------------------------------------------------------#
# Project: plantyp
# Script purpose: Extract CTD data from ecoparts using psampleid
# Date: 01/02/2021
# Author: Thelma Panaiotis
#--------------------------------------------------------------------------#

source("lib/set_up.R")
library(ecotaxar)
library(parallel)

load("data/01.bio_data.Rdata")


## Extract CTD data ----
#--------------------------------------------------------------------------#
# Connect to EcoTaxa database
db <- db_connect_ecotaxa()

# List samples
smp <- zoo_conc %>% 
  select(title, profile, sampleid, psampleid, lat, lon) %>% 
  unique() %>% 
  left_join(tbl(db, "part_samples") %>% select(psampleid, datetime = sampledate, ctd_desc), copy=TRUE) %>% 
  collect()

# Parallel extraction of ctd data
ctdl <- mclapply(unique(smp$psampleid), function(id) {
  x <- smp %>% filter(psampleid == id)
  
  # Get CTD fields description
  if (is.na(x$ctd_desc[1]) | x$ctd_desc[1] == "") {
    # if no extra columns are defined, to not extract anything more
    desc <- NULL
  } else {
    # when extra columns are defined, get their specification
    desc <- parse_mapping(iconv(x$ctd_desc[1], from="latin1", to="utf8"))
    # remove special characters in the variable names
    clean_names <- names(desc) %>% 
      str_replace("µ", "u") %>%
      make.names() %>%
      str_replace_all("\\.{2,}", ".") %>% # remove double dots
      str_replace_all("\\.$" , "")        # remove dots at end f name
    # prepare the variable description for select()
    desc <- str_c("extrames", desc)
    names(desc) <- clean_names
  }
  
  # get CTD data
  ctd_id <- tbl(db, "part_ctd") %>% 
    # for the samples of interest
    #filter(psampleid %in% !!x$psampleid) %>%
    filter(psampleid == as.integer(id)) %>%
    select(
      psampleid,
      # get usual variables
      depth,
      temp=temperature, sal=practical_salinity, dens=in_situ_density_anomaly,
      oxygen_mass, oxygen_vol, nitrate, chla=chloro_fluo,
      # and add potential extra ones
      desc
    ) %>%
    collect()
  
  return(ctd_id)
}, mc.cores=1) # on 1 core because error otherwise
# this returns a list, recombine it into a tibble
ctd_raw <- do.call(bind_rows, ctdl)
rm(ctdl)

message(
  "ctd data found for ", length(unique(ctd_raw$psampleid)), " psamples",
  " out of ", length(unique(smp$psampleid)), " in zoo data"
)

## Add metadata
ctd_raw <- ctd_raw %>%
  # keep information in one of two columns for the redundant variables
  mutate(
    oxy_mass=ifelse(is.na(oxygen_mass), Oxygen.umol.Kg, oxygen_mass),
    oxy_vol=ifelse(is.na(oxygen_vol), Oxygen.ml.l, oxygen_vol)
  ) %>% 
  # keep only the relevant columns
  select(
    psampleid,
    depth,
    temp, sal,
    oxy_mass, oxy_vol, chla
  ) %>% 
  # add project and profile data
  left_join(smp %>% select(psampleid, title, profile, datetime, lat, lon)) %>% 
  relocate(title, profile, lat, lon, datetime)

# Disconnect from database
db_disconnect_ecotaxa(db)


## Save raw CTD data ----
#--------------------------------------------------------------------------#
save(ctd_raw, file = "data/02.ctd_raw.Rdata")

# write a csv.gz file for data storage 
#write.csv(ctd_raw, file=gzfile("data/export/02.ctd_raw.csv.gz"))
