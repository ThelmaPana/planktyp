#--------------------------------------------------------------------------#
# Project: plankty
# Script purpose: Repo set-up
# Date: 03/02/2021
# Author: Thelma Panaiotis
#--------------------------------------------------------------------------#


## Load packages ----
#--------------------------------------------------------------------------#
suppressMessages(library(tidyverse))
library(lubridate)


## Create directories for plots ----
#--------------------------------------------------------------------------#
dir.create("plots/zoo", showWarnings=FALSE, recursive=TRUE) # for zooplankton data
dir.create("plots/env", showWarnings=FALSE, recursive=TRUE) # for environmental data
dir.create("plots/provinces", showWarnings=FALSE, recursive=TRUE) # for provinces
dir.create("plots/analysis/epi", showWarnings=FALSE, recursive=TRUE) # for epipelagic analysis
dir.create("plots/analysis/mesosup", showWarnings=FALSE, recursive=TRUE) # for mesopelagic sup analysis
dir.create("plots/analysis/mesoinf", showWarnings=FALSE, recursive=TRUE) # for mesopelagic inf analysis
dir.create("plots/analysis/bathy", showWarnings=FALSE, recursive=TRUE) # for bathypelagic analysis


## World background for maps ----
#--------------------------------------------------------------------------#
world <- fortify(map_data('world')) %>% rename(lon=long)


## Custom functions ----
#--------------------------------------------------------------------------#
# Proportion of NA in a vector
na_prop <- function (x){sum(is.na(x))/length(x)}

# Sum of NA in a vector
na_sum <- function (x){sum(is.na(x))}

# Convert degrees to radians
deg2rad <- function(deg){(deg * pi) / (180)}

# Compute AOU
aou <- function(sal, temp, oxy){
  ## Compute apparent oxygen utilization
  # AOU is computed as the difference between O2 concentration at saturation and observed O2 concentration
  # Input
  #     - sal : salinity
  #     - temp : temperature (°C)
  #     - oxy : oxygen concentration (µmol/Kg)
  #
  # Output
  #     - aou : apparent oxygen utilisation (µmol/Kg)
  
  # Convert temp to kelvin
  temp1 <-  (temp + 273.15) / 100
  
  # Compute oxygen concentration (mL/Kg) at saturation 
  osat <-  -177.7888 + 255.5907 / temp1 + 146.4813 * log(temp1) - 22.2040 * temp1
  osat <- osat + sal * (-0.037362 + temp * (0.016504 - 0.0020564 * temp))
  osat <-  exp(osat)
  
  # Convert from mL/Kg to µmol/kg
  osat <-  osat * 1000 / 22.392
  
  # Compute AOU
  aou <- osat - oxy
  
  return(aou)
} 

# Standardize a vector to mean = 0 and sd = 1
scale2 <- function(x, na.rm = T){(x - mean(x, na.rm = na.rm)) / sd(x, na.rm)}

# Generate points of a circle
circle_fun <- function(center = c(0,0), diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

# Count modalities of a qualitative variable
count_mods <- function(x){length(unique(x))}
