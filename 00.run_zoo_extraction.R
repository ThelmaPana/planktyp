#--------------------------------------------------------------------------#
# Project: planktyp
# Script purpose: run JO's extraction of EcoTaxa objects (see https://github.com/jiho/UVP5_images_dataset)
# Date: 01/02/2021
# Author: Thelma Panaiotis
#--------------------------------------------------------------------------#

# JO's scripts:
# - 0.setup.R
# - 1.list_projects.R
# - 2.get_samples.R
# - 3.get_objects.R
# - 4.get_volumes.R
# - 5.reformat_objects.R
# - 6.regroup_taxa.R

source("lib/extract_zoo/0.setup.R")
source("lib/extract_zoo/1.list_projects.R")
source("lib/extract_zoo/2.get_samples.R")
source("lib/extract_zoo/3.get_objects.R")
source("lib/extract_zoo/4.get_volumes.R")
source("lib/extract_zoo/5.reformat_objects.R")
source("lib/extract_zoo/6.regroup_taxa.R")