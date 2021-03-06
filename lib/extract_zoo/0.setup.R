#
# Perform common actions for all scripts
#
# (c) 2020 Jean-Olivier Irisson, GNU General Public License v3

# load packages
suppressMessages(library("tidyverse"))
library("ecotaxar")
db <- db_connect_ecotaxa()

# create directories for large data that should live outside the repository
data_dir <- "data/extract_zoo"
proj_dir <- "data/extract_zoo/projects"
#img_dir <- "~/UVP5_images_dataset/images/"

dir.create(data_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(proj_dir, showWarnings=FALSE, recursive=TRUE)
#dir.create(img_dir, showWarnings=FALSE, recursive=TRUE)