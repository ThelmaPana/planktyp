#
# List EcoPart projects and determine which to use
#
# (c) 2020 Jean-Olivier Irisson, GNU General Public License v3

source("lib/extract_zoo/0.setup.R")

# get ids and titles from the database
projects <- tbl(db, "part_projects") %>%
  select(pprojid, ptitle, data_owner=do_email, projid) %>%
  # add EcoTaxa title
  left_join(tbl(db, "projects") %>% select(projid, title)) %>%
  collect() %>%
  arrange(pprojid) %>%
  # add empty column
  mutate(use="")

# write it into a file
#write_tsv(projects, "data/projects/UVP5_projects.tsv")
write_tsv(projects, file.path(data_dir, "UVP5_projects.tsv"))
# now select which ones can be used in the global dataset

