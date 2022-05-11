#
# Regroup taxa into consistently sorted and mutually exclusive groups
#
# (c) 2020 Jean-Olivier Irisson, GNU General Public License v3

source("lib/extract_zoo/0.setup.R")
library("feather")

# read all objects information
o <- read_feather("data/00.all_zoo.feather")

#g <- read_csv("https://docs.google.com/spreadsheets/d/1cILxbd4BN3Qez0bSf2Up69Gh2-Jdbzp4aXGsZ31RMUg/export?format=csv") %>%
#g <- read_csv("https://docs.google.com/spreadsheets/d/1NFgpzkFXVBEuobaggmApIYgQ2HE837zwsStu7WTU5hQ/export?format=csv", col_types=cols()) %>%
#  drop_na(taxon)
# TODO replace by a static .csv file once this is stabilised
g <- read_csv("data/extract_zoo/UVP5_taxo_regrouped.csv", col_types=cols()) %>% drop_na(taxon)


# Try to infer a lineage for the new grouping
# get all groups
groups <- distinct(g, group_thelma)
# get the available taxonomy
taxo <- distinct(o, lineage, taxon)
# start with perfect matches
group_lineages <- groups %>%
  mutate(
    group_for_match=group_thelma %>%
      str_replace("_others", "") %>%
      str_replace("_cavo_or_creseis", "")
  ) %>%
  left_join(taxo, by=c("group_for_match"="taxon"))
# and fill the rest manually
group_lineages$lineage[which(group_lineages$group_for_match=="Nostocales")] <- "living/Bacteria/Proteobacteria/Cyanobacteria/Nostocales"
# filter(taxo, str_detect(lineage, "Nostocales"))
group_lineages$lineage[which(group_lineages$group_for_match=="Cnidaria")] <- "living/Eukaryota/Opisthokonta/Holozoa/Metazoa/Cnidaria"
group_lineages$lineage[which(group_lineages$group_for_match=="Collodaria_colonial")] <- "	living/Eukaryota/Harosa/Rhizaria/Retaria/Polycystinea/Collodaria/colonial"
# filter(taxo, str_detect(lineage, "Cnidaria"))
group_lineages$lineage[which(group_lineages$group_for_match=="misc?")] <- "living"
group_lineages$lineage[which(group_lineages$group_for_match=="?")] <- "living"

# add lineage to the grouping
g <- left_join(g, group_lineages)

# get a few stats
# g %>% drop_na(group) %>% group_by(lineage, group) %>% summarise(n=sum(n), .groups="drop") %>% View()

# and add to the objects file
all_zoo <- left_join(o, select(g, taxon, group=group_thelma, group_lineage=lineage) %>% na.omit())
sort(unique(all_zoo$group))

all_zoo %>% pull(group) %>% unique() %>% sort()
all_zoo %>% filter(str_detect(group, "colonial"))
all_zoo %>% filter(group == "Collodaria")
o %>% filter(str_detect(taxon, "colonial"))

write_feather(all_zoo, "data/00.all_zoo.feather")

# write a csv.gz file for data storage 
#write.csv(all_zoo, file=gzfile("data/export/00.all_zoo.csv.gz"))
