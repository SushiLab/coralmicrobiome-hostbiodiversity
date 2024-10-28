## Project: Aquarium Microbiome
## Script purpose: Zoom in based on taxonomy and examine pathway compositions
## Date: August 2024
## Author: Fabienne Wiederkehr

# Load packages
require(data.table)
require(tidyverse)
require(ComplexUpset)

# Load KEGG info
path_list_cat<-fread("../data/kegg/pathway_list_cat.csv",header=TRUE,data.table=FALSE) %>%
  filter(grepl("Metabolism",pathway_category)) %>% #keep only pathways associated with metabolism
  filter(!pathway_class=="Global and overview maps") #remove global and overview maps
ko_pathway_link<-fread("../data/kegg/ko_pathway_link.tsv",header=TRUE,data.table=FALSE) %>%
  filter(pathway_id %in% path_list_cat$pathway_id)

# Load data
ko_abtab<-fread("../data/ko_abtab_asv.tsv",sep="\t",header=TRUE,data.table=FALSE)

# Select ASVs
source<-c("asv_000202","asv_000307")

# Summarise KOs to pathways
path_abtab<-ko_abtab %>%
  left_join(ko_pathway_link,by="ko_id") %>%
  group_by(pathway_id) %>%
  summarise(across(where(is.numeric),~sum(.x,na.rm=TRUE))) %>%
  na.omit() %>%
  column_to_rownames("pathway_id") %>%
  select(all_of(source))

# Remove all-0 pathways
path_abtab<-path_abtab[which(apply(path_abtab,1,sum)>0),] #remove all-0

# Annotate pathways and convert to presence-absence
path_abtab_id<-path_abtab %>%
  rownames_to_column("pathway_id") %>%
  left_join(path_list_cat,by="pathway_id") %>%
  mutate_if(is.numeric,as.logical)

# Generate upset plot
upset(path_abtab_id,source,
      width_ratio=0.3,min_degree=1,
      base_annotations=list(
        "Number of pathways"=intersection_size(counts=TRUE)),
      set_sizes=FALSE,
      stripes="transparent",
      name=NULL,
      sort_sets=FALSE,
      guides="over")
ggsave(filename="../data/40_FunPred_upset_Sinularia.png",height=unit(6,"cm"),width=unit(12,"cm"))
