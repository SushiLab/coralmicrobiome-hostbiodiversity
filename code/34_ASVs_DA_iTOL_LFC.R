## Project: Aquarium Microbiome
## Script purpose: Pull the significantly differentially abundant ASVs
## Date: March 2024
## Author: Fabienne Wiederkehr

# Load packages
require(data.table)
require(tidyverse)

# Load data
asv_da<-fread("../data/asv_diffab.tsv",sep="\t",header=TRUE,data.table=FALSE)

# Transform dataframe into wider format
asv_da_long<-asv_da %>%
  select(-q_treatment_binarydegraded) %>%
  pivot_wider(names_from=species,values_from=lfc_treatment_binarydegraded,values_fill=0) %>%
  mutate(`Haliclona cnidata`=0)

# Arrange dataframe as matrix
asv_da_mat<-asv_da_long %>%
  select(asv,contains(" ")) %>%
  column_to_rownames("asv") %>% as.matrix()

# Get inverse matrix to compare biodiverse vs degraded
asv_da_inv<-asv_da_mat * -1

# Export lfc table
asv_da_lfc<-asv_da_inv %>% as.data.frame() %>%
  rownames_to_column("asv")
fwrite(asv_da_lfc,"../data/asv_diffab_lfc.tsv",sep="\t",na="NA")
