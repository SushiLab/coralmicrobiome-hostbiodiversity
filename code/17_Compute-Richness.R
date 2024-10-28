## Project: Aquarium Microbiome
## Script purpose: Compute richness
## Date: May 2024
## Author: Fabienne Wiederkehr

# Load packages
require(data.table)
require(tidyverse)
require(vegan)
require(hillR)

# Load paths
asv_abtab_paths<-list.files("../data/rarefied/",pattern="asv_abtab_rarefied",full.names=TRUE)

# Load metadata
metadat<-fread("../data/rarefied/metadat_rarefied.tsv",sep="\t",header=TRUE,data.table=FALSE)

# Compute richness
asv_richtab_list<-list(NULL)
for (i in 1:length(asv_abtab_paths)){
  # Load data
  asv_abtab<-fread(asv_abtab_paths[i],sep="\t",header=TRUE,data.table=FALSE) %>%
    column_to_rownames("id") %>%
    filter(rownames(.) %in% metadat$id)
  
  # Compute Hill 0 (richness)
  hill_zero<-hill_taxa(asv_abtab,q=0) %>% as.data.frame() %>%
    rename("hill_zero"=".")
  
  # Compute Hill 1 (exponential of Shannon’s entropy index)
  hill_one<-hill_taxa(asv_abtab,q=1) %>% as.data.frame() %>%
    rename("hill_one"=".")
  
  # Compute Hill 2 (inverse of Simpson’s concentration index)
  hill_two<-hill_taxa(asv_abtab,q=2) %>% as.data.frame() %>%
    rename("hill_two"=".")
  
  # Bind Hill dataframe
  asv_richtab_tmp<-hill_zero %>%
    bind_cols(hill_one) %>%
    bind_cols(hill_two)
  
  # Save richness in list
  asv_richtab_list[[i]]<-asv_richtab_tmp
}

# Reduce list to a single dataframe
asv_richtab_combined<-Reduce("+",asv_richtab_list) #sum up the corresponding elements
asv_richtab<-asv_richtab_combined/length(asv_richtab_list) #calculate mean of each element

# Write richness dataframe
asv_richtab<-asv_richtab %>%
  rownames_to_column("id")
fwrite(asv_richtab,"../data/asv_richtab.tsv",sep="\t",na="NA")
