## Project: Aquarium Microbiome
## Script purpose: Compute distances
## Date: January 2024
## Author: Fabienne Wiederkehr

# Load packages
require(data.table)
require(tidyverse)
require(vegan)

# Load paths
asv_abtab_paths<-list.files("../data/rarefied/",pattern="asv_abtab_rarefied",full.names=TRUE)

# Load metadata
metadat<-fread("../data/rarefied/metadat_rarefied.tsv",sep="\t",header=TRUE,data.table=FALSE)

# Compute BC distances
asv_bctab_list<-list(NULL)
for (i in 1:length(asv_abtab_paths)){
  # Load data
  asv_abtab<-fread(asv_abtab_paths[i],sep="\t",header=TRUE,data.table=FALSE) %>%
    column_to_rownames("id") %>%
    filter(rownames(.) %in% metadat$id)
  
  # Compute BC distances
  asv_bctab_tmp<-vegdist(sqrt(asv_abtab)) %>% as.matrix() %>% as.data.frame()
  diag(asv_bctab_tmp)<-NA #set diagonal (same-to-same comparison) to NA
  
  # Save BC distances in list
  asv_bctab_list[[i]]<-asv_bctab_tmp
}

# Reduce list to a single dataframe
asv_bctab_combined<-Reduce("+",asv_bctab_list) #sum up the corresponding elements
asv_bctab<-asv_bctab_combined/length(asv_bctab_list) #calculate mean of each element

# Write BC distance dataframe
asv_bctab<-asv_bctab %>%
  rownames_to_column("id")
fwrite(asv_bctab,"../data/asv_bctab.tsv",sep="\t",na="NA")
