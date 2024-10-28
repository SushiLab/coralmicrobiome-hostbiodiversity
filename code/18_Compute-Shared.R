## Project: Aquarium Microbiome
## Script purpose: Compute number of shared ASVs
## Date: May 2024
## Author: Fabienne Wiederkehr

# Load packages
require(data.table)
require(tidyverse)
require(vegan)

# Load paths
asv_abtab_paths<-list.files("../data/rarefied/",pattern="asv_abtab_rarefied",full.names=TRUE)

# Load metadata
metadat<-fread("../data/rarefied/metadat_rarefied.tsv",sep="\t",header=TRUE,data.table=FALSE)

# Compute number of shared ASVs
asv_shared_list<-list(NULL)
for (i in 1:length(asv_abtab_paths)){
  # Load data
  asv_abtab<-fread(asv_abtab_paths[i],sep="\t",header=TRUE,data.table=FALSE) %>%
    column_to_rownames("id") %>%
    filter(rownames(.) %in% metadat$id)
  
  # Compute number of shared ASVs
  asv_shared_tmp<-designdist(asv_abtab,method="J",terms="binary") %>% as.matrix() %>% as.data.frame() #diagonal 0 by default
  diag(asv_shared_tmp)<-NA
  
  # Save number of shared ASVs in list
  asv_shared_list[[i]]<-asv_shared_tmp
}

# Reduce list to a single dataframe
asv_shared_combined<-Reduce("+",asv_shared_list) #sum up the corresponding elements
asv_shared<-asv_shared_combined/length(asv_shared_list) #calculate mean of each element

# Write richness dataframe
asv_shared<-asv_shared %>%
  rownames_to_column("id")
fwrite(asv_shared,"../data/asv_shared.tsv",sep="\t",na="NA")
