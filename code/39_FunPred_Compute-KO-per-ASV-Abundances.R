## Project: Aquarium Microbiome
## Script purpose: Compute KO abundances averaged over the 50 picrust output files (~rarefaction files) for the differentially abundant ASVs
## Date: August 2024
## Author: Fabienne Wiederkehr

# Load packages
require(data.table)
require(tidyverse)
require(vegan)

# Load paths
ko_abtab_paths<-list.files("../data/picrust_output/",pattern="asv_abtab_rarefied_",full.names=TRUE)

# Load data
asv_da<-fread("../data/asv_diffab_lfc.tsv",sep="\t",header=TRUE,data.table=FALSE) %>%
  select(asv)

# Get a list of all KOs
ko_list<-fread(paste(ko_abtab_paths[1],"/KO_predicted.tsv.gz",sep=""),sep="\t",header=TRUE,data.table=FALSE) %>%
  filter(sequence %in% asv_da$asv) %>%
  column_to_rownames("sequence") %>% as.matrix() %>% t() %>% as.data.frame() %>%
  rownames_to_column("ko_id") %>%
  select(ko_id)
for (i in 2:length(ko_abtab_paths)){
  # Load data and extract KOs
  ko_abtab<-fread(paste(ko_abtab_paths[i],"/KO_predicted.tsv.gz",sep=""),sep="\t",header=TRUE,data.table=FALSE) %>%
    filter(sequence %in% asv_da$asv) %>%
    column_to_rownames("sequence") %>% as.matrix() %>% t() %>% as.data.frame() %>%
    rownames_to_column("ko_id") %>%
    select(ko_id)
  
  # Save KOs in list
  ko_list<-intersect(ko_abtab,ko_list)
}

# Compute KO abundances
ko_abtab_list<-list(NULL)
for (i in 1:length(ko_abtab_paths)){
  # Load data
  ko_abtab<-fread(paste(ko_abtab_paths[i],"/KO_predicted.tsv.gz",sep=""),sep="\t",header=TRUE,data.table=FALSE) %>%
    filter(sequence %in% asv_da$asv) %>%
    column_to_rownames("sequence") %>% as.matrix() %>% t() %>% as.data.frame() %>%
    rownames_to_column("ko_id")
  
  # Ensure KOs are the same in each file
  ko_dat<-ko_list %>%
    left_join(ko_abtab,by="ko_id") %>%
    column_to_rownames("ko_id")
  
  # Save KOs in list
  ko_abtab_list[[i]]<-ko_dat
}

# Reduce list to a single dataframe
ko_abtab_combined<-Reduce("+",ko_abtab_list) #sum up the corresponding elements
ko_abtab_mean<-ko_abtab_combined/length(ko_abtab_list) #calculate mean of each element

# Write file
ko_abtab_mean<-ko_abtab_mean %>%
  rownames_to_column("ko_id")
fwrite(ko_abtab_mean,"../data/ko_abtab_asv.tsv",sep="\t",na="NA")
