## Project: Aquarium Microbiome
## Script purpose: Rarefy data
## Date: October 2023
## Author: Fabienne Wiederkehr

# Load packages
require(data.table)
require(tidyverse)
require(TAF)
require(vegan)

# Load data
asv_abtab<-fread("../data/asv_abtab_cleaned.tsv",sep="\t",header=TRUE,data.table=FALSE) %>%
  filter(!grepl("Dic",id)) %>%
  column_to_rownames("id")
metadat<-fread("../data/metadat_cleaned.tsv",sep="\t",header=TRUE,data.table=FALSE) %>%
  filter(id %in% rownames(asv_abtab))

# Set rarefaction cutoff and generate 50 files
rr_cutoff<-2500

mkdir("../data/rarefied")
for (i in 1:50){
  # Rarefy
  asv_abtab_rr<-asv_abtab %>%
    filter(rowSums(.)>=rr_cutoff)
  asv_abtab_rr<-rrarefy(asv_abtab_rr,sample=rr_cutoff) #subsample
  asv_abtab_rr<-asv_abtab_rr[,which(apply(asv_abtab_rr,2,sum)>0)] %>% as.data.frame() %>% #remove all-0 (non-present when downsampling) ASVs
    rownames_to_column("id")
  fwrite(asv_abtab_rr,file=paste("../data/rarefied/asv_abtab_rarefied_",i,".tsv",sep=""),sep="\t")
}

# Generate reduced metadata
metadat_rr<-metadat %>%
  filter(id %in% asv_abtab_rr$id)
fwrite(metadat_rr,"../data/rarefied/metadat_rarefied.tsv",sep="\t",na="NA")
