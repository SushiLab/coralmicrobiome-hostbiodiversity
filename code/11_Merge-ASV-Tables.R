## Project: Aquarium Microbiome
## Script purpose: Merge ASV tables from host-associated and free-living samples
## Date: August 2023
## Author: Fabienne Wiederkehr

# Load packages
require(data.table)
require(tidyverse)
require(dada2)

# Load data and filter out eukaryotic, chloroplast, and mitochondrial reads
asv_dat1<-fread("../data/asv_dat_cleaned_free-living.tsv",sep="\t",header=TRUE,data.table=FALSE) %>%
  select(-asv) %>%
  column_to_rownames("seq") %>%
  as.matrix() %>% t()

asv_dat2<-fread("../data/asv_dat_cleaned_host-associated.tsv",sep="\t",header=TRUE,data.table=FALSE) %>%
  select(-asv) %>%
  column_to_rownames("seq") %>%
  as.matrix() %>% t()

# Merge the two ASV tables
asv_dat<-mergeSequenceTables(asv_dat1,asv_dat2)
asv_dat_collapsed<-collapseNoMismatch(asv_dat) #combine together sequences that are identical up to shifts and/or length

# Remove singletons
asv_dat_nosing<-asv_dat_collapsed[,which(apply(asv_dat_collapsed,2,sum)>1)]
asv_dat_nosing<-asv_dat_nosing %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("seq")

# Rename columns as ASVs
rownames(asv_dat_nosing)<-paste("asv",str_pad(1:nrow(asv_dat_nosing),width=6,pad=0),sep="_")
asv_dat_renamed<-asv_dat_nosing %>%
  rownames_to_column("asv")
fwrite(asv_dat_renamed,"../data/asv_dat_cleaned.tsv",sep="\t",na="NA")

# Generate ASV abundance table
asv_abtab<-asv_dat_renamed %>%
  column_to_rownames("asv") %>%
  select(-seq) %>%
  t() %>% as.data.frame() %>% #transpose table
  rownames_to_column("id") #move rownames to column
fwrite(asv_abtab,"../data/asv_abtab_cleaned.tsv",sep="\t",na="NA")
