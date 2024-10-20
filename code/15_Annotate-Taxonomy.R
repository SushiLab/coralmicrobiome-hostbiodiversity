## Project: Aquarium Microbiome
## Script purpose: Add taxonomic information
## Date: October 2023
## Author: Fabienne Wiederkehr

# Load packages
require(data.table)
require(tidyverse)
require(DECIPHER)
require(dada2)

# Load Data
asv_dat<-fread("../data/asv_dat_cleaned.tsv",sep="\t",header=TRUE,data.table=FALSE)

# Read FASTA file
seqs<-DNAStringSet(x=asv_dat$seq)
names(seqs)<-asv_dat$asv

# Run assignTaxonomy and parse
taxa<-assignTaxonomy(seqs,"../data/silva_nr99_v138.1_train_set.fa",minBoot=80) #download from SILVA DB
unname(taxa)

# Transform dataframe and write taxinfo file
asv_taxinfo<-taxa %>% as.data.frame() %>%
  rownames_to_column("seq") %>%
  left_join(select(asv_dat,asv,seq),by="seq")
fwrite(asv_taxinfo,"../data/asv_dat_taxinfo.tsv",sep="\t",na="NA")
