## Project: Aquarium Microbiome
## Script purpose: Compile sequence files
## Date: July 2024
## Author: Fabienne Wiederkehr

# Load packages
require(data.table)
require(tidyverse)
require(TAF)
require(vegan)

# Load paths
asv_abtab_paths<-list.files("../data/rarefied/",pattern="asv_abtab_rarefied",full.names=TRUE)

# Load data
metadat<-fread("../data/rarefied/metadat_rarefied.tsv",sep="\t",header=TRUE,data.table=FALSE)
asv_dat<-fread("../data/asv_dat_cleaned.tsv",sep="\t",header=TRUE,data.table=FALSE)

mkdir("../data/picrust_input")
# Write fasta files for all 
for (i in 1:length(asv_abtab_paths)){
  # Load data
  asv_abtab<-fread(asv_abtab_paths[i],sep="\t",header=TRUE,data.table=FALSE) %>%
    column_to_rownames("id") %>%
    filter(rownames(.) %in% metadat$id) %>% as.matrix() %>% t() %>% as.data.frame()
  
  # Append sequence information
  asv_seq<-data.frame(asv=rownames(asv_abtab)) %>%
    left_join(asv_dat,by="asv")

  # Create and export fasta file
  asv_headers<-paste0(">",asv_seq$asv)
  asv_seqs<-asv_seq$seq
  asv_fasta<-c(rbind(asv_headers,asv_seqs))
  write.table(asv_fasta,file=paste("../data/picrust_input/asv_abtab_rarefied_",i,".fasta",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
}
