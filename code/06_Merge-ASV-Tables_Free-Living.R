## Project: Aquarium Microbiome
## Script purpose: Merge ASV tables from raw_free-living_1 and raw_free-living_2
## Date: August 2023
## Author: Fabienne Wiederkehr

# Load packages
require(data.table)
require(tidyverse)
require(dada2)

# Load data and filter out eukaryotic, chloroplast, and mitochondrial reads
asv_dat1<-fread("../data/raw_free-living_1.asvs.tsv",sep="\t",header=TRUE,data.table=FALSE) %>%
  filter(!grepl("domain;unclassified|domain;Eukaryota|family;Mitochondria|order;Chloroplast",tax)) %>% #remove reads unclassified at domain level (removes short ASVs), eukaryotic reads, mitochondria and chloroplasts
  select(-c(asv,otu,uparse_info,tax)) %>%
  column_to_rownames("seq") %>%
  as.matrix() %>% t()

asv_dat2<-fread("../data/raw_free-living_2.asvs.tsv",sep="\t",header=TRUE,data.table=FALSE) %>%
  filter(!grepl("domain;unclassified|domain;Eukaryota|family;Mitochondria|order;Chloroplast",tax)) %>% #remove reads unclassified at domain level (removes short ASVs), eukaryotic reads, mitochondria and chloroplasts
  select(-c(asv,otu,uparse_info,tax)) %>%
  column_to_rownames("seq") %>%
  as.matrix() %>% t()

# Merge the two ASV tables
asv_dat<-mergeSequenceTables(asv_dat1,asv_dat2)
asv_dat_collapsed<-collapseNoMismatch(asv_dat) #combine together sequences that are identical up to shifts and/or length

asv_dat_collapsed_df<-asv_dat_collapsed %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("seq")

# Rename columns as ASVs
rownames(asv_dat_collapsed_df)<-paste("asv",str_pad(1:nrow(asv_dat_collapsed_df),width=6,pad=0),sep="_")
asv_dat_renamed<-asv_dat_collapsed_df %>%
  rename_all(~stringr::str_replace(.,"^WIED23-._","")) %>% #remove prefix
  rename_all(~stringr::str_replace(.,"_METAB","")) %>% #remove suffix
  rownames_to_column("asv")
fwrite(asv_dat_renamed,"../data/asv_dat_original_free-living.tsv",sep="\t",na="NA")

# Generate ASV abundance table
asv_abtab<-asv_dat_renamed %>%
  column_to_rownames("asv") %>%
  select(-seq) %>%
  t() %>% as.data.frame() %>% #transpose table
  rownames_to_column("id") #move rownames to column
fwrite(asv_abtab,"../data/asv_abtab_original_free-living.tsv",sep="\t",na="NA")
