## Project: Aquarium Microbiome
## Script purpose: Download KEGG information from the API
## Date: July 2024
## Author: Fabienne Wiederkehr

# Load packages
require(tidyverse)
require(TAF)

# Load KEGG info (pathway, ko, and their link file)
pathway_list<-read.table("https://rest.kegg.jp/list/pathway",header=FALSE,sep="\t",stringsAsFactors=FALSE) %>%
  mutate(V1=sub("^map","ko",V1)) #substitute map with ko (to match PICRUSt2 output)
ko_list<-read.table("https://rest.kegg.jp/list/ko",header=FALSE,sep="\t",stringsAsFactors=FALSE)
link<-read.table("https://rest.kegg.jp/link/ko/pathway",header=FALSE,sep="\t",stringsAsFactors=FALSE)

# Assign column names
colnames(pathway_list)<-c("pathway_id","pathway_description")
colnames(ko_list)<-c("ko_id","ko_description")
colnames(link)<-c("pathway_id","ko_id")

# Remove prefixes in the link file
ko_pathway_link<-link %>%
  mutate(pathway_id=sub("^path:","",pathway_id)) %>%
  mutate(ko_id=sub("^ko:","",ko_id)) %>%
  filter(!grepl("map",pathway_id))

# Write files
mkdir("../data/kegg")
fwrite(pathway_list,"../data/kegg/pathway_list.tsv",sep="\t",na="NA")
fwrite(ko_list,"../data/kegg/ko_list.tsv",sep="\t",na="NA")
fwrite(ko_pathway_link,"../data/kegg/ko_pathway_link.tsv",sep="\t",na="NA")
