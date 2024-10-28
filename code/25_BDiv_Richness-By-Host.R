## Project: Aquarium Microbiome
## Script purpose: Compute richness across treatment and source by host
## Date: April 2024
## Author: Fabienne Wiederkehr

# Load packages
require(data.table)
require(tidyverse)
require(ggpubr)

# Load data
metadat<-fread("../data/rarefied/metadat_rarefied.tsv",sep="\t",header=TRUE,data.table=FALSE) %>%
  filter(!is.na(species)) %>%
  filter(phase=="P2") %>%
  mutate(source=factor(source,levels=c("host-associated","exuded"))) %>%
  mutate(treatment_binary=factor(treatment_binary,levels=c("biodiverse","degraded"))) %>%
  mutate(species=factor(species,levels=c("Montipora digitata","Pocillopora verrucosa",
                                         "Sinularia sp.","Xenia sp.",
                                         "Caulerpa sp.","Peyssonnelia sp.",
                                         "Haliclona cnidata")))
asv_richtab<-fread("../data/asv_richtab.tsv",sep="\t",header=TRUE,data.table=FALSE) %>%
  column_to_rownames("id") %>%
  select(richness=hill_one)

# Get sample size
table(metadat$treatment_binary,metadat$species,metadat$source) #%>% as.data.frame() %>% rename(treatment_binary=Var1,species=Var2,source=Var3,n=Freq)

# Initialise variables
stats_rich<-NULL

# Loop according to experimental design
for (k in levels(metadat$source)){
  for (i in levels(metadat$treatment_binary)){
    cat("Treatment",i,"\n")
    
    # Subset metadata
    metadat_red<-metadat %>%
      filter(treatment_binary==i&source==k)
    
    # Subset BC distances
    asv_richtab_red<-asv_richtab %>%
      filter(rownames(.) %in% metadat_red$id)
    
    # Summarise (mean) distance between species
    asv_rich<-asv_richtab_red %>%
      rownames_to_column("id") %>%
      left_join(select(metadat,id,species),by="id") %>% #append metadata
      group_by(species) %>%
      summarise(across(where(is.numeric),~mean(.x,na.rm=TRUE))) %>%
      column_to_rownames("species")
    
    # Write file
    va<-asv_rich %>%
      rownames_to_column("species")
    fwrite(va,paste("../data/25_BDiv_va_",k,i,".tsv",sep=""),sep="\t",na="NA")
    
    # Prepare stats dataframe
    tmp<-asv_richtab_red %>%
      rownames_to_column("id") %>%
      left_join(select(metadat,id,species,treatment_binary),by="id")

    stats_rich<-stats_rich %>%
      bind_rows(tmp)
  }
  
  # Compare treatments
  node_stats<-compare_means(richness~treatment_binary,data=stats_rich,group.by="species",p.adjust.method="holm") %>% as.data.frame()
  fwrite(node_stats,paste("../data/25_BDiv_stats_va_",k,".tsv",sep=""),sep="\t",na="NA")
} 
