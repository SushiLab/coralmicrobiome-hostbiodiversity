## Project: Aquarium Microbiome
## Script purpose: Compute shared across treatment and source by host
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
asv_shared<-fread("../data/asv_shared.tsv",sep="\t",header=TRUE,data.table=FALSE) %>%
  column_to_rownames("id")

# Initialise variables
stats_dist<-NULL

# Loop according to experimental design
for (k in levels(metadat$source)){
  for (i in levels(metadat$treatment_binary)){
    cat("Treatment",i,"\n")
    
    # Subset metadata
    metadat_red<-metadat %>%
      filter(treatment_binary==i&source==k)
    
    # Subset BC distances
    asv_shared_red<-asv_shared %>%
      filter(rownames(.) %in% metadat_red$id) %>% as.matrix() %>% t() %>% as.data.frame() %>%
      filter(rownames(.) %in% metadat_red$id)
    
    # Summarise (mean) distance between species
    asv_dist<-asv_shared_red %>%
      rownames_to_column("id") %>%
      left_join(select(metadat,id,species),by="id") %>% #append metadata
      group_by(species) %>%
      summarise(across(where(is.numeric),~mean(.x,na.rm=TRUE))) %>%
      column_to_rownames("species") %>% t() %>% as.data.frame() %>% #transform dataframe and repeat process
      rownames_to_column("id") %>%
      left_join(select(metadat,id,species),by="id") %>%
      group_by(species) %>%
      summarise(across(where(is.numeric),~mean(.x,na.rm=TRUE))) %>%
      column_to_rownames("species")
    diag(asv_dist)<-NA
    
    # Define edges (distance between species)
    edge_df<-asv_dist %>%
      mutate(from=rownames(.)) %>%
      gather(to,dist,rownames(.)) %>%
      mutate(interaction=interaction(from,to,sep=";")) %>%
      column_to_rownames("interaction") %>%
      select(dist)
    
    # Write file
    ed<-edge_df %>%
      rownames_to_column("interaction") %>%
      separate(interaction,into=c("from","to"),sep=";") %>%
      rename(shared=dist)
    fwrite(ed,paste("../data/28_ASVs_ed_",k,i,".tsv",sep=""),sep="\t",na="NA")
    
    # Transform into 3-column dataframe
    tmp<-asv_shared_red %>%
      rownames_to_column("sample_1") %>%
      pivot_longer(-sample_1,names_to="sample_2",values_to="shared_dist") %>% #names_to puts the column names in a second variable, values_to puts the values in a third variable
      filter(!is.na(shared_dist)) %>% #remove NAs
      left_join(rename_all(dplyr::select(metadat,id,species,treatment_binary),~paste("sample_1",.,sep="_")),by=c("sample_1"="sample_1_id")) %>% #append metadata by sample_1
      left_join(rename_all(dplyr::select(metadat,id,species),~paste("sample_2",.,sep="_")),by=c("sample_2"="sample_2_id")) %>% #append metadata by sample_2
      filter(!sample_1_species==sample_2_species) %>%
      unite(comparison,sample_1_species,sample_2_species,sep=" vs ") %>%
      dplyr::select(shared_dist,comparison,sample_1_treatment_binary) %>%
      rename_all(~stringr::str_replace(.,"^sample_1_","")) #remove prefix

    stats_dist<-stats_dist %>%
      bind_rows(tmp)
  }
  
  # Compare treatments
  edge_stats<-compare_means(shared_dist~treatment_binary,data=stats_dist,group.by="comparison",p.adjust.method="holm") %>% as.data.frame()
  fwrite(edge_stats,paste("../data/28_ASVs_stats_ed_",k,".tsv",sep=""),sep="\t",na="NA")
} 
