## Project: Aquarium Microbiome
## Script purpose: Compute number of unique ASVs across treatment and source by host
## Date: April 2024
## Author: Fabienne Wiederkehr

# Load packages
require(data.table)
require(tidyverse)
require(vegan)
require(ComplexUpset)

# Load paths
asv_abtab_paths<-list.files("../data/rarefied/",pattern="asv_abtab_rarefied",full.names=TRUE)

# Load data
metadat<-fread("../data/rarefied/metadat_rarefied.tsv",sep="\t",header=TRUE,data.table=FALSE) %>%
  filter(!is.na(species)) %>%
  filter(phase=="P2") %>%
  mutate(source=factor(source,levels=c("host-associated","exuded"))) %>%
  mutate(treatment_binary=factor(treatment_binary,levels=c("biodiverse","degraded"))) %>%
  mutate(species=factor(species,levels=c("Caulerpa sp.","Haliclona cnidata","Montipora digitata","Sinularia sp.","Pocillopora verrucosa","Peyssonnelia sp.","Xenia sp.")))

# Get set names
set<-rev(levels(metadat$species)) #pull species names and reverse order (upset() fills rows bottom-to-top)

# Initialise variables
node_list<-list(NULL)

# Loop according to experimental design
for (k in levels(metadat$source)){
  for (i in levels(metadat$treatment_binary)){
    cat("Treatment",i,"\n")
    
    # Iterate over each rarefaction file
    for (m in 1:length(asv_abtab_paths)){
      cat("\tIteration",m,"/",length(asv_abtab_paths),"\n")
      
      # Load data
      asv_abtab<-fread(asv_abtab_paths[m],sep="\t",header=TRUE,data.table=FALSE) %>%
        filter(id %in% metadat$id)
      
      # Subsample
      metadat_red<-metadat %>%
        filter(treatment_binary==i&source==k)
      asv_abtab_red<-metadat_red %>%
        select("id") %>%
        left_join(asv_abtab,by="id")
      
      # Transform dataframe
      asv_abtab_long<-asv_abtab_red %>%
        left_join(dplyr::select(metadat_red,id,species),by="id") %>% #append metadata
        group_by(species) %>%
        summarise(across(where(is.numeric),~sum(.x,na.rm=TRUE))) %>%
        column_to_rownames("species") %>% t() %>% as.data.frame()
      
      # Convert to presence-absence
      df_combined<-asv_abtab_long %>%
        mutate_if(is.numeric,as.logical)
      
      # Remove all-0 ASVs
      toplot<-df_combined %>%
        as.data.frame() %>%
        mutate_if(is.logical,as.numeric) %>% #transform TRUE/FALSE to 1/0
        filter(rowSums(.)>0) %>% #remove all-0 ASVs
        mutate_if(is.numeric,as.logical) #transform 1/0 to TRUE/FALSE
      
      # Define nodes (unqiue ASVs)
      node_list[[m]]<-upset_data(toplot,set)$sizes$exclusive_intersection %>% as.data.frame() %>%
        rename("unique_asvs"=".") %>%
        filter(rownames(.) %in% set) %>%
        arrange(match(rownames(.),set))
    }
    
    # Summarise (mean) over permutations
    node_combined<-Reduce(`+`,node_list)/length(node_list)
    va<-node_combined %>%
      rownames_to_column("species")
    fwrite(va,paste("../data/29_ASVs_va_",k,i,".tsv",sep=""),sep="\t",na="NA")
  }
}
