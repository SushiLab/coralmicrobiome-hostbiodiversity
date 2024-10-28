## Project: Aquarium Microbiome
## Script purpose: Analyse the communities of pathway predictions by species for host-associated data only
## Date: July 2024
## Author: Fabienne Wiederkehr

# Load packages
require(data.table)
require(tidyverse)
require(vegan)
require(ape)
require(pairwiseAdonis)

# Load KEGG info
path_list_cat<-fread("../data/kegg/pathway_list_cat.csv",header=TRUE,data.table=FALSE) %>%
  filter(pathway_category=="Metabolism") %>% #keep only pathways associated with metabolism
  filter(!pathway_class=="Global and overview maps") #remove global and overview maps

for (i in c("exuded","host-associated")){
  # Load data
  metadat<-fread("../data/rarefied/metadat_rarefied.tsv",sep="\t",header=TRUE,data.table=FALSE) %>%
    filter(!is.na(species)) %>%
    filter(source==i)
  path_dat<-fread("../data/picrust_output/path_abun_unstrat.tsv.gz",sep="\t",header=TRUE,data.table=FALSE) %>%
    column_to_rownames("pathway") %>% t() %>% as.data.frame() %>%
    filter(rownames(.) %in% metadat$id) %>%
    select(all_of(intersect(colnames(.),path_list_cat$pathway_id)))
  colours<-readRDS("../data/colour_palette.RDS")
  
  # Compare communities based on Bray-Curtis distance
  path_bctab<-vegdist(sqrt(path_dat))
  
  # Get principal coordinates
  pcoa_res<-pcoa(D=path_bctab)
  pcoa_df<-pcoa_res$vectors %>% as.data.frame() %>%
    rownames_to_column("id") %>%
    left_join(metadat,by="id") %>%
    mutate(species=factor(species,levels=c("Montipora digitata","Pocillopora verrucosa",
                                           "Sinularia sp.","Xenia sp.",
                                           "Caulerpa sp.","Peyssonnelia sp.",
                                           "Haliclona cnidata")))
  
  # PERMDISP
  permutest(betadisper(path_bctab,metadat$species),permutations=how(nperm=1000))
  
  # PERMANOVA
  perm_res<-adonis2(path_bctab~species,data=metadat,permutations=1000,by="terms")
  
  # Perform post-hoc test for dissimilarity (ADONIS)
  posthoc_adonis<-pairwise.adonis2(path_bctab~species,data=metadat,nperm=1000,p.adjust.method="holm")
  
  # Plot community distances (PCoA)
  ggplot(data=pcoa_df,aes(x=Axis.1,y=Axis.2)) +
    geom_hline(yintercept=0,linetype=2) +
    geom_vline(xintercept=0,linetype=2) +
    geom_point(size=4,alpha=0.9,aes(colour=species)) +
    scale_colour_manual(values=colours) +
    scale_fill_manual(values=colours,guide="none") +
    labs(title=paste("PCoA based on square-root transformed Bray-Curtis distances on",i,"data",sep=" "),subtitle=bquote("PERMANOVA: "~italic(p)~"value ="~.(signif(perm_res$`Pr(>F)`[1],3))),colour="Species") +
    xlab(bquote("Axis 1 ("~.(round(pcoa_res$values$Relative_eig[1]*100,2))~"%)")) +
    ylab(bquote("Axis 2 ("~.(round(pcoa_res$values$Relative_eig[2]*100,2))~"%)")) +
    theme_bw()
  
  ggsave(filename=paste("../data/38_FunPred_pathway_species_",i,".png",sep=""),width=unit(12,"cm"),height=unit(6,"cm"))
}
