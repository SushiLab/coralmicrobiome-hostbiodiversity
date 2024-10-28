## Project: Aquarium Microbiome
## Script purpose: Compare the host-associated and free-living microbial communities of hosts living in complex vs reduced reef communities
## Date: November 2023
## Author: Fabienne Wiederkehr

# Load packages
require(data.table)
require(tidyverse)
require(vegan)
#require(ape)

# Load data
metadat<-fread("../data/rarefied/metadat_rarefied.tsv",sep="\t",header=TRUE,data.table=FALSE) %>%
  filter(!is.na(species)) %>%
  filter(phase=="P2") %>%
  mutate(source=factor(source,levels=c("host-associated","exuded"))) %>%
  mutate(species=factor(species,levels=c("Montipora digitata","Pocillopora verrucosa",
                                         "Sinularia sp.","Xenia sp.",
                                         "Caulerpa sp.","Peyssonnelia sp.",
                                         "Haliclona cnidata")))
asv_bctab<-fread("../data/asv_bctab.tsv",sep="\t",header=TRUE,data.table=FALSE) %>%
  column_to_rownames("id")

# Set diagonal (same-to-same comparison) to 0
diag(asv_bctab)<-0

# Analyse by species and source
spp<-levels(metadat$species)
source<-levels(metadat$source)
adonis_df<-NULL

for (k in source){
  for (i in spp){ #for all species, do:
    # Subset metadat
    metadat_red<-metadat %>%
      filter(species==i) %>%
      filter(source==k)
    
    # Subset BC distances
    asv_bctab_red<-asv_bctab %>%
      filter(rownames(.) %in% metadat_red$id) %>% as.matrix() %>% t() %>% as.data.frame() %>%
      filter(rownames(.) %in% metadat_red$id)
    
    # Get principal coordinates
    dist<-asv_bctab_red %>% as.matrix() %>% as.dist()
    
    # PERMANOVA
    adonis_res<-adonis2(dist~treatment_binary,data=metadat_red,permutations=1000,by="terms")
    tmp<-adonis_res %>%
      as.data.frame() %>%
      dplyr::select(R2,p_value=`Pr(>F)`) %>%
      rownames_to_column("term") %>%
      mutate(species=i) %>%
      mutate(source=k)
    adonis_df<-adonis_df %>% bind_rows(tmp)
  }
}

# Partition dissimilarity among sources of variation
toplot<-adonis_df %>%
  filter(term %in% c("treatment_binary")) %>%
  mutate(term=case_when(term=="treatment_binary"~"Treatment")) %>%
  mutate(term=fct_relevel(term,"Treatment")) %>%
  mutate(p_value_dich=ifelse(p_value<0.05,"sign.","n.s.")) %>%
  mutate(source=factor(source,levels=c("host-associated","exuded"))) %>%
  mutate(species=factor(species,levels=c("Montipora digitata","Pocillopora verrucosa",
                                         "Sinularia sp.","Xenia sp.",
                                         "Caulerpa sp.","Peyssonnelia sp.",
                                         "Haliclona cnidata")))

# Plot PERMANOVA results
ggplot(data=toplot) +
  geom_bar(stat="identity",position="dodge",aes(x=species,y=R2*100,fill=source,colour=p_value_dich),alpha=0.5,size=2) +
  geom_text(stat="identity",position=position_dodge(width=0.9),aes(x=species,y=R2*100,label=round(p_value,3),group=source),vjust=-0.7,size=3) + 
  scale_colour_manual(values=c("#FFBB78","#41b6c4","#cb181d","#EAEAEA"),aesthetics=c("fill","colour"),breaks=c("host-associated","exuded","sign.","n.s.")) +
  theme_bw() +
  ylab("Variance explained (%)") +
  labs(fill=bquote("Source and"~italic(p)~"value"),
       title="Microbiomes of hosts living in complex vs reduced reef communities",
       subtitle=bquote("Permutational MANOVA (sign.:"~italic(p)~"value < 0.05)")) +
  guides(colour="none") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        axis.title.x=element_blank())
ggsave(filename="../data/21_BDiv_PERMANOVA_complexity.png",width=10,height=7,dpi=500)
