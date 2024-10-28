## Project: Aquarium Microbiome
## Script purpose: Compare communities
## Date: May 2024
## Author: Fabienne Wiederkehr

# Load packages
require(data.table)
require(tidyverse)
require(vegan)
require(ape)
require(pairwiseAdonis)

# Load data
metadat<-fread("../data/rarefied/metadat_rarefied.tsv",sep="\t",header=TRUE,data.table=FALSE)
asv_bctab<-fread("../data/asv_bctab.tsv",sep="\t",header=TRUE,data.table=FALSE) %>%
  column_to_rownames("id") %>%
  filter(rownames(.) %in% metadat$id) %>% as.matrix() %>% t() %>% as.data.frame() %>%
  filter(rownames(.) %in% metadat$id)
colours<-readRDS("../data/colour_palette.RDS")

# Set diagonal (same-to-same comparison) to 0
diag(asv_bctab)<-0

# Get principal coordinates
pcoa_res<-pcoa(D=asv_bctab %>% as.matrix() %>% as.dist()) #sqrt-transform to lessen impact of abundant ASVs
pcoa_df<-pcoa_res$vectors %>% as.data.frame() %>%
  rownames_to_column("id") %>%
  left_join(metadat,by="id") %>%
  mutate(species=ifelse(is.na(species),sample,species)) %>%
  mutate(species=factor(species,levels=c("Montipora digitata","Pocillopora verrucosa",
                                         "Sinularia sp.","Xenia sp.",
                                         "Caulerpa sp.","Peyssonnelia sp.",
                                         "Haliclona cnidata",
                                         "tank water","incubation control","microscopy-slide control",
                                         "extraction mock","PCR mock",
                                         "extraction blank","PCR blank","PBS control","Tris-NaCl control"))) %>%
  mutate(source=factor(source,levels=c("host-associated","exuded","free-living")))

# Plot community distances (PCoA) by sample type
ggplot(data=pcoa_df,aes(x=Axis.1,y=Axis.2,colour=species)) +
  geom_hline(yintercept=0,linetype=2) +
  geom_vline(xintercept=0,linetype=2) +
  geom_point(size=4,alpha=0.9) +
  scale_colour_manual(values=colours) +
  facet_wrap(~source) +
  labs(title="PCoA based on square-root transformed Bray-Curtis distances") +
  xlab(bquote("Axis 1 ("~.(round(pcoa_res$values$Relative_eig[1]*100,2))~"%)")) +
  ylab(bquote("Axis 2 ("~.(round(pcoa_res$values$Relative_eig[2]*100,2))~"%)")) +
  theme_bw()
ggsave(filename="../data/19_BDiv_PCoA_by-source.png",width=15,height=7,dpi=500)

# PERMANOVA
for (i in c("host-associated","exuded")){
  cat(i,"\n")
  # Subset
  metadat_perm<-metadat %>%
    filter(!is.na(species)) %>%
    filter(source==i)
  
  asv_bctab_perm<-asv_bctab %>%
    filter(rownames(.) %in% metadat_perm$id) %>% as.matrix() %>% t() %>% as.data.frame() %>%
    filter(rownames(.) %in% metadat_perm$id) %>% as.matrix() %>% as.dist()
  
  adonis_res<-adonis2(asv_bctab_perm~species,data=metadat_perm,permutations=1000,by="terms")
  print(adonis_res)
  
  # Perform post-hoc test for dissimilarity (ADONIS)
  pair_res<-pairwise.adonis2(asv_bctab_perm~species,data=metadat_perm,nperm=1000,p.adjust.method="holm")
  print(pair_res)
}

# PERMDISP
permutest(betadisper(asv_bctab_perm,metadat_perm$species),permutations=how(nperm=1000))
