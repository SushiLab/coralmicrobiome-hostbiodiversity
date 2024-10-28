## Project: Aquarium Microbiome
## Script purpose: Explore host-specificity and exudation as a function of the complexity of the reef community
## Date: January 2024
## Author: Fabienne Wiederkehr

# Load packages
require(data.table)
require(tidyverse)
require(ggforce)
require(ggpubr)
require(patchwork)

# Load data
metadat<-fread("../data/rarefied/metadat_rarefied.tsv",sep="\t",header=TRUE,data.table=FALSE) %>%
  mutate(species=factor(species,levels=c("Montipora digitata","Pocillopora verrucosa",
                                         "Sinularia sp.","Xenia sp.",
                                         "Caulerpa sp.","Peyssonnelia sp.",
                                         "Haliclona cnidata")))
asv_bctab<-fread("../data/asv_bctab.tsv",sep="\t",header=TRUE,data.table=FALSE) %>%
  column_to_rownames("id") %>%
  filter(rownames(.) %in% metadat$id)
colours<-readRDS("../data/colour_palette.RDS")

# Transform into 3-column dataframe
asv_bc_dist<-asv_bctab %>%
  rownames_to_column("sample_1") %>%
  pivot_longer(-sample_1,names_to="sample_2",values_to="bc_dissim") %>% #names_to puts the column names in a second variable, values_to puts the values in a third variable
  filter(!is.na(bc_dissim)) %>% #remove NAs
  left_join(rename_all(dplyr::select(metadat,id,fragment_id,species,phase,tank,sample,sampling_date,treatment_binary),~paste("sample_1",.,sep="_")),by=c("sample_1"="sample_1_id")) %>% #append metadata by sample_1
  left_join(rename_all(dplyr::select(metadat,id,fragment_id,species,phase,tank,sample,sampling_date),~paste("sample_2",.,sep="_")),by=c("sample_2"="sample_2_id")) #append metadata by sample_2

# Prepare host to tank distance dataframe
dist_hosttotank<-asv_bc_dist %>%
  filter(sample_1_phase=="P2"&sample_2_phase=="P2") %>% #retain phase 2 distances
  filter(sample_1_tank==sample_2_tank) %>% #retain distances if tank is the same
  filter(sample_1_sample %in% c("mucus","biofilm","tissue")&sample_2_sample=="tank water") %>%
  dplyr::select(bc_dissim,contains("sample_1")) %>%
  rename_all(~stringr::str_replace(.,"^sample_1_","")) #remove prefix

# Prepare exudates to incubation control distance dataframe
dist_exudatetocontrol<-asv_bc_dist %>%
  filter(sample_1_phase=="P2"&sample_2_phase=="P2") %>% #retain phase 2 distances
  filter(sample_1_sampling_date==sample_2_sampling_date) %>% #retain distances if sampling_date is the same
  filter(sample_1_sample=="exudates"&sample_2_sample=="incubation control") %>%
  dplyr::select(bc_dissim,contains("sample_1")) %>%
  rename_all(~stringr::str_replace(.,"^sample_1_","")) #remove prefix

# Compare BC host-to-tank distances
g<-ggplot(data=dist_hosttotank,aes(x=species,y=bc_dissim,fill=treatment_binary,colour=treatment_binary)) +
  geom_violin(alpha=0.5,draw_quantiles=0.5) +
  geom_sina(aes(colour=species),size=0.8) +
  ylab("Bray-Curtis distance\nhost to tank") +
  ylim(0.85,1) + #align with script 53 
  #labs(fill="Reef community",colour="Reef community") +
  scale_colour_manual(values=colours) +
  scale_fill_manual(values=colours) +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        legend.position="none")
g_stats<-compare_means(bc_dissim~treatment_binary,data=dist_hosttotank,group.by="species",p.adjust.method="holm") %>% as.data.frame()
#wilcox.test(bc_dissim~treatment_binary,data=dist_hosttotank %>% filter(species=="Haliclona cnidata")) #for group means

# Compare BC exudate-to-control distances
h<-ggplot(data=dist_exudatetocontrol,aes(x=species,y=bc_dissim,fill=treatment_binary,colour=treatment_binary)) +
  geom_violin(alpha=0.5,draw_quantiles=0.5) +
  geom_sina(aes(colour=species),size=0.8) +
  ylab("Bray-Curtis distance\nexudate to control") +
  ylim(0.48,1) + #align with script 53 
  labs(fill="Reef community",colour="Reef community") +
  scale_colour_manual(values=colours) +
  scale_fill_manual(values=colours) +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        legend.position="none")
h_stats<-compare_means(bc_dissim~treatment_binary,data=dist_exudatetocontrol,group.by="species",p.adjust.method="holm") %>% as.data.frame()

# Combine plots
patchwork::wrap_plots(list(g,h),guides="collect",ncol=1) +
  plot_annotation(title="Differential distances as a function of the complexity of the reef community")
ggsave(filename="../data/23_BDiv_complexity_specificity-exudation.png",width=8,height=8,dpi=500)

# Stats
stats<-levels(metadat$species) %>% as.data.frame() %>%
  rename("species"=".") %>%
  left_join(dplyr::select(g_stats,species,p.adj),by="species") %>%
  left_join(dplyr::select(h_stats,species,p.adj),by="species") %>%
  rename("p.adj.hosttotank"="p.adj.x") %>%
  rename("p.adj.exudatetocontrol"="p.adj.y")
