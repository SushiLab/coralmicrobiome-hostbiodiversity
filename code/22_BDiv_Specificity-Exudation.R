## Project: Aquarium Microbiome
## Script purpose: Explore host-specificity and exudation
## Date: November 2023
## Author: Fabienne Wiederkehr

# Load packages
require(data.table)
require(tidyverse)
require(ggforce)
require(ggpubr)
require(patchwork)
require(biotools)

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
  left_join(rename_all(dplyr::select(metadat,id,fragment_id,species,phase,tank,sample,sampling_date),~paste("sample_1",.,sep="_")),by=c("sample_1"="sample_1_id")) %>% #append metadata by sample_1
  left_join(rename_all(dplyr::select(metadat,id,fragment_id,species,phase,tank,sample,sampling_date),~paste("sample_2",.,sep="_")),by=c("sample_2"="sample_2_id")) #append metadata by sample_2

# Prepare host to tank distance dataframe
dist_hosttotank<-asv_bc_dist %>%
  filter(sample_1_phase==sample_2_phase) %>% #retain distances if phase is the same
  filter(sample_1_tank==sample_2_tank) %>% #retain distances if tank is the same
  filter(sample_1_sample %in% c("mucus","biofilm","tissue")&sample_2_sample=="tank water") %>%
  group_by(sample_1) %>%
  summarise(mean(bc_dissim),sd(bc_dissim)) %>%
  mutate(ymin=`mean(bc_dissim)`-`sd(bc_dissim)`,ymax=`mean(bc_dissim)`+`sd(bc_dissim)`) %>%
  mutate(ymax=ifelse(ymax>1,1,ymax)) %>%
  rename(bc_hosttotank=`mean(bc_dissim)`) %>%
  mutate(sample_1=gsub("h_","",sample_1))

# Prepare exudates to incubation control distance dataframe
dist_exudatetocontrol<-asv_bc_dist %>%
  filter(sample_1_sampling_date==sample_2_sampling_date) %>% #retain distances if sampling_date is the same
  filter(sample_1_sample=="exudates"&sample_2_sample=="incubation control") %>%
  group_by(sample_1) %>%
  summarise(mean(bc_dissim),sd(bc_dissim)) %>%
  mutate(xmin=`mean(bc_dissim)`-`sd(bc_dissim)`,xmax=`mean(bc_dissim)`+`sd(bc_dissim)`) %>%
  mutate(xmax=ifelse(xmax>1,1,xmax)) %>%
  rename(bc_exudatetocontrol=`mean(bc_dissim)`) %>%
  left_join(dplyr::select(metadat,id,species,treatment,treatment_binary,phase),by=c("sample_1"="id")) %>%
  mutate(sample_1=gsub("w_","",sample_1))

# Prepare dataframe to plot
toplot<-dist_exudatetocontrol %>%
  left_join(dist_hosttotank,by="sample_1") %>% as.data.frame() %>%
  na.omit()

# MANOVA
man_res<-manova(cbind(bc_exudatetocontrol,bc_hosttotank)~species,data=toplot)
man_sum<-summary(man_res)

# Pairwise comparisons
mvpaircomp(man_res,"species",test="Wilks",adjust="holm")

# Plot
g<-ggplot(data=toplot,aes(x=bc_exudatetocontrol,y=bc_hosttotank)) +
  stat_ellipse(geom="polygon",alpha=0.5,aes(fill=species),level=0.5) +
  geom_point(size=3,aes(colour=species)) +
  geom_vline(aes(xintercept=mean(bc_exudatetocontrol)),linetype="dashed") +
  geom_hline(aes(yintercept=mean(bc_hosttotank)),linetype="dashed") +
  xlab("Bray-Curtis distance\nexudate to control") +
  ylab("Bray-Curtis distance\nhost to tank") +
  xlim(min(toplot$xmin),1) +
  ylim(min(toplot$ymin),1) +
  scale_colour_manual(values=colours) +
  scale_fill_manual(values=colours,guide="none") +
  #coord_fixed() +
  labs(title="Testing for differences between species",subtitle=bquote("MANOVA: "~italic(p)~"value ="~.(signif(man_sum$stats[1,"Pr(>F)"],3))),colour="Species",fill="Species") +
  theme_bw() +
  theme(aspect.ratio=1)
ggsave(g,filename="../data/22_BDiv_specificity-exudation.png",width=12,height=6,dpi=500)

# Adding the standard deviation
g +
  geom_linerange(aes(x=bc_exudatetocontrol,ymin=ymin,ymax=ymax,colour=species)) +
  geom_linerange(aes(y=bc_hosttotank,xmin=xmin,xmax=xmax,colour=species))
ggsave(filename="../data/22_BDiv_specificity-exudation_sd.png",width=12,height=6,dpi=500)
