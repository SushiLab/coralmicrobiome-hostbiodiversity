## Project: Aquarium Microbiome
## Script purpose: Correlate alpha-diversity changes of exuded vs host-associated microbial communities
## Date: January 2024
## Author: Fabienne Wiederkehr

# Load packages
require(data.table)
require(tidyverse)
require(vegan)
require(ggpubr)
require(patchwork)

# Load data
metadat<-fread("../data/rarefied/metadat_rarefied.tsv",sep="\t",header=TRUE,data.table=FALSE) %>%
  filter(grepl("exudates|mucus|biofilm|tissue",sample)) %>%
  mutate(species=factor(species,levels=c("Montipora digitata","Pocillopora verrucosa",
                                         "Sinularia sp.","Xenia sp.",
                                         "Caulerpa sp.","Peyssonnelia sp.",
                                         "Haliclona cnidata")))
asv_richtab<-fread("../data/asv_richtab.tsv",sep="\t",header=TRUE,data.table=FALSE) %>%
  filter(id %in% metadat$id)
colours<-readRDS("../data/colour_palette.RDS")

# Compute alpha-diversity dataframe
hill_df<-asv_richtab %>%
  left_join(dplyr::select(metadat,id,source,species,taxonomic_group,treatment_binary),by="id") %>%
  mutate(id=gsub("^._","",id)) #remove prefix

# Prepare dataframe to plot
toplot<-hill_df %>%
  pivot_wider(names_from=source,values_from=c(hill_zero,hill_one,hill_two)) %>%
  na.omit()

# Correlation plot
plot_zero<-ggplot(data=toplot,aes(x=`hill_zero_exuded`,y=`hill_zero_host-associated`)) +
  geom_point(size=3,aes(colour=species)) +
  stat_cor(p.accuracy=0.001,r.accuracy=0.01,size=3) +
  scale_colour_manual(values=colours) +
  facet_wrap(~species,ncol=7) +
  geom_abline(slope=1,intercept=0,linetype="dashed") +
  xlab("q = 0") +
  ylab("q = 0") +
  xlim(0,max(toplot$`hill_zero_exuded`,toplot$`hill_zero_host-associated`)) +
  ylim(0,max(toplot$`hill_zero_exuded`,toplot$`hill_zero_host-associated`)) +
  coord_fixed() +
  theme_bw()
plot_one<-ggplot(data=toplot,aes(x=`hill_one_exuded`,y=`hill_one_host-associated`)) +
  geom_point(size=3,aes(colour=species)) +
  stat_cor(p.accuracy=0.001,r.accuracy=0.01,size=3) +
  scale_colour_manual(values=colours) +
  facet_wrap(~species,ncol=7) +
  geom_abline(slope=1,intercept=0,linetype="dashed") +
  xlab("q = 1") +
  ylab("host-associated\nq = 1") +
  xlim(0,max(toplot$`hill_one_exuded`,toplot$`hill_one_host-associated`)) +
  ylim(0,max(toplot$`hill_one_exuded`,toplot$`hill_one_host-associated`)) +
  coord_fixed() +
  theme_bw()
plot_two<-ggplot(data=toplot,aes(x=`hill_two_exuded`,y=`hill_two_host-associated`)) +
  geom_point(size=3,aes(colour=species)) +
  stat_cor(p.accuracy=0.001,r.accuracy=0.01,size=3) +
  scale_colour_manual(values=colours) +
  facet_wrap(~species,ncol=7) +
  geom_abline(slope=1,intercept=0,linetype="dashed") +
  xlab("q = 2\nexuded") +
  ylab("q = 2") +
  xlim(0,max(toplot$`hill_two_exuded`,toplot$`hill_two_host-associated`)) +
  ylim(0,max(toplot$`hill_two_exuded`,toplot$`hill_two_host-associated`)) +
  coord_fixed() +
  theme_bw()

# Combine plots
patchwork::wrap_plots(list(plot_zero,plot_one,plot_two),guides="collect",ncol=1) +
  plot_annotation(title="Correlation of alpha diversity metrices of exuded vs host-associated microbial communities",)

ggsave(filename="../data/20_ADiv_host-exudates.png",width=unit(15,"cm"),height=unit(7,"cm"))
