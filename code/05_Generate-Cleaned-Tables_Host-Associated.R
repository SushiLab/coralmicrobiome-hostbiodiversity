## Project: Aquarium Microbiome
## Script purpose: Generate cleaned tables for subsequent analyses
## Date: May 2023
## Author: Fabienne Wiederkehr

# Load packages
require(data.table)
require(tidyverse)
require(ggforce)

# Load data
metadat<-fread("../data/metadata_host-associated.csv",header=TRUE,data.table=FALSE) %>%
  filter(!is.na(sequencing_run)) %>%
  filter(sample %in% c("mucus","biofilm","tissue"))
fwrite(metadat,"../data/metadat_cleaned_host-associated.tsv",sep="\t",na="NA")

asv_dat<-fread("../data/decontam_host-associated/asv_dat_decontam_thresh-0.1.tsv",sep="\t",header=TRUE,data.table=FALSE) %>%
  select(asv,seq,all_of(metadat$id))
fwrite(asv_dat,"../data//asv_dat_cleaned_host-associated.tsv",sep="\t",na="NA")

colours<-readRDS("../data/colour_palette.RDS")

# Generate ASV abundance table
asv_abtab<-asv_dat %>%
  select(asv,where(is.numeric)) %>%
  column_to_rownames("asv") %>%
  t()

# Compute sequencing depth
seq_depth_df<-apply(asv_abtab,1,sum) %>% as.data.frame() %>% #sum all ASVs per sample
  rename("seq_depth"=".") %>%
  rownames_to_column("id") %>%
  left_join(metadat,by="id") %>% #append calculated sequencing depth to metadat
  mutate(species=ifelse(is.na(species),sample.y,species)) %>% #add control type to column species
  mutate(species=factor(species,levels=c("Montipora digitata","Pocillopora verrucosa",
                                         "Sinularia sp.","Xenia sp.",
                                         "Caulerpa sp.","Peyssonnelia sp.",
                                         "Haliclona cnidata","Dictyoceratida",
                                         "extraction mock","PCR mock",
                                         "extraction blank","PCR blank","PBS control","Tris-NaCl control")))

# Plot sequencing depth
ggplot(data=seq_depth_df,aes(x=fct_reorder(species,seq_depth,median),y=seq_depth,fill=species,colour=species)) +
  geom_violin(draw_quantiles=0.5,scale="width",alpha=0.5) +
  geom_sina(scale="width") +
  scale_fill_manual(values=colours) +
  scale_colour_manual(values=colours) +
  scale_y_log10() +
  ylab("Sequencing depth") +
  theme_bw() +
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="none")
ggsave(filename="../data/05_seq-depth.png",width=10,height=7,dpi=500)

# Write ASV abundance table
asv_abtab<-asv_abtab %>% as.data.frame() %>%
  rownames_to_column("id")
fwrite(asv_abtab,"../data/asv_abtab_cleaned_host-associated.tsv",sep="\t",na="NA")
