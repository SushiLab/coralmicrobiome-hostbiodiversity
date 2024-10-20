## Project: Aquarium Microbiome
## Script purpose: Combine metadata and explore sequencing depth
## Date: October 2023
## Author: Fabienne Wiederkehr

# Load packages
require(data.table)
require(tidyverse)
require(ggforce)

# Load data
metadat_free<-fread("../data/metadata_free-living.csv",header=TRUE,data.table=FALSE) %>%
  filter(!is.na(sequencing_run)) %>%
  filter(!grepl("w_P0_MdiPcPcMo2",id)) %>% #this sample dropped out after sequencing
  filter(sample %in% c("exudates","tank water","incubation control")) %>%
  mutate(source=ifelse(sample=="exudates","exuded","free-living")) %>%
  select(id,source,sample,type,fragment_id,species,taxonomic_group,phase,treatment,treatment_binary,replicate,tank,sampling_date)
metadat_host<-fread("../data/metadata_host-associated.csv",header=TRUE,data.table=FALSE) %>%
  filter(!is.na(sequencing_run)) %>%
  filter(sample %in% c("mucus","biofilm","tissue")) %>%
  mutate(source="host-associated") %>%
  select(id,source,sample,type,fragment_id,species,taxonomic_group,phase,treatment,treatment_binary,replicate,tank,sampling_date)
metadat<-metadat_free %>%
  bind_rows(metadat_host)
fwrite(metadat,"../data/metadat_cleaned.tsv",sep="\t",na="NA")

asv_abtab<-fread("../data/asv_abtab_cleaned.tsv",sep="\t",header=TRUE,data.table=FALSE) %>%
  column_to_rownames("id")

colours<-readRDS("../data/colour_palette.RDS")

# Compute sequencing depth
seq_depth_df<-apply(asv_abtab,1,sum) %>% as.data.frame() %>% #sum all ASVs per sample
  rename("seq_depth"=".") %>%
  rownames_to_column("id") %>%
  left_join(metadat,by="id") %>% #append calculated sequencing depth to metadat
  mutate(species=ifelse(is.na(species),sample,species)) %>% #add control type to column species
  mutate(species=factor(species,levels=c("Montipora digitata","Pocillopora verrucosa",
                                         "Sinularia sp.","Xenia sp.",
                                         "Caulerpa sp.","Peyssonnelia sp.",
                                         "Haliclona cnidata","Dictyoceratida",
                                         "tank water","incubation control"))) %>%
  mutate(source=factor(source,levels=c("host-associated","exuded","free-living")))

# Plot sequencing depth
ggplot(data=seq_depth_df,aes(x=fct_reorder(species,seq_depth,median),y=seq_depth)) +
  geom_violin(draw_quantiles=0.5,scale="width",alpha=0.5,aes(fill=species,colour=species)) +
  geom_sina(scale="width",aes(colour=source)) +
  scale_fill_manual(values=colours) +
  scale_colour_manual(values=colours) +
  scale_y_log10() +
  ylab("Sequencing depth") +
  theme_bw() +
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="none")
ggsave(filename="../data/12_seq-depth.png",width=10,height=7,dpi=500)
