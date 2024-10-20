## Project: Aquarium Microbiome
## Script purpose: Explore data, decide whether merging before decontam is ok (do all neg controls cluster together?)
## Date: August 2023
## Author: Fabienne Wiederkehr

# Load packages
require(data.table)
require(tidyverse)
require(vegan)
require(ggforce)
require(ape)

# Load data
asv_abtab<-fread("../data/asv_abtab_original_free-living.tsv",sep="\t",header=TRUE,data.table=FALSE) %>%
  column_to_rownames("id") %>% as.matrix()
metadat<-fread("../data/metadata_free-living.csv",header=TRUE,data.table=FALSE) %>%
  mutate(sample=factor(sample,levels=c("exudates","tank water",
                                       "incubation control","microscopy-slide control",
                                       "extraction mock","PCR mock",
                                       "extraction blank","PCR blank")))
colours<-readRDS("../data/colour_palette.RDS")

# Compute sequencing depth
seq_depth_df<-apply(asv_abtab,1,sum) %>% as.data.frame() %>% #sum all ASVs per sample
  rename("seq_depth"=".") %>%
  rownames_to_column("id") %>%
  left_join(metadat,by="id") #append calculated sequencing depth to metadat

# Plot sequencing depth
ggplot(data=seq_depth_df,aes(x=fct_reorder(sample,seq_depth,median),y=seq_depth,fill=sample,colour=sample)) +
  geom_violin(draw_quantiles=0.5,scale="width",alpha=0.5,show.legend=FALSE) +
  scale_fill_manual(values=colours) +
  scale_colour_manual(values=colours) +
  geom_sina(scale="width") +
  scale_y_log10() +
  ylab("Sequencing depth") +
  theme_bw() +
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="none")
ggsave(filename="../data/07_seq-depth.png",width=10,height=7,dpi=500)

# Rarefy at 1000 reads
rr_cutoff<-1000
asv_abtab_rr<-asv_abtab[which(apply(asv_abtab,1,sum)>=rr_cutoff),] #only keep samples >=1000 reads
asv_abtab_rr<-rrarefy(asv_abtab_rr,sample=rr_cutoff) #subsample to 1000 reads
asv_abtab_rr<-asv_abtab_rr[,which(apply(asv_abtab_rr,2,sum)>0)] #remove all-0 (non-present when downsampling) ASVs

# Compare communities (PCoA)
pcoa_res<-pcoa(D=vegdist(sqrt(asv_abtab_rr))) #sqrt-transform to lessen impact of abundant ASVs
pcoa_df<-pcoa_res$vectors %>% as.data.frame() %>%
  rownames_to_column("id") %>%
  left_join(metadat,by="id")

# Plot community distances (PCoA)
ggplot(data=pcoa_df,aes(x=Axis.1,y=Axis.2,colour=sample)) +
  geom_hline(yintercept=0,linetype=2) +
  geom_vline(xintercept=0,linetype=2) +
  geom_point(size=4,alpha=0.9) +
  scale_colour_manual(values=colours) +
  xlab(bquote(.(round(pcoa_res$values$Relative_eig[1]*100,2))~"%")) +
  ylab(bquote(.(round(pcoa_res$values$Relative_eig[2]*100,2))~"%")) +
  theme_bw()

# Plot community distances (PCoA) by sample type
ggplot(data=pcoa_df,aes(x=Axis.1,y=Axis.2,colour=sample)) +
  geom_hline(yintercept=0,linetype=2) +
  geom_vline(xintercept=0,linetype=2) +
  geom_point(size=4,alpha=0.9) +
  scale_colour_manual(values=colours) +
  xlab(bquote(.(round(pcoa_res$values$Relative_eig[1]*100,2))~"%")) +
  ylab(bquote(.(round(pcoa_res$values$Relative_eig[2]*100,2))~"%")) +
  facet_wrap(~sample) +
  guides(colour="none") +
  theme_bw()
ggsave(filename="../data/07_PCoA_by-sample-type.png",width=10,height=7,dpi=500)

# Plot community distances (PCoA) by species
ggplot(data=pcoa_df,aes(x=Axis.1,y=Axis.2,colour=sample)) +
  geom_hline(yintercept=0,linetype=2) +
  geom_vline(xintercept=0,linetype=2) +
  geom_point(size=4,alpha=0.9) +
  scale_colour_manual(values=colours) +
  xlab(bquote(.(round(pcoa_res$values$Relative_eig[1]*100,2))~"%")) +
  ylab(bquote(.(round(pcoa_res$values$Relative_eig[2]*100,2))~"%")) +
  facet_wrap(~species) +
  theme_bw()
ggsave(filename="../data/07_PCoA_by-species.png",width=10,height=7,dpi=500)
