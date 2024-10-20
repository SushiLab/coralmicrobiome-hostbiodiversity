## Project: Aquarium Microbiome
## Script purpose: Explore how sample-to-control distance changes after decontamination
## Date: September 2023
## Author: Fabienne Wiederkehr

# Load packages
require(data.table)
require(tidyverse)
require(vegan)
require(ggforce)
require(ape)
require(pairwiseAdonis)
require(patchwork)

# Load data
asv_abtab<-fread("../data/asv_abtab_original_free-living.tsv",sep="\t",header=TRUE,data.table=FALSE) %>%
  column_to_rownames("id") %>% as.matrix()
asv_abtab_decontam<-fread("../data/decontam_free-living/asv_dat_decontam_thresh-0.2.tsv",sep="\t",header=TRUE,data.table=FALSE) %>%
  select(asv,where(is.numeric)) %>%
  column_to_rownames("asv") %>%
  t()
metadat<-fread("../data/metadata_free-living.csv",header=TRUE,data.table=FALSE)
colours<-readRDS("../data/colour_palette.RDS")

# Choose rarefaction cutoff
rr_cutoff<-1000

# Rarefy original (= no decontam) abundance table
asv_abtab_nodecontam_rr<-asv_abtab[which(apply(asv_abtab,1,sum)>=rr_cutoff),] #only keep samples >=1000 reads
asv_abtab_nodecontam_rr<-rrarefy(asv_abtab_nodecontam_rr,sample=rr_cutoff) #subsample to 1000 reads
asv_abtab_nodecontam_rr<-asv_abtab_nodecontam_rr[,which(apply(asv_abtab_nodecontam_rr,2,sum)>0)] #remove all-0 (non-present when downsampling) ASVs

# Analyse sample-to-blank distances (original data)
bc_dist_nodecontam<-vegdist(sqrt(asv_abtab_nodecontam_rr)) %>% as.matrix() %>% as.data.frame() %>%
  rownames_to_column("id_1") %>%
  pivot_longer(-id_1,names_to="id_2",values_to="distance") %>%
  left_join(select(metadat,id,sample,species,aDNA_smear_conc_nmol_l),by=c("id_1"="id")) %>%
  left_join(select(metadat,id,sample),by=c("id_2"="id")) %>%
  filter(sample.x %in% c("tank water","exudates","incubation control","microscopy-slide control") & sample.y %in% c("PCR blank","extraction blank")) %>%
  mutate(sample.x=factor(sample.x))

# Plot sample-to-blank distances (original data)
g_nodecontam<-ggplot(data=bc_dist_nodecontam,aes(x=as.numeric(aDNA_smear_conc_nmol_l),y=distance,colour=sample.x)) +
  geom_point(alpha=0.9) +
  scale_x_sqrt() +
  scale_colour_manual(values=colours) +
  ylim(0,1) +
  xlab("DNA concentration of the sample (nmol/l)") +
  ylab("Bray-Curtis distance to a control") +
  labs(title="Original data") +
  guides(colour="none") +
  theme_bw()

# Rarefy cleaned abundance table (threshold = 0.1)
asv_abtab_decontam_rr<-asv_abtab_decontam[which(apply(asv_abtab_decontam,1,sum)>=rr_cutoff),]
asv_abtab_decontam_rr<-rrarefy(asv_abtab_decontam_rr,sample=rr_cutoff)
asv_abtab_decontam_rr<-asv_abtab_decontam_rr[,which(apply(asv_abtab_decontam_rr,2,sum)>0)]

# Analyse sample-to-blank distances (cleaned data)
bc_dist_decontam<-vegdist(sqrt(asv_abtab_decontam_rr)) %>% as.matrix() %>% as.data.frame() %>%
  rownames_to_column("id_1") %>%
  pivot_longer(-id_1,names_to="id_2",values_to="distance") %>%
  left_join(select(metadat,id,sample,species,aDNA_smear_conc_nmol_l),by=c("id_1"="id")) %>%
  left_join(select(metadat,id,sample),by=c("id_2"="id")) %>%
  filter(sample.x %in% c("tank water","exudates","incubation control","microscopy-slide control") & sample.y %in% c("PCR blank","extraction blank")) %>%
  mutate(sample.x=factor(sample.x,levels=c("tank water","exudates","incubation control","microscopy-slide control")))

# Plot sample-to-blank distances (cleaned data)
g_decontam<-ggplot(data=bc_dist_decontam,aes(x=as.numeric(aDNA_smear_conc_nmol_l),y=distance,colour=sample.x)) +
  geom_point(alpha=0.9) +
  scale_x_sqrt() +
  scale_colour_manual(values=colours) +
  ylim(0,1) +
  xlab("DNA concentration of the sample (nmol/l)") +
  ylab("Bray-Curtis distance to a control") +
  labs(title="Contaminants removed",colour="Sample type") +
  theme_bw()

# Plot
g_nodecontam | g_decontam
ggsave(filename="../data/09_sample-to-blank-original-vs-cleaned.png",width=10,height=7,dpi=500)

# Create dataframe with factor decontam/nodecontam
df<-bc_dist_nodecontam %>%
  mutate(status="not decontaminated") %>%
  bind_rows(bc_dist_decontam) %>%
  mutate(status=replace_na(status,"decontaminated")) %>%
  mutate(status=fct_relevel(status,"not decontaminated","decontaminated"))

# Plot BC distance sample-to-control before and after decontamination by species
ggplot(data=df,aes(x=fct_reorder(species,distance,median),y=distance,fill=status)) +
  geom_violin(draw_quantiles=0.5,scale="width",alpha=0.5) +
  geom_sina(scale="width",aes(colour=status)) +
  scale_fill_manual(values=c("#000000","#d9d9d9")) +
  scale_colour_manual(values=c("#000000","#d9d9d9")) +
  ylim(0,1) +
  ylab("Bray-Curtis distance to a control") +
  labs(fill="Decontaminated\nor not?",colour="Decontaminated\nor not?") +
  theme_bw() +
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
ggsave(filename="../data/09_sample-to-blank-original-vs-cleaned-by-species.png",width=10,height=7,dpi=500)

# Rarefy at 1000 reads
rr_cutoff<-1000
asv_abtab_rr<-asv_abtab_decontam[which(apply(asv_abtab_decontam,1,sum)>=rr_cutoff),] #only keep samples >=1000 reads
asv_abtab_rr<-rrarefy(asv_abtab_rr,sample=rr_cutoff) #subsample to 1000 reads
asv_abtab_rr<-asv_abtab_rr[,which(apply(asv_abtab_rr,2,sum)>0)] #remove all-0 (non-present when downsampling) ASVs

# Compare communities (PCoA)
pcoa_res<-pcoa(D=vegdist(sqrt(asv_abtab_rr))) #sqrt-transform to lessen impact of abundant ASVs
pcoa_df<-pcoa_res$vectors %>% as.data.frame() %>%
  rownames_to_column("id") %>%
  left_join(metadat,by="id") %>%
  mutate(species=ifelse(is.na(species),sample,species)) %>%
  mutate(sample=factor(sample,levels=c("exudates","tank water","incubation control",
                                       "extraction mock","PCR mock",
                                       "microscopy-slide control","extraction blank","PCR blank"))) %>%
  mutate(species=factor(species,levels=c("Montipora digitata","Pocillopora verrucosa",
                                         "Sinularia sp.","Xenia sp.",
                                         "Caulerpa sp.","Peyssonnelia sp.",
                                         "Haliclona cnidata","Dictyoceratida",
                                         "tank water","incubation control","microscopy-slide control",
                                         "extraction mock","PCR mock",
                                         "extraction blank","PCR blank")))

# Plot community distances (PCoA) by sample type
ggplot(data=pcoa_df,aes(x=Axis.1,y=Axis.2,colour=species)) +
  geom_hline(yintercept=0,linetype=2) +
  geom_vline(xintercept=0,linetype=2) +
  geom_point(size=4,alpha=0.9) +
  scale_colour_manual(values=colours) +
  xlab(bquote(.(round(pcoa_res$values$Relative_eig[1]*100,2))~"%")) +
  ylab(bquote(.(round(pcoa_res$values$Relative_eig[2]*100,2))~"%")) +
  facet_wrap(~sample) +
  theme_bw()
ggsave(filename="../data/09_PCoA_by-species.png",width=10,height=7,dpi=500)

# PERMANOVA
metadat_perm<-metadat %>%
  filter(!is.na(sequencing_run)) %>%
  filter(sample %in% c("exudates")) %>%
  filter(!grepl("w_P0_MdiPcPcMo2",id)) #this sample dropped out after sequencing

asv_abtab_perm<-asv_abtab_decontam %>% as.data.frame() %>%
  rownames_to_column("id") %>%
  filter(id %in% metadat_perm$id) %>%
  column_to_rownames("id")

metadat_perm<-data.frame(id=rownames(asv_abtab_perm)) %>% #make sure to keep order of samples
  left_join(metadat_perm,by="id")

adonis2(vegdist(sqrt(asv_abtab_perm))~species,data=metadat_perm,permutations=1000,by="terms")

# PERMDISP
permutest(betadisper(vegdist(sqrt(asv_abtab_perm)),metadat_perm$species),permutations=how(nperm=1000))

# Perform post-hoc test for dissimilarity (ADONIS)
posthoc_adonis<-pairwise.adonis2(vegdist(sqrt(asv_abtab_perm))~species,data=metadat_perm,nperm=1000,p.adjust.method="fdr")

# Explore communities by taxonomic group
comm_tax_group<-function(x) {
  metadat_tax_group<-metadat %>%
    filter(grepl(x,taxonomic_group))
  
  asv_abtab_tax_group<-asv_abtab_rr %>% as.data.frame() %>%
    filter(rownames(.) %in% metadat_tax_group$id)
  
  pcoa_tax_group<-pcoa(D=vegdist(sqrt(asv_abtab_tax_group)))
  pcoa_df_tax_group<-pcoa_tax_group$vectors %>% as.data.frame() %>%
    rownames_to_column("id") %>%
    left_join(metadat,by="id") %>%
    mutate(species=ifelse(is.na(species),sample,species)) %>%
    mutate(species=factor(species,levels=c("Montipora digitata","Pocillopora verrucosa",
                                           "Sinularia sp.","Xenia sp.",
                                           "Caulerpa sp.","Peyssonnelia sp.",
                                           "Haliclona cnidata","Dictyoceratida",
                                           "tank water","incubation control","microscopy-slide control",
                                           "extraction mock","PCR mock",
                                           "extraction blank","PCR blank")))
  
  ggplot(data=pcoa_df_tax_group,aes(x=Axis.1,y=Axis.2,colour=species,shape=treatment)) +
    geom_hline(yintercept=0,linetype=2) +
    geom_vline(xintercept=0,linetype=2) +
    geom_point(size=4,alpha=0.9) +
    scale_colour_manual(values=colours) +
    scale_shape_manual(values=c(15,16,17)) +
    labs(colour=x) +
    xlab(bquote(.(round(pcoa_tax_group$values$Relative_eig[1]*100,2))~"%")) +
    ylab(bquote(.(round(pcoa_tax_group$values$Relative_eig[2]*100,2))~"%")) +
    facet_wrap(~sample) +
    theme_bw()
}

plots<-list(comm_tax_group("sponge"),comm_tax_group("coral"),comm_tax_group("macroalgae"))
patchwork::wrap_plots(plots,guides="collect",ncol=3)
ggsave(filename="../data/09_PCoA_by-organismal-group.png",width=unit(12,"cm"),height=unit(5,"cm"))
