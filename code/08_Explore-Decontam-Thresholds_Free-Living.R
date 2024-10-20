## Project: Aquarium Microbiome
## Script purpose: Run decontam with various thresholds to inform decision
## Date: September 2023
## Author: Fabienne Wiederkehr

# Load packages
require(data.table)
require(tidyverse)
require(TAF)
require(vegan)
require(decontam)

# Load data
asv_abtab<-fread("../data/asv_abtab_original_free-living.tsv",sep="\t",header=TRUE,data.table=FALSE) %>%
  column_to_rownames("id") %>% as.matrix()
asv_dat<-fread("../data/asv_dat_original_free-living.tsv",sep="\t",header=TRUE,data.table=FALSE)
metadat<-fread("../data/metadata_free-living.csv",header=TRUE,data.table=FALSE) %>%
  filter(!is.na(sequencing_run))
colours<-readRDS("../data/colour_palette.RDS")

# Prepare list of negatives
is_negative<-data.frame(id=rownames(asv_abtab)) %>%
  left_join(select(metadat,id,sample,aDNA_smear_conc_nmol_l),by="id") %>%
  mutate(is_neg=ifelse(sample %in% c("PCR blank","extraction blank"),TRUE,FALSE)) %>%
  pull(is_neg)

# Prepare concentrations
conc<-data.frame(id=rownames(asv_abtab)) %>%
  left_join(select(metadat,id,sample,aDNA_smear_conc_nmol_l),by="id") %>%
  mutate(dna_conc=as.numeric(aDNA_smear_conc_nmol_l)) %>%
  mutate(dna_conc=ifelse(dna_conc<=0,0.0001,dna_conc)) %>%
  pull(dna_conc)

# Run decontam (prevalence method) with different thresholds and generate decontaminated ASV files
mkdir("../data/decontam_free-living")
for (thresh in c(0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5)){
  cat("Running decontam with P threshold of",thresh,"...")
  res<-isContaminant(asv_abtab,conc=conc,neg=is_negative,threshold=thresh) %>%
    rownames_to_column(var="asv")
  asv_dat_decontam<-asv_dat
  asv_dat_decontam[match(res$asv[res$contaminant==TRUE],asv_dat_decontam$asv),3:ncol(asv_dat_decontam)]<-0
  fwrite(res,file=paste("../data/decontam_free-living/decontam_stats_thresh-",thresh,".tsv",sep=""),sep="\t")
  fwrite(asv_dat_decontam,file=paste("../data/decontam_free-living/asv_dat_decontam_thresh-",thresh,".tsv",sep=""),sep="\t")
  cat("\n")
}

# Explore percentage of original reads/ASVs kept after decontam using various thresholds
# Load decontaminated ASV files
asv_dat_paths<-list.files("../data/decontam_free-living/",pattern="asv_dat_decontam_thresh",full.names=TRUE)

# Calculate number of reads/ASVs kept for each of the decontam thresholds (plus full = no decontam) and bind to dataframe
res<-NULL
for (i in asv_dat_paths){
  tmp<-fread(i,sep="\t",header=TRUE,data.table=FALSE)
  tmp<-data.frame(nreads=colSums(tmp[,-c(1:2)]),
                  nspec=specnumber(t(tmp[,-c(1:2)])),
                  file=basename(i)) %>%
    rownames_to_column("sample")
  res<-res %>%
    bind_rows(tmp)
}
tmp<-data.frame(nreads=colSums(asv_dat[,-c(1:2)]),
                nspec=specnumber(t(asv_dat[,-c(1:2)])),
                file="full") %>%
  rownames_to_column("sample")
res<-res %>%
  bind_rows(tmp)

# Compute percentage of original reads kept for each of the decontam thresholds (plus full = no decontam)
res_nreads<-res %>%
  select(-nspec) %>%
  pivot_wider(names_from="file",values_from="nreads") %>%
  mutate(across(is.numeric,~./full)) %>% #calculate percentage of original reads kept (reads after decontam divided by reads before decontam)
  pivot_longer(-sample,names_to="file",values_to="perc_original_reads")

# Compute percentage of original ASVs kept for each of the decontam thresholds (plus full = no decontam)
res_nspec<-res %>%
  select(-nreads) %>%
  pivot_wider(names_from="file",values_from="nspec") %>%
  mutate(across(is.numeric,~./full)) %>% #calculate percentage of original ASVs kept (ASVs after decontam divided by ASVs before decontam)
  pivot_longer(-sample,names_to="file",values_to="perc_original_spp")

# Prepare dataframe
df_plot_species<-res_nreads %>%
  left_join(res_nspec,by=c("file","sample")) %>%
  left_join(metadat,by=c("sample"="id")) %>%
  mutate(species=ifelse(is.na(species),sample.y,species)) %>% #add control type to column species
  group_by(file,species, sample.y) %>%
  summarise(perc_original_reads_median=median(perc_original_reads), #compute median of percentage of original reads kept
            perc_original_reads_25CI=quantile(perc_original_reads,probs=0.25),
            perc_original_reads_75CI=quantile(perc_original_reads,probs=0.75),
            perc_original_spp_median=median(perc_original_spp), #compute median of percentage of original ASVs kept
            perc_original_spp_25CI=quantile(perc_original_spp,probs=0.25),
            perc_original_spp_75CI=quantile(perc_original_spp,probs=0.75)) %>%
  mutate(file=gsub("asv_dat_decontam_thresh-","",file)) %>% #keep only threshold value in column file
  mutate(file=gsub(".tsv","",file)) %>% #keep only threshold value in column file
  mutate(file=ifelse(file=="full",0,as.numeric(file))) %>% #replace "full" with 0 in column file
  mutate(species=factor(species,levels=c("Montipora digitata","Pocillopora verrucosa",
                                         "Sinularia sp.","Xenia sp.",
                                         "Caulerpa sp.","Peyssonnelia sp.",
                                         "Haliclona cnidata","Dictyoceratida",
                                         "tank water","incubation control","microscopy-slide control",
                                         "extraction blank","PCR blank",
                                         "extraction mock","PCR mock")))
  
# Plot percentage of original reads kept by species
ggplot(data=df_plot_species,aes(x=file,y=perc_original_reads_median,ymin=perc_original_reads_25CI,ymax=perc_original_reads_75CI,fill=species,colour=species)) +
  geom_ribbon(alpha=0.6) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept=0.1,linetype="dashed") +
  scale_fill_manual(values=colours) +
  scale_colour_manual(values=colours) +
  facet_wrap(~species) +
  xlab("Decontam threshold") +
  ylab("Percentage of original reads kept\n(median, 25th and 75th percentile)") +
  guides(colour="none",fill="none") +
  theme_bw()
ggsave(filename="../data/08_perc-reads-kept-by-species.png",width=10,height=7,dpi=500)

# Plot percentage of original ASVs kept by species
ggplot(data=df_plot_species,aes(x=file,y=perc_original_spp_median,ymin=perc_original_spp_25CI,ymax=perc_original_spp_75CI,fill=species,colour=species)) +
  geom_ribbon(alpha=0.6) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept=0.1,linetype="dashed") +
  scale_fill_manual(values=colours) +
  scale_colour_manual(values=colours) +
  facet_wrap(~species) +
  xlab("Decontam threshold") +
  ylab("Percentage of original ASVs kept\n(median, 25th and 75th percentile)") +
  guides(colour="none",fill="none") +
  theme_bw()
ggsave(filename="../data/08_perc-asvs-kept-by-species.png",width=10,height=7,dpi=500)
