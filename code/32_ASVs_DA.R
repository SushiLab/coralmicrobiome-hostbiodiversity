## Project: Aquarium Microbiome
## Script purpose: Analyse the differentially abundant ASVs between host-associated and free-living microbial communities
## Date: February 2024
## Author: Fabienne Wiederkehr

# Load packages
require(data.table)
require(tidyverse)
require(TAF)
require(TreeSummarizedExperiment)
require(foreach)
require(doRNG)
require(ape)
require(patchwork)

# Load function to source .R files from local directory
sourceDir<-function(path,trace=TRUE,...) {
  op<-options();on.exit(options(op)) # to reset after each 
  for (nm in list.files(path, pattern="[.][R]$")) {
    if(trace) cat(nm,":")
    source(file.path(path,nm),...)
    if(trace) cat("\n")
    options(op)
  }
}

# Load functions (downloaded from https://github.com/FrederickHuangLin/ANCOMBC/tree/bugfix/R to subfolder /code/ancombc/R)
sourceDir("ancombc/R/")
#to automate the removal taxa with zero variance errors, I changed the following scripts and code lines:
#ancombc2.R lines 532-538:
# if (nrow(bias1) == length(fix_eff)){ #added by FW
#   bias1 = data.frame(bias1, row.names = fix_eff, check.names = FALSE)
#   colnames(bias1) = c("delta_em", "delta_wls", "var_delta")
# } #added by FW
# else{ #added by FW
#   return() #added by FW
# } #added by FW

#ancombc_bias_correct.R lines 534-537:
# redflag_taxids = as.list(names(which(nu0 == 0))) #added by FW
# fwrite(redflag_taxids, file = paste("../data/ancombc/redflag_taxids_", str_replace(i, " ", "-"), "_", k, ".tsv", sep = ""), sep = "\t", na = "NA") #added by FW
# return() #added by FW
# # stop(stop_txt, call. = FALSE) #disabled by FW

# Load data
asv_abtab<-fread("../data/asv_abtab_cleaned.tsv",sep="\t",header=TRUE,data.table=FALSE) #no rarefaction needed
metadat<-fread("../data/metadat_cleaned.tsv",sep="\t",header=TRUE,data.table=FALSE) %>%
  filter(phase=="P2") %>%
  filter(!sample %in% c("incubation control","tank water")) %>%
  mutate(treatment_binary=factor(treatment_binary,levels=c("biodiverse","degraded"))) %>%
  mutate(source=factor(source,levels=c("host-associated","exuded"))) %>%
  mutate(species=factor(species,levels=c("Montipora digitata","Pocillopora verrucosa",
                                         "Sinularia sp.","Xenia sp.",
                                         "Caulerpa sp.","Peyssonnelia sp.",
                                         "Haliclona cnidata")))

# Analyse by species
spp<-levels(metadat$species)
source<-levels(metadat$source)
mkdir("../data/ancombc")
mkdir("../data/differential-abundance")

for (k in source){
  for (i in spp){
    cat("Differential abundance analysis of the",k,"microbial community of",i,"\n")
    # Subset and format metadata
    metadat_red<-metadat %>%
      filter(species==i) %>%
      filter(source==k) %>%
      column_to_rownames("id") %>%
      select(source,treatment_binary)
    
    # Subset and transform abundance table
    asv_abtab_red<-metadat_red %>%
      rownames_to_column("id") %>%
      dplyr::select(id) %>%
      left_join(asv_abtab,by="id") %>% #ensure order of samples is the same as in metadat
      column_to_rownames("id") %>% t() #features need to be rows, samples columns
    
    # Construct TreeSummarizedExperiment
    tse<-TreeSummarizedExperiment(assays=S4Vectors::SimpleList(counts=asv_abtab_red),colData=metadat_red)
    
    # Run ancombc2
    ancombc_res<-ancombc2(data=tse,assay_name="counts",fix_formula="treatment_binary",p_adj_method="holm",pseudo_sens=TRUE,prv_cut=0.25,group="treatment_binary",struc_zero=TRUE,verbose=TRUE)
    
    # If zero variances occur, remove taxa and rerun ancombc2()
    if (file.exists(paste("../data/ancombc/redflag_taxids_",str_replace(i," ","-"),"_",k,".tsv",sep=""))){
      # Load redflag taxids (zero variances)
      redflag_taxids<-fread(file=paste("../data/ancombc/redflag_taxids_",str_replace(i," ","-"),"_",k,".tsv",sep=""),sep="\t",header=FALSE,data.table=FALSE)
      
      # Remove taxa for which zero variances where found
      asv_abtab_red_var<-asv_abtab_red %>% as.data.frame() %>%
        filter(!rownames(.) %in% redflag_taxids) %>% as.matrix()
      
      # Construct new TreeSummarizedExperiment
      tse_var<-TreeSummarizedExperiment(assays=S4Vectors::SimpleList(counts=asv_abtab_red_var),colData=metadat_red)
      
      # Run new ancombc2
      ancombc_res<-ancombc2(data=tse_var,assay_name="counts",fix_formula="treatment_binary",p_adj_method="holm",pseudo_sens=TRUE,prv_cut=0.25,group="treatment_binary",struc_zero=TRUE,verbose=TRUE)
    }
    
    # Save ANCOMBC2 output
    fwrite(ancombc_res$res,file=paste("../data/differential-abundance/diffab_",str_replace(i," ","-"),"_",k,".tsv",sep=""),sep="\t",na="NA")
  }
}
