## Project: Aquarium Microbiome
## Script purpose: Generate a maximum likelihood tree for the significantly differentially abundant ASVs
## Date: April 2024
## Author: Fabienne Wiederkehr

# Load packages
require(data.table)
require(tidyverse)
require(msa)
require(phangorn)

# Load data
asv_da<-fread("../data/asv_diffab_lfc.tsv",sep="\t",header=TRUE,data.table=FALSE) %>%
  select(asv) %>%
  left_join(fread("../data/asv_diffab_taxinfo.tsv",sep="\t",header=TRUE,data.table=FALSE),by="asv")

# Align sequences through a multiple sequence aligner
seq<-DNAStringSet(asv_da$seq)
names(seq)<-asv_da$asv #names sequences
msa<-msa(seq,method="Muscle")

# Convert to phyDat object
phyDat<-msaConvert(msa,type="phangorn::phyDat")

# Select model
mt<-modelTest(phyDat,control=pml.control(trace=0))

# Compute the two best models
fit_TVM<-pml_bb(mt,model="TVM+G(4)+I")
fit_GTR<-pml_bb(mt,model="GTR+G(4)+I")

# Bootstrap
bs<-bootstrap.pml(fit_TVM,bs=100,optNni=TRUE,control=pml.control(trace=0))
bs<-bootstrap.pml(fit_GTR,bs=100,optNni=TRUE,control=pml.control(trace=0))

# Assign standard bootstrap values to the tree
tree_TVM<-plotBS(fit_TVM$tree,bs,type="n")
tree_GTR<-plotBS(fit_GTR$tree,bs,type="n")

# Export tree
write.tree(tree_TVM,"../data/74_ASVs_DA_TVM.tree")
write.tree(tree_GTR,"../data/74_ASVs_DA_GTR.tree")
