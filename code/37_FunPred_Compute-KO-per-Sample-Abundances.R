## Project: Aquarium Microbiome
## Script purpose: From K ids, compute pathway abundances averaged over the 50 picrust output files (~rarefaction files)
## Date: July 2024
## Author: Fabienne Wiederkehr

# Load packages
require(data.table)
require(tidyverse)
require(vegan)

# Load paths
path_abtab_paths<-list.files("../data/picrust_output/",pattern="asv_abtab_rarefied_",full.names=TRUE)

# Get a list of all functions
path_list<-fread(paste(path_abtab_paths[1],"/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz",sep=""),sep="\t",header=TRUE,data.table=FALSE) %>%
  select(`function`)
for (i in 2:length(path_abtab_paths)){
  # Load data and extract function
  path_abtab<-fread(paste(path_abtab_paths[i],"/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz",sep=""),sep="\t",header=TRUE,data.table=FALSE) %>%
    select(`function`)
  
  # Save functions in list
  path_list<-intersect(path_abtab,path_list)
}

# Compute function abundances
path_abtab_list<-list(NULL)
for (i in 1:length(path_abtab_paths)){
  # Load data
  path_abtab<-fread(paste(path_abtab_paths[i],"/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz",sep=""),sep="\t",header=TRUE,data.table=FALSE)
  
  # Ensure functions are the same in each file
  path_dat<-path_list %>%
    left_join(path_abtab,by="function") %>%
    column_to_rownames("function")
  
  # Save function in list
  path_abtab_list[[i]]<-path_dat
}

# Function to check if all column names are identical
check_column_names<-function(dfs) {
  col_names<-lapply(dfs,colnames)
  all_identical<-all(sapply(col_names,function(x) identical(x,col_names[[1]])))
  return(all_identical)
}

# Check if column names are identical
check_column_names(path_abtab_list)

# Reduce list to a single dataframe
path_abtab_combined<-Reduce("+",path_abtab_list) #sum up the corresponding elements
path_abtab_mean<-path_abtab_combined/length(path_abtab_list) #calculate mean of each element

# Write file
path_abtab_mean<-path_abtab_mean %>%
  rownames_to_column("ko_id")
fwrite(path_abtab_mean,"../data/ko_abtab_sample.tsv",sep="\t",na="NA")
