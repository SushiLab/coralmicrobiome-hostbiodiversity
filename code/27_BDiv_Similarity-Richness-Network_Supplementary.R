## Project: Aquarium Microbiome
## Script purpose: Plot similarity between and richness of hosts across treatment and source
## Date: May 2024
## Author: Fabienne Wiederkehr

# Load packages
require(data.table)
require(tidyverse)
require(ggplot2)
require(igraph)
require(ggraph)
require(patchwork)

# Load paths
ed_paths<-list.files("../data/",pattern="24_BDiv_ed",full.names=TRUE)
va_paths<-list.files("../data/",pattern="25_BDiv_va",full.names=TRUE)

# Initiate variables
plots<-list(NULL)
counter<-1

for (i in 1:length(va_paths)){
  # Load files
  va<-fread(va_paths[i],sep="\t",header=TRUE,data.table=FALSE) %>%
    arrange(match(species,c("Xenia sp.","Peyssonnelia sp.","Pocillopora verrucosa","Sinularia sp.","Montipora digitata","Haliclona cnidata","Caulerpa sp.")))
  ed<-fread(ed_paths[i],sep="\t",header=TRUE,data.table=FALSE)
  
  # Initiate plot
  ig<-igraph::graph_from_data_frame(d=ed,vertices=va,directed=FALSE)
  tg<-tidygraph::as_tbl_graph(ig) %>% 
    tidygraph::activate(nodes) %>% 
    dplyr::mutate(label=name)
  richness<-V(tg)$richness
  
  # Edge size shows distance between hosts
  g<-tg %>%
    ggraph(layout="auto") +
    geom_edge_arc(lineend="round",
                  strength=0.1,
                  aes(edge_width=similarity,edge_colour=similarity)) +
    geom_node_point(size=richness/100,
                    colour=c("#cbc9e2","#bae4b3","#fcae91","#9e9ac8","#fb6a4a","#6baed6","#74c476")) +
    geom_node_text(aes(label=name), 
                   repel=TRUE,
                   check_overlap=TRUE, 
                   colour="#000000") +
    scale_edge_width(breaks=c(0.01,0.05,0.1,0.15,0.2,0.25,0.3),
                     limits=c(0.01,0.3),
                     range=c(0,7)) +
    scale_edge_colour_gradient(breaks=c(0.01,0.05,0.1,0.15,0.2,0.25,0.3),
                               limits=c(0.01,0.3),
                               low="#d9d9d9",
                               high="#000000") +
    labs(title=gsub(".*_(.*?)\\.tsv","\\1",basename(va_paths[i]))) +
    theme_graph(base_family="sans",background="#FFFFFF") +
    theme(legend.position="top",
          plot.title=element_text(size=11))
  
  plots[[counter]]<-g
  counter<-counter+1
}

# Create patchwork plot
patchwork::wrap_plots(plots,ncol=2,guides="collect")
ggsave(filename="../data/27_BDiv_similarity-richness-network.png",width=unit(15,"cm"),height=unit(12,"cm"))
