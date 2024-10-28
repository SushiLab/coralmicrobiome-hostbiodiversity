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

for (i in c("exuded","host-associated")){
  # Load data
  ed<-fread(paste("../data/24_BDiv_ed_",i,"biodiverse.tsv",sep=""),sep="\t",header=TRUE,data.table=FALSE) %>%
    rename("biodiverse"="similarity") %>%
    left_join(fread(paste("../data/24_BDiv_ed_",i,"degraded.tsv",sep=""),sep="\t",header=TRUE,data.table=FALSE)) %>%
    rename("degraded"="similarity") %>%
    mutate(similarity=0.5-degraded/(biodiverse+degraded)) %>% #compute ratio
    filter(!is.na(similarity))
  
  va<-fread(paste("../data/25_BDiv_va_",i,"biodiverse.tsv",sep=""),sep="\t",header=TRUE,data.table=FALSE) %>%
    rename("biodiverse"="richness") %>%
    left_join(fread(paste("../data/25_BDiv_va_",i,"degraded.tsv",sep=""),sep="\t",header=TRUE,data.table=FALSE)) %>%
    rename("degraded"="richness") %>%
    mutate(richness=0.5-degraded/(biodiverse+degraded)) %>%
    arrange(match(species,c("Xenia sp.","Peyssonnelia sp.","Pocillopora verrucosa","Sinularia sp.","Montipora digitata","Haliclona cnidata","Caulerpa sp.")))
  
  # Initiate plot
  ig<-igraph::graph_from_data_frame(d=ed,vertices=va,directed=FALSE)
  tg<-tidygraph::as_tbl_graph(ig) %>% 
    tidygraph::activate(nodes) %>% 
    dplyr::mutate(label=name)
  richness<-V(tg)$richness
  
  # Edge and node plot
  tg %>%
    ggraph(layout="auto") +
    geom_edge_arc(lineend="round",
                  strength=0.1,
                  aes(edge_width=abs(similarity),edge_colour=similarity>0)) + #check if similarity >0
    geom_node_point(size=abs(richness*25),
                    colour=ifelse(richness>0,"#D0ABA4","#6B6969")) +
    geom_node_text(aes(label=name), 
                   repel=TRUE,
                   check_overlap=TRUE, 
                   colour="#000000") +
    scale_edge_width(breaks=c(0,0.1,0.2,0.3),
                     limits=c(0,0.35),
                     range=c(0,4)) +
    scale_edge_colour_manual(values=c("TRUE"="#D0ABA4","FALSE"="#6B6969"), 
                             breaks=c(TRUE,FALSE),
                             guide=NULL) +
    theme_graph(base_family="sans",background="#FFFFFF") +
    theme(aspect.ratio=1,
          legend.position="top",legend.title=element_blank())
  ggsave(filename=paste("../data/26_BDiv_similarity-richness-network_",i,".png",sep=""),width=15,height=20,units="cm")
}
