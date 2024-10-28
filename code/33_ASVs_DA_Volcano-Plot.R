## Project: Aquarium Microbiome
## Script purpose: Plot the DA output
## Date: March 2024
## Author: Fabienne Wiederkehr

# Load packages
require(data.table)
require(tidyverse)
require(patchwork)

# Initiate variables
asv_da<-NULL

# Merge output from DA analysis
for (i in c("exuded","host-associated")){
  plots<-list(NULL)
  counter<-1
  for (k in c("Montipora-digitata","Pocillopora-verrucosa",
              "Sinularia-sp.","Xenia-sp.",
              "Caulerpa-sp.","Peyssonnelia-sp.",
              "Haliclona-cnidata")){
    # Load data
    ancombc_res<-fread(paste("../data/differential-abundance/diffab_",k,"_",i,".tsv",sep=""),sep="\t",header=TRUE,data.table=FALSE) %>%
      mutate(species=str_replace(k,"-"," ")) %>%
      mutate(source=i) %>%
      select(asv=taxon,species,source,lfc_treatment_binarydegraded,q_treatment_binarydegraded,diff_treatment_binarydegraded)
    
    # Append DA ASVs to dataframe
    asv_da<-asv_da %>%
      bind_rows(ancombc_res)
    
    # Plot log-fold change of treatment effect
    g<-ggplot(ancombc_res,aes(x=lfc_treatment_binarydegraded*-1,y=q_treatment_binarydegraded,col=diff_treatment_binarydegraded)) + #median difference in clr values
      geom_point(alpha=0.7,show.legend=FALSE) +
      geom_hline(yintercept=0.05,lty=2,col="#C7C7C7") +
      scale_y_log10() +
      scale_colour_manual(values=c("TRUE"="#cb181d","FALSE"="#C7C7C7")) +
      labs(subtitle=str_replace(k,"-"," ")) +
      xlab("degraded vs biodiverse") +
      ylab(bquote("Holm corrected "~italic(p)~"value")) +
      theme_bw()
    
    plots[[counter]]<-g
    counter<-counter+1
  }
  
  # Wrap plots per source
  patchwork::wrap_plots(plots,guides="collect",ncol=2) +
    plot_annotation(title=i)
  
  ggsave(filename=paste("../data/71_ASVs_DA_volcano-plot_",i,".png",sep=""),width=5,height=9,dpi=500)
}

# Extract significant changes
asv_da_sign<-asv_da %>%
  select(-diff_treatment_binarydegraded) %>%
  filter(q_treatment_binarydegraded<=0.05)

# Export merged file
fwrite(asv_da_sign,"../data/asv_diffab.tsv",sep="\t",na="NA")
