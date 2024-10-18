## Project: Aquarium Microbiome
## Script purpose: Colour palette
## Author: Fabienne Wiederkehr
## Date: 14/07/2023

# Created from both https://colorbrewer2.org/#type=sequential&scheme=Blues&n=4 and https://community.tableau.com/thread/235115

colour_palette<-c(
  # Species
  "Montipora digitata"       = "#fb6a4a",
  "Pocillopora verrucosa"    = "#fcae91",
  "Sinularia sp."            = "#9e9ac8",
  "Xenia sp."                = "#cbc9e2",
  "Caulerpa sp."             = "#74c476",
  "Peyssonnelia sp."         = "#bae4b3",
  "Haliclona cnidata"        = "#6baed6",
  
  # Sample
  "mucus"                    = "#ae017e",
  "biofilm"                  = "#238b45",
  "tissue"                   = "#2171b5",
  "tank water"               = "#253494",
  "exudates"                 = "#41b6c4",
  "incubation control"       = "#edf8b1",
  
  # Reef complexity
  "biodiverse"               = "#D0ABA4",
  "degraded"                 = "#6B6969",
  
  # Sample source
  "host-associated"          = "#FFBB78",
  "exuded"                   = "#41b6c4",
  "free-living"              = "#EAEAEA"
)
