#----------------------10. Environmental Data-----------------------------------
#Contents-----
#1. Load libraries and themes
#2. Load data
#3. Mandatory residue monitoring
#4. Environmental Fungicides
#5.

#1. Load libraries and themes----------
library(tidyverse)
library(patchwork)
library(scales)
library(extrafont)
library(extrafontdb)


BELMAP_colourscheme <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c",
"#fdbf6f","#ff7f00","#cab2d6")

moiras_graph_theme<- function(..., base_size = 12){theme(
  panel.background = element_rect(fill="#f7f7f7", colour = "#001f3f"), #transparent panel bg
  plot.background = element_rect(fill='#f7f7f7', color=NA), #transparent plot bg
  panel.grid.major.y =  element_line(colour = "#ccc6c4"), #grey y major gridlines
  panel.grid.major.x =  element_blank(), #remove major gridlines
  panel.grid.minor = element_blank(), #remove minor gridlines
  # legend.background = element_rect(fill='#001f3f'), #transparent legend bg
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 12),
  #  legend.box.background = element_rect(fill='#001f3f'),
  text = element_text( family = "Gill Sans MT"),
  axis.text = element_text(size = 14),
  #plot.background=element_blank(),#, , size = 10, family = "calibri"
  axis.title = element_text(size = 16),
  axis.text.x = element_text(angle = 90,vjust = 1, hjust = -1)
)}

#2. Load data--------------
FAO_fungicide<- read.csv("\\\\sciensano.be/FS/1220_BactDis_Scientist/Scientists/Moira/BELMAP 2023/Data/3.Environment/FAOSTAT_data_en_7-6-2023.csv",
                         sep = ",", header = TRUE) %>%
  dplyr::select(Item, Year, Value)

#3. Mandatory residue monitoring--------------

#4. Environmental Fungicides--------------


(FAO_fungicide_barchart <- FAO_fungicide %>%
  filter(!grepl("ungicide",Item)) %>%
  ggplot(aes(x=Year,y=Value))+
  geom_bar(aes(fill=Item), stat="identity", position = "stack")+
  scale_fill_manual(values = BELMAP_colourscheme)+
  scale_x_continuous(breaks = seq(2011,2020,by=1))+
   labs(y="Agricultural use (tonnes)")
)




#5. --------------