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
library(lubridate)
library(patchwork)
library(viridis)


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


# function for graphs - legend into space

#shift legend to fill blanks of facet_grid
shift_legend <- function(p){
  
  # check if p is a valid object
  if(!"gtable" %in% class(p)){
    if("ggplot" %in% class(p)){
      gp <- ggplotGrob(p) # convert to grob
    } else {
      message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
      return(p)
    }
  } else {
    gp <- p
  }
  
  # check for unfilled facet panels
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  if(length(empty.facet.panels) == 0){
    message("There are no unfilled facet panels to shift legend into. Returning original plot.")
    return(p)
  }
  
  # establish extent of unfilled facet panels (including any axis cells in between)
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  empty.facet.panels <- list(min(empty.facet.panels[["t"]]), min(empty.facet.panels[["l"]]),
                             max(empty.facet.panels[["b"]]), max(empty.facet.panels[["r"]]))
  names(empty.facet.panels) <- c("t", "l", "b", "r")
  
  # extract legend & copy over to location of unfilled facet panels
  guide.grob <- which(gp[["layout"]][["name"]] == "guide-box")
  if(length(guide.grob) == 0){
    message("There is no legend present. Returning original plot.")
    return(p)
  }
  gp <- gtable_add_grob(x = gp,
                        grobs = gp[["grobs"]][[guide.grob]],
                        t = empty.facet.panels[["t"]],
                        l = empty.facet.panels[["l"]],
                        b = empty.facet.panels[["b"]],
                        r = empty.facet.panels[["r"]],
                        name = "new-guide-box")
  
  # squash the original guide box's row / column (whichever applicable)
  # & empty its cell
  guide.grob <- gp[["layout"]][guide.grob, ]
  if(guide.grob[["l"]] == guide.grob[["r"]]){
    gp <- gtable_squash_cols(gp, cols = guide.grob[["l"]])
  }
  if(guide.grob[["t"]] == guide.grob[["b"]]){
    gp <- gtable_squash_rows(gp, rows = guide.grob[["t"]])
  }
  gp <- gtable_remove_grobs(gp, "guide-box")
  
  return(gp)
}





#2. Load data--------------
FAO_fungicide<- read.csv("\\\\sciensano.be/FS/1220_BactDis_Scientist/Scientists/Moira/BELMAP 2023/Data/3.Environment/FAOSTAT_data_en_7-6-2023.csv",
                         sep = ",", header = TRUE) %>%
  dplyr::select(Item, Year, Value)

#load data

water_raw_data<- read_csv2("\\\\sciensano.be/FS/1220_BactDis_Scientist/Scientists/Moira/BELMAP 2023/Data/3.Environment/water_antimicrobial_residues.csv")


flanders_raw_data1<- read_csv2("\\\\sciensano.be/FS/1220_BactDis_Scientist/Scientists/Moira/BELMAP 2023/Data/3.Environment/VMM_water_data_raw.csv")


flanders_measurement_points<- read_csv2("\\\\sciensano.be/FS/1220_BactDis_Scientist/Scientists/Moira/BELMAP 2023/Data/3.Environment/VMM_meetpunt_locations.csv")%>%
  filter(regio == "Vlaanderen") %>%
  select(meetplaats,provincie) %>%
  mutate(meetplaats = as.character(meetplaats))


flanders_raw_data<- left_join(flanders_raw_data1,flanders_measurement_points, by = "meetplaats")

#3. Mandatory residue monitoring--------------

#tidy data

# names(flanders_raw_data)
# 
# names(flanders_measurement_points)
# 
# unique(flanders_measurement_points$regio)


water_data<-  water_raw_data %>%
  pivot_longer(cols = -c(Date, Antibiotic), names_to = "Location", values_to = "concentration") %>%
  mutate(concentration = str_replace(concentration,",","."))%>%
    mutate(exact_value = if_else(grepl("<", concentration), "< Value Recorded", "Value Recorded"),
           Date = dmy(Date),
           Organisation = if_else(grepl("_",Location), "Brussels Environment", "SPW")) %>%
  mutate(`Concentration (ug)` = as.numeric(str_replace(concentration,"<","")),
         Antibiotic = case_when(
           grepl("Amoxicilline", Antibiotic) ~ "Amoxicillin",
           grepl("Azitromycine ", Antibiotic) ~  "Azithromycin",
           grepl("Ciprofloxacine ", Antibiotic) ~  "Ciprofloxacin",
           grepl("Claritromycine", Antibiotic) ~  "Clarithromycin",
           grepl("Erytromycine ", Antibiotic) ~  "Erythromycin",
           .default = Antibiotic))  %>%
  dplyr::select(-concentration)


flanders_water_data <- flanders_raw_data %>%
  filter(!is.na(waterlichaam))%>%
  unite("Amoxicillin", Amoxicillin.ng.L.sign:Amoxicillin.ng.L., remove = TRUE,sep = " ") %>%
  unite("Azithromycin", Azithromyc.ng.L.sign:Azithromyc.ng.L., remove = TRUE,sep = " ") %>%
  unite("Ciprofloxacin", Ciprofloxacin.ng.L.sign:Ciprofloxacin.ng.L., remove = TRUE,sep = " ") %>%
  unite("Clarithromycin", Clarithromyc.ng.L.sign:Clarithromyc.ng.L., remove = TRUE,sep = " ") %>%
  unite("Erythromycin", Erythromyc.ng.L.sign:Erythromyc.ng.L, remove = TRUE,sep = " ") %>%
  pivot_longer(cols = -c(meetplaats , waterlichaam,X,Y,Deelmonster,Datum,provincie), names_to = "Antibiotic", values_to = "concentration") %>%
  filter(!is.na(provincie)) %>%
  mutate(exact_value = if_else(grepl("< ", concentration), "< Value Recorded", "Value Recorded"),
         Date = dmy(Datum),
         Location =  case_when(
           provincie == "Antwerpen" ~ "Antwerp",
           provincie == "Limburg" ~ "Limburg",
           provincie == "Oost-Vlaanderen" ~ "East Flanders",
           provincie == "West-Vlaanderen" ~ "West Flanders",
           provincie == "Vlaams-Brabant" ~ "Flemish Brabant"),
         Organisation = "VMM",
        `Concentration (ug)` = as.numeric(str_replace(concentration,"< ",""))/1000) %>%
  dplyr::select(-c(meetplaats, waterlichaam, X, Y, Deelmonster,Datum,concentration, provincie))
  
long_env_data <- rbind(water_data,flanders_water_data) %>%
  filter(!is.na(Date))



#add PNEC

PNECS <- tibble(
  Antibiotic = c("Clarithromycin", "Azithromycin", "Amoxicillin", "Ciprofloxacin", "Erythromycin"),
  PNEC = c(0.12, 0.019, 0.078, 0.089, 0.2)
)


#add in EU averages
Europe_data<- read.csv("\\\\sciensano.be/FS/1220_BactDis_Scientist/Scientists/Moira/BELMAP 2023/Data/3.Environment/Europe_data.csv", sep = ";") %>%
  select(Antibiotic, MEC)


#make range and mean:

Europe_data_clean <- Europe_data %>%
  group_by(Antibiotic) %>%
  summarise(mean_MEC = mean(MEC), min_MEC = min(MEC), max_MEC = max(MEC)) %>%
  ungroup() %>%
  mutate(Date = "01-01-2015")%>%
  mutate(Date_read = dmy(Date)) %>%
  filter(!is.na(mean_MEC))


plot_list= list()

plot_colours <-  c("#a6cee3","#1f78b4","#fb9a99","#b2df8a","#33a02c")

for (i in unique(long_env_data$Organisation)) {
  data_plot<- long_env_data %>%
    filter(Organisation == i)

#   if(grepl("VMM", unique(data_plot$Organisation))){
#     plot_colours <- viridis(43)
#   }else{
#     plot_colours <-  c("#a6cee3","#1f78b4","#fb9a99","#b2df8a","#33a02c")    
#   }
 
  
  
  first_plot_trial <- data_plot %>%
    ggplot()+
    geom_jitter(aes(x=Date, y=`Concentration (ug)`, colour =  Location, fill =Location, shape = exact_value),size = 2, width = 20, height = 0.002)+
    geom_errorbar(data = Europe_data_clean, aes(x = Date_read, ymin=min_MEC, ymax=max_MEC), size =1, width=0, colour = "blue")+
    geom_point(data = Europe_data_clean, aes(x = Date_read, y=mean_MEC), size =2, colour = "red")+
    scale_colour_manual(values = plot_colours)+
    scale_fill_manual(values = plot_colours)+
    scale_shape_manual(values=c(1, 19))+
    scale_y_continuous( trans='log10', limits = c(0.0001,100))+
    facet_wrap(~Antibiotic)+
    labs(x= "Time", y="Concentration (ug/l)", shape = "",title = i)+
    geom_hline(data = PNECS, aes(yintercept=PNEC),linetype='dotted')+
    theme_classic()+
    theme(strip.text = element_text(size = 12)
          #strip.background = element_blank()
    )
  plot_list[[i]] <- first_plot_trial
  
}


plot_list[[1]] / plot_list[[2]] / plot_list[[3]] +
  plot_layout(guides = "collect")


##Print plots to tiff for report ---------

tiff("Brussel_env1.tiff", width = 12, height = 7, units = "in", res = 800)
# Call plot
grid.draw(shift_legend(plot_list[[1]]))
# Closing the graphical device
dev.off()

tiff("SPW_1.tiff", width = 12, height = 7, units = "in", res = 800)
# Call plot
grid.draw(shift_legend(plot_list[[2]]))
# Closing the graphical device
dev.off()

tiff("VMM_1.tiff", width = 12, height = 7, units = "in", res = 800)
# Call plot
grid.draw(shift_legend(plot_list[[3]]))
# Closing the graphical device
dev.off()




#4. Environmental Fungicides--------------


(FAO_fungicide_barchart <- FAO_fungicide %>%
  filter(!grepl("ungicide",Item)) %>%
  ggplot(aes(x=Year,y=Value))+
  geom_bar(aes(fill=Item), stat="identity", position = "stack")+
  scale_fill_manual(values = BELMAP_colourscheme)+
  scale_x_continuous(breaks = seq(2011,2020,by=1))+
   labs(y="Agricultural use (tonnes)")
)


# human Fungicide for comparison

human_AMC_raw <- read.csv("\\\\sciensano.be/FS/1220_BactDis_Scientist/Scientists/Moira/BELMAP 2023/Data/1.AMC/Combined_AMC/ESAC_OUTPUT_VOLUME.csv") %>%
  dplyr::select(-X)


#azole atc codes:
Azole_codes <- c("D01AC", "J02AB", "J02AC")

human_AMC_azole_raw <-human_AMC_raw %>%
  filter(grepl(paste(Azole_codes,collapse = "|"),ATCCode)) %>%
  mutate(Year = DateUsedForStatisticsYear ) %>%
  group_by(Year) %>%
  summarise(Volume.mg = sum(TOTAL_VOLUME)) %>%
  ungroup %>%
  mutate(Volume.tonnes = round(Volume.mg/1000000000)) %>%
  dplyr::select(Year, Volume.tonnes) %>%
  mutate(Sector = "Human medicine")

# Azole compared to human :


FAO_fungicide_azole <- FAO_fungicide %>%
  filter(Item == "Fung & Bact – Triazoles, diazoles") %>%
  mutate(Volume.tonnes = round(Value)) %>%
  dplyr::select(Year, Volume.tonnes) %>%
  mutate(Sector = "Agricultural")

Azoles_uncount <- rbind(human_AMC_azole_raw,FAO_fungicide_azole) %>%
  uncount(data = ., weights = Volume.tonnes) %>%
  filter(as.numeric(Year) > 2011) %>%
  filter(as.numeric(Year) < 2021)


Azoles <- rbind(human_AMC_azole_raw,FAO_fungicide_azole) %>%
  filter(as.numeric(Year) > 2011) %>%
  filter(as.numeric(Year) < 2021)
# 
# (Azole_figure1 <- Azoles_uncount  %>%
#   ggplot() +
#     geom_point(data = Azoles_uncount[Azoles_uncount$Sector == "Human medicine",],aes(x = Sector, y = Year, colour = Sector), size = 1)+
#   geom_jitter(data = Azoles_uncount[Azoles_uncount$Sector == "Agricultural",],
#               aes(x = Sector, y = Year, colour = Sector), 
#               size = 1,height = 0.3)+
#     scale_y_continuous(limits = c(2011.5,2020.5), breaks =seq(2012,2020,by = 1))+
#     coord_flip())


(Azole_figure2 <- Azoles  %>%
    ggplot() +
     geom_bar(aes(x = Year, y = Volume.tonnes, fill = Sector), position = "stack", stat = "identity")+
    scale_x_continuous(limits = c(2011.5,2020.5), breaks =seq(2012,2020,by = 1))+
    ylab("Consumption (tonnes)"))

# (Azole_figure3 <- Azoles  %>%
#     ggplot() +
#     geom_point(aes(x = Year, Sector, size = Volume.tonnes, colour = Sector))+
#     scale_x_continuous(limits = c(2011.5,2020.5), breaks =seq(2012,2020,by = 1))+
#     coord_flip())
# 
# 
# Azole_figure1+Azole_figure2





