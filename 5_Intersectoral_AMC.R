#-------------3. Intersector AMC comparison ------------------------------------
#Contents-----
#1. Load libraries
#2. Calculate Human biomass and make human mg/kg dataframe
#3. Load data
#4. Comparative AMC - analysis
#5. Comparative AMC - graphs

#1. Load libraries ------------------------------------------------------------
library(tidyverse)
library(patchwork)
library(jpeg)
library(scales)
library(extrafont)
library(remotes)
library(extrafontdb)
library(magrittr)
library(AER)
library(AICcmodavg)
library(pscl)
library(MASS)
library(COMPoissonReg)
library(performance)
library(grid)
#library(emojis)
library(ggtext)
library(ggrepel)
library(showtext)
#load.fontawesome()
# install.packages("showtext")

BELMAP_colourscheme <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c",
                         "#fdbf6f","#ff7f00","#cab2d6")



BELMAP_intersect_AMC_colourscheme<-c("Food producing animal" = "#99d594",
                                     "Veterinary" = "#99d594",
                                     "Human - Hospital" = "#9ecae1",
                                     "Human - Ambulatory" = "#4292c6",
                                     "Human -Ambulatory" = "#4292c6",
                                     "Human" = "#08519c",
                                     "Human - Total" = "#08519c")

BELMAP_location_colourscheme <- c("Belgium"="#fc8d59",
                                  "Neighbours" = "#99d594",
                                  "Hospital Lab"="#fc8d59",
                                  "Private lab" = "#99d594",
                                  "Europe" = "#3288bd")


BELMAP_Human_AntiB_colourscheme <- c("Glycopeptides" = "#33a02c",
                                     "3rd and 4th generation\ncephalosporins" = "#1a1a1a",
                                     "Fluoroquinolones" = "#fdbf6f",
                                     "Polymixins" = "#1a1a1a",
                                     "Piperacillin + INH" = "#1f78b4",
                                     "Linezolid" = "#ff7f00",
                                     "Carbapenems" = "#a6cee3")


# "Glycopeptides" = "#33a02c",
# "3rd and 4th generation\ncephalosporins" = "#1a1a1a",
# "Fluoroquinolones" = "#fdbf6f",
# "Polymixins" = "#1a1a1a",
# "Piperacillin + INH" = "#1f78b4",
# "Linezolid" = "#ff7f00",
# "Carbapenems" = "#a6cee3"
#"#33a02c","#e31a1c","#fdbf6f","#ff7f00","#cab2d6", "#6a3d9a","#ffff99") #



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

# Important step to enable showtext font rendering!
showtext_auto()


# adding icons to figures 
font_add("Gill Sans MT", "GIL_____.TTF")

font_add('fa-solid', '../font_awesome_font_files/fontawesome-free-6.4.0-desktop/otfs/Font Awesome 6 Free-Solid-900.otf')

upward_arrow <- "<span style='font-family:fa-solid'>&#xf062;</span>"
downward_arrow <- "<span style='font-family:fa-solid'>&#xf063;</span>" 
equals <- "<span style='font-family:fa-solid'>&#xf52c;</span>"   #" = " 
oscillate <- "<span style='font-family:fa-solid'>&#xf83e;</span>"   #" ~ "




# #2. Calculate Human biomass and make human mg/kg dataframe-----
# #human population size
# Human_pop_size <- read.csv("Data/1.AMC/Combined_AMC/Eurostat_population_data.csv",
#                               sep = ";", header = TRUE) %>%
#   mutate(Age = case_when(
#     age == "Y_LT1" ~ 0.5,
#     grepl("Y[1-9]",age) ~  as.numeric(str_replace(age,"Y","")),
#     age == "Y_OPEN" ~ 100))
# 
# 
# #human weight
# human_BW_av <- read.csv("Data/1.AMC/Combined_AMC/EFSA_av_human_bodyweight.csv",
#                         sep = ";", header = TRUE)
# 
# 
# # NB. childrens weights taken from JIACRA childrens table
# #as eurostat age is only per year we took the mean for children under 1 i.e. weight
# # (0.-3 months, + weight 3-6 months + 2*weight 6-12 months)/4 
# # = (4.8+6.7+8.8+8.8) / 4 =7.275
# 
# #AMC by mg/kg
# 
# 
# # make human AMC datasheet - mg for antimicrobials
# 
# # load data
# human_AMC_raw <- read.csv("Data/1.AMC/Combined_AMC/ESAC_OUTPUT_VOLUME.csv") %>%
#   dplyr::select(-X)
# 
# #total volume
# 
# names(human_AMC_raw)
# 
# unique(human_AMC_raw$TOTAL_VOLUME_UNIT)
# # all mg so can just sum volumes
# 
# human_AMC_total <- human_AMC_raw %>%
#   filter(grepl("^J01", ATCCode)) %>%
#   mutate(Year = DateUsedForStatisticsYear) %>%
#   group_by(Year) %>%
#   mutate(National = sum(TOTAL_VOLUME)) %>%
#   ungroup() %>%
#   group_by(Year, Sector, National) %>%
#   summarise(Volume = sum(TOTAL_VOLUME)) %>%
#   ungroup() %>%
#   mutate(Measurement = "Total") %>%
#   pivot_wider(names_from = Sector, values_from = Volume) %>%
#   pivot_longer(cols = c(HC,AC,National), names_to = "Sector", values_to = "Volume (mg)")
#   
#   
#   
# #quinolones
# 
# human_AMC_quinolones <- human_AMC_raw %>%
#   filter(grepl("^J01M", ATCCode)) %>%
#   mutate(Year = DateUsedForStatisticsYear) %>%
#   group_by(Year) %>%
#   mutate(National = sum(TOTAL_VOLUME)) %>%
#   ungroup() %>%
#   group_by(Year, Sector, National) %>%
#   summarise(Volume = sum(TOTAL_VOLUME)) %>%
#   ungroup() %>%
#   mutate(Measurement = "Quinolones") %>%
#   pivot_wider(names_from = Sector, values_from = Volume) %>%
#   pivot_longer(cols = c(HC,AC,National), names_to = "Sector", values_to = "Volume (mg)")
# 
# 
# #polymyxins
# 
# human_AMC_polymixins <- human_AMC_raw %>%
#   filter(grepl("^J01DH", ATCCode)) %>%
#   mutate(Year = DateUsedForStatisticsYear) %>%
#   group_by(Year) %>%
#   mutate(National = sum(TOTAL_VOLUME)) %>%
#   ungroup() %>%
#   group_by(Year, Sector, National) %>%
#   summarise(Volume = sum(TOTAL_VOLUME)) %>%
#   ungroup() %>%
#   mutate(Measurement = "Polymyxins") %>%
#   pivot_wider(names_from = Sector, values_from = Volume) %>%
#   pivot_longer(cols = c(HC,AC,National), names_to = "Sector", values_to = "Volume (mg)")
# 
# #penicillins
# 
# human_AMC_penicillins <- human_AMC_raw %>%
#   filter(grepl("^J01C", ATCCode)) %>%
#   mutate(Year = DateUsedForStatisticsYear) %>%
#   group_by(Year) %>%
#   mutate(National = sum(TOTAL_VOLUME)) %>%
#   ungroup() %>%
#   group_by(Year, Sector, National) %>%
#   summarise(Volume = sum(TOTAL_VOLUME)) %>%
#   ungroup() %>%
#   mutate(Measurement = "Penicillins") %>%
#   pivot_wider(names_from = Sector, values_from = Volume) %>%
#   pivot_longer(cols = c(HC,AC,National), names_to = "Sector", values_to = "Volume (mg)")
# 
# #macrolides
# human_AMC_macrolides <- human_AMC_raw %>%
#   filter(grepl("^J01FA", ATCCode)) %>%
#   mutate(Year = DateUsedForStatisticsYear) %>%
#   group_by(Year) %>%
#   mutate(National = sum(TOTAL_VOLUME)) %>%
#   ungroup() %>%
#   group_by(Year, Sector, National) %>%
#   summarise(Volume = sum(TOTAL_VOLUME)) %>%
#   ungroup() %>%
#   mutate(Measurement = "Macrolides") %>%
#   pivot_wider(names_from = Sector, values_from = Volume) %>%
#   pivot_longer(cols = c(HC,AC,National), names_to = "Sector", values_to = "Volume (mg)")
# 
# #tetracyclines
# 
# human_AMC_tetracyclines <- human_AMC_raw %>%
#   filter(grepl("^J01A", ATCCode)) %>%
#   mutate(Year = DateUsedForStatisticsYear) %>%
#   group_by(Year) %>%
#   mutate(National = sum(TOTAL_VOLUME)) %>%
#   ungroup() %>%
#   group_by(Year, Sector, National) %>%
#   summarise(Volume = sum(TOTAL_VOLUME)) %>%
#   ungroup() %>%
#   mutate(Measurement = "Tetracyclines") %>%
#   pivot_wider(names_from = Sector, values_from = Volume) %>%
#   pivot_longer(cols = c(HC,AC,National), names_to = "Sector", values_to = "Volume (mg)")
# 
# #3rd and 4th generation cephalosporins
# 
# human_AMC_3GC <- human_AMC_raw %>%
#   filter(grepl("^J01DD|^J01DE", ATCCode)) %>%
#   mutate(Year = DateUsedForStatisticsYear) %>%
#   group_by(Year) %>%
#   mutate(National = sum(TOTAL_VOLUME)) %>%
#   ungroup() %>%
#   group_by(Year, Sector, National) %>%
#   summarise(Volume = sum(TOTAL_VOLUME)) %>%
#   ungroup() %>%
#   mutate(Measurement = "3rd and 4th gen cephalosporins") %>%
#   pivot_wider(names_from = Sector, values_from = Volume) %>%
#   pivot_longer(cols = c(HC,AC,National), names_to = "Sector", values_to = "Volume (mg)")
# 
# 
# #carbapenems
# 
# human_AMC_carbapenems <- human_AMC_raw %>%
#   filter(grepl("^J01DH", ATCCode)) %>%
#   mutate(Year = DateUsedForStatisticsYear) %>%
#   group_by(Year) %>%
#   mutate(National = sum(TOTAL_VOLUME)) %>%
#   ungroup() %>%
#   group_by(Year, Sector, National) %>%
#   summarise(Volume = sum(TOTAL_VOLUME)) %>%
#   ungroup() %>%
#   mutate(Measurement = "Carbapenems") %>%
#   pivot_wider(names_from = Sector, values_from = Volume) %>%
#   pivot_longer(cols = c(HC,AC,National), names_to = "Sector", values_to = "Volume (mg)")
# 
# 
# Human_AMC_volumes <- rbind(human_AMC_total,human_AMC_3GC,human_AMC_carbapenems,human_AMC_macrolides,
#                            human_AMC_penicillins,human_AMC_polymixins,human_AMC_quinolones,human_AMC_tetracyclines)
# 
# 



# #1st and 2nd generation cephalosporins
# 
# human_AMC_1GC <- human_AMC_raw %>%
#   filter(grepl("^J01DB|^J01DC", ATCCode)) %>%
#   mutate(Year = DateUsedForStatisticsYear) %>%
#   group_by(Year) %>%
#   mutate(National = sum(TOTAL_VOLUME)) %>%
#   ungroup() %>%
#   group_by(Year, Sector, National) %>%
#   summarise(Volume = sum(TOTAL_VOLUME)) %>%
#   ungroup() %>%
#   mutate(Measurement = "1st and 2nd gen cephalosporins") %>%
#   pivot_wider(names_from = Sector, values_from = Volume) %>%
#   pivot_longer(cols = c(HC,AC,National), names_to = "Sector", values_to = "Volume (mg)")
# 
# human_AMC_trimeth <- human_AMC_raw %>%
#   filter(grepl("^J01E", ATCCode)) %>%
#   mutate(Year = DateUsedForStatisticsYear) %>%
#   group_by(Year) %>%
#   mutate(National = sum(TOTAL_VOLUME)) %>%
#   ungroup() %>%
#   group_by(Year, Sector, National) %>%
#   summarise(Volume = sum(TOTAL_VOLUME)) %>%
#   ungroup() %>%
#   mutate(Measurement = "Trimethoprim & Sulphonam") %>%
#   pivot_wider(names_from = Sector, values_from = Volume) %>%
#   pivot_longer(cols = c(HC,AC,National), names_to = "Sector", values_to = "Volume (mg)")
# 
# 
# human_AMC_amino <- human_AMC_raw %>%
#   filter(grepl("^J01G", ATCCode)) %>%
#   mutate(Year = DateUsedForStatisticsYear) %>%
#   group_by(Year) %>%
#   mutate(National = sum(TOTAL_VOLUME)) %>%
#   ungroup() %>%
#   group_by(Year, Sector, National) %>%
#   summarise(Volume = sum(TOTAL_VOLUME)) %>%
#   ungroup() %>%
#   mutate(Measurement = "Aminoglycosides") %>%
#   pivot_wider(names_from = Sector, values_from = Volume) %>%
#   pivot_longer(cols = c(HC,AC,National), names_to = "Sector", values_to = "Volume (mg)")
# 
# 
# 
# Human_AMC_volumes2 <- rbind(human_AMC_1GC,human_AMC_trimeth,human_AMC_amino)
# 




# #2b. Calculate Human biomass------
# 
# human_biomass_annual <-  left_join(Human_pop_size, human_BW_av, by = c("Age","sex")) %>%
#   mutate(age_biomass = OBS_VALUE*weight.kg) %>%
#   mutate(Year = TIME_PERIOD) %>%
#   group_by(Year) %>%
#   summarise(annual_bodymass = sum(age_biomass))%>%
#   ungroup()
# 
# #2b. Make human database------
# 
# human_mg_per_kg2 <- left_join(Human_AMC_volumes2,human_biomass_annual, by = c("Year")) %>%
#   mutate(`Volume (mg/kg)` = `Volume (mg)` / annual_bodymass,
#           Unit = "mg/kg biomass",
#          `Biomass (tonnes)` = annual_bodymass/1000,
#          `Volume (kg)` = `Volume (mg)`/1000000)%>%
#   arrange(Sector,Measurement,Year) %>%
#   dplyr::select(Year,Measurement,`Volume (mg/kg)`,Unit,Sector, `Volume (kg)`, `Biomass (tonnes)`)



#write.csv2(human_mg_per_kg,"Data/1.AMC/Combined_AMC/human_mg_per_kg.csv")


#write.csv2(human_mg_per_kg2,"Data/1.AMC/Combined_AMC/human_mg_per_kg2.csv")


# 3. Load data ---------------------

Human_vet_AMC <- read.csv2("Data/1.AMC/Combined_AMC/AMC_vet_human_mg_kg.csv") %>%
  mutate(Category = paste(Sector,Measurement,sep = "_"))


# #4. Comparative AMC - for loop for totals------
# names(Human_vet_AMC)
# 
# 
# start_data_frame <- tibble(
#   "Year" = "",
#   "Measurement" = "",
#   "Value"= "",
#   "Unit" = "",
#   "Sector" = "",
#   "kg.active.product" = "",
#   "tonnes.of.biomass" = "",
#   "Category" = "",
#   "Year_simple"= "",
#   "icon" = "",
#   "signif" = ""
# )
# 
# write_csv(start_data_frame,file = "5_Intersectoral_AMC_analysis_outcomes.csv")
# 
# # list AMC indicators
# 
# AMC_categories<- unique(Human_vet_AMC[Human_vet_AMC$Measurement == "Total",]$Category)
# 
# #mkae output test results file
# capture.output(print("Analysis of Comparative AMC data"),file = "5_Intersectoral_AMC_correlations/AMC_data_analysis.txt")
# 
# 
# #make list for overdispersion plots
# qq_plot_list<- list()
# #open for loop AMC ---------------------------
# 
# for(i in AMC_categories){
# 
#   # 1. check assumptions for Pearsons
#   dataset_AMC<- Human_vet_AMC %>%
#     filter(grepl(i,Category,  fixed = TRUE))
#   min_year = min(dataset_AMC$Year)
#   max_year = max(dataset_AMC$Year)
# 
#   dataset_AMC <- dataset_AMC  %>%
#   mutate(Year_simple = Year - min_year) %>%
#     distinct()
# 
# 
# 
#   capture.output(print(paste("Analysis for",i,sep = " : ")),file = "5_Intersectoral_AMC_correlations/AMC_data_analysis.txt", append = TRUE)
# 
#   #save qq plot
# 
#   #Normality check
#   # qqnorm(dataset_AMC$Year_simple, pch =1, frame = FALSE)
#   # qqline(dataset_AMC$Year_simple, col = "steelblue", lwd = 2)
#   #
#   # qqnorm(dataset_AMC$Value, pch =1, frame = FALSE)
#   # qqline(dataset_AMC$Value, col = "steelblue", lwd = 2)
# 
# 
#   # qqPlot(dataset_AMC$Value)
#   # qqPlot_value<- recordPlot()
#   #
#   # qqPlot(dataset_AMC$Year)
#   # qqPlot_year<- recordPlot()
#   #
#   # #assign to local environment
#   # assign(paste(i,"qqplot_value"),print(qqPlot_value))
#   # assign(paste(i,"qqplot_year"),print(qqPlot_year))
#   #
#   # #Shapiro test: if the  p-values are greater than the significance level 0.05 => the distribution of the data are not significantly different from normal distribution
#   # # Shapiro-Wilk normality test for Year_simple
#   #
#   Shapiro_year <- shapiro.test(dataset_AMC$Year_simple) #p-value = 0.8698
#   capture.output(print("Testing Normality of year"),file = "5_Intersectoral_AMC_correlations/AMC_data_analysis.txt", append = TRUE)
#   capture.output(print(Shapiro_year),file = "5_Intersectoral_AMC_correlations/AMC_data_analysis.txt", append = TRUE)
# 
# 
#   # Shapiro-Wilk normality test for Value
#   Shapiro_value <- shapiro.test(dataset_AMC$Value) #p-value = 0.4784
# 
#   capture.output(print("Testing Normality of Values"),file = "5_Intersectoral_AMC_correlations/AMC_data_analysis.txt", append = TRUE)
#   capture.output(print(Shapiro_value),file = "5_Intersectoral_AMC_correlations/AMC_data_analysis.txt", append = TRUE)
# 
# 
#   # 2. calculate correlation tests
# 
#   #Pearson (parametric test, assumes linearity and normality) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <cor> and a p-value of <p-value>.
#   Pcor_Value <- cor(dataset_AMC$Year_simple, dataset_AMC$Value, method = "pearson")
#   Pcortest_Value <- cor.test(dataset_AMC$Year_simple, dataset_AMC$Value, method = "pearson")
#   Pcor_Value
#   #        cor -0.9856587; p-value = 2.905e-08
#   Pcortest_Value
#   #        cor -0.9856587; p-value = 2.905e-08
# 
#   capture.output(print("#Pearson (parametric test, assumes linearity and normality) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <cor> and a p-value of <p-value>.
# "),file = "5_Intersectoral_AMC_correlations/AMC_data_analysis.txt", append = TRUE)
#   capture.output(print(Pcortest_Value),file = "5_Intersectoral_AMC_correlations/AMC_data_analysis.txt", append = TRUE)
# 
# 
#   #Kendall (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <tau> and a p-value of <p-value>.
#   Kcor_Value <- cor(dataset_AMC$Year_simple, dataset_AMC$Value, method = "kendall")
#   Kcortest_Value <- cor.test(dataset_AMC$Year_simple, dataset_AMC$Value, method = "kendall")
#   Kcor_Value
#   # [1] -0.9272727
# 
#   Kcortest_Value
#   #  tau -0.9272727;  p-value = 3.257e-06
# 
#   capture.output(print("#Kendall (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <tau> and a p-value of <p-value>.
# "),file = "5_Intersectoral_AMC_correlations/AMC_data_analysis.txt",  append = TRUE)
#   capture.output(print(Kcortest_Value),file = "5_Intersectoral_AMC_correlations/AMC_data_analysis.txt", append = TRUE)
# 
# 
# 
#   #Spearman (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <rho> and a p-value of <p-value>.
#   Scor_Value <- cor(dataset_AMC$Year_simple, dataset_AMC$Value, method = "spearman")
#   Scortest_Value <- cor.test(dataset_AMC$Year_simple, dataset_AMC$Value, method = "spearman")
#   Scor_Value
#   #[1] -0.9818182
# 
#   Scortest_Value
# 
#   capture.output(print("#Spearman (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <rho> and a p-value of <p-value>.
# "),file = "5_Intersectoral_AMC_correlations/AMC_data_analysis.txt",  append = TRUE)
#   capture.output(print(Scortest_Value),file = "5_Intersectoral_AMC_correlations/AMC_data_analysis.txt", append = TRUE)
# 
# 
#   # manually check that results of correlation tests are aligned
# 
# 
#   # 3. select correlation and save icon-------------
#   #select model
#   if(Shapiro_value$p.value >  0.05 ){   # if Shapiro > 0.05 - can be treated as normally distributed
#     select_correlation_value <-  Pcortest_Value$estimate
#     select_correlation_pvalue <-  Pcortest_Value$p.value
#   }else if(Shapiro_value$p.value <  0.05){ # if Shapiro < 0.05 - not normally distributed
#     select_correlation_value <-  Scortest_Value$estimate
#     select_correlation_pvalue <-  Scortest_Value$p.value
#   }else{
#     print(paste("ERROR no model selected for",i))
#   }
# 
#   #if signif - extract direction of model
#   if(select_correlation_pvalue<0.05){
#     if(select_correlation_value > 0){
#       model_icon <- "upward_arrow"
#       signif = case_when(
#         select_correlation_pvalue > 0.01 && select_correlation_pvalue < 0.05 ~ "*",
#         select_correlation_pvalue > 0.001 && select_correlation_pvalue < 0.01 ~ "**",
#         select_correlation_pvalue < 0.001 ~ "***"
#       )
#     }else if(select_correlation_value < 0){
#       model_icon = "downward_arrow"
#       signif = case_when(
#         select_correlation_pvalue > 0.01 && select_correlation_pvalue < 0.05 ~ "*",
#         select_correlation_pvalue > 0.001 && select_correlation_pvalue < 0.01 ~ "**",
#         select_correlation_pvalue < 0.001 ~ "***"
#       )
#     }
#   } else if(select_correlation_pvalue>0.05){
#     if((max(dataset_AMC$Value)-min(dataset_AMC$Value))/mean(dataset_AMC$Value) > 0.25){
#       model_icon <- "oscilate"
#       signif = ""
#     }else if((max(dataset_AMC$Value)-min(dataset_AMC$Value))/mean(dataset_AMC$Value) < 0.25){
#       model_icon <- "equals"
#       signif <- ""
#     }
#   }
# 
# 
#   # Add labels-------------------------
#   dataset_AMC_analysis <- dataset_AMC %>%
#     mutate(icon = model_icon) %>%
#     mutate(signif = signif)
# 
#   # Save output to dataframe for graphs ---------------------
# 
#   write_csv(dataset_AMC_analysis,file = "5_Intersectoral_AMC_analysis_outcomes.csv", append = TRUE)
# 
# 
#   #close for loop  ---------------------------
# 
# }




#5. Comparative AMC - graphs-----


comparative_graph_data <- read_csv("5_Intersectoral_AMC_analysis_outcomes.csv") %>%
  mutate(label_icon = case_when(
    icon == "upward_arrow" ~ upward_arrow,
    icon == "downward_arrow" ~ downward_arrow,
    icon == "equals" ~ equals,
    icon == "oscilate" ~ oscillate,
    .default = ""
  )) %>%
  mutate(label = if_else((is.na(signif) | signif == ""), label_icon, paste(label_icon,signif,sep=" "))) %>%
  mutate(label = if_else(Year == 2021, label,NA)) %>%
  filter(Year > 2011) %>%
  filter( Year < 2022) %>%
  mutate(y_coord = case_when(
    Sector == "Veterinary" ~ 60,
    Sector == "Human -Ambulatory" ~ 90,
    Sector == "Human - Hospital" ~ 30,
    Sector == "Human" ~ 120)) %>%
  mutate(Species = if_else(grepl("Veterinary", Sector), "Veterinary", "Human")) %>%
  mutate(Sector = if_else(Sector == "Human", "Human - Total", Sector)) %>%
  mutate(Sector = factor(Sector, 
                         levels = c("Human - Total", "Human -Ambulatory","Human - Hospital", "Veterinary")))

#make line for 2021 to add icon

add_line <-  tibble(
  "Year" = c(2021,2021),
  "Measurement" = c("none","none"),
  "Value" = c(0,0),
  "Unit" = c("mg","mg"),
  "Sector"  = c("Zoo","Zoo"),
  "kg.active.product" = c(NA,NA),
  "tonnes.of.biomass" = c(NA,NA),
  "Category"= c(NA,NA),
  "Year_simple"= c(NA,NA),
  "icon"= c(NA,NA),
  "signif"= c(NA,NA),
  "label_icon"= c(NA,NA),
  "label" = c(NA,NA),
  "y_coord"= c(NA,NA),
  "Species"= c("You","Zoo"))


comparative_graph_data1<-rbind(comparative_graph_data,add_line)

(Comparative_graph_AMC <- comparative_graph_data1 %>%
    ggplot()+
    geom_bar(data = comparative_graph_data1[!(comparative_graph_data1$Sector == "Human"),], aes(x=Species, y = as.numeric(Value), fill = Sector, colour = Species), stat="identity", position = "stack")+
    scale_fill_manual(values = BELMAP_intersect_AMC_colourscheme )+
    scale_colour_manual(values = BELMAP_intersect_AMC_colourscheme,breaks = NULL)+
    facet_wrap(~Year, nrow=1, switch = "x", scales = "free_x")+
    #  scale_alpha_manual(values = c("Primary" = 1, "Secondary" = 0.5))+
    # geom_richtext( size = 16, hjust = 0, label.colour = NA) +
    #  scale_x_continuous(limits = c(0,4), breaks ="none")+
    scale_y_continuous(limits = c(-1,175), breaks = seq(0,150,50))+
    geom_richtext(aes(x = 3.5, y = y_coord,
                      label = label,
                      fill = Sector), stat = "unique",
                  colour= "white", show.legend = FALSE, size = 4)+
    # fill = NA, label.color = NA, # remove background and outline
    # label.padding = grid::unit(rep(0, 4), "pt")) +
    #  geom_text(aes(x = Year, y = 2, label = Sample_size.x), na.rm = TRUE, size = 2.5)+  #make this additional for interactive report
    labs(y= "Antimicrobial Consumption (mg/kg biomass)", x = "Year")&
    moiras_graph_theme()+
    theme(axis.text.x = element_blank(),
          strip.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text.x.bottom = element_text(size = 14)
    ))





# by category 
#-bar chart for 2021 and box plots of variation--------

category_AMC_2021 <- Human_vet_AMC %>%
  filter(!grepl("Total",Measurement)) %>%
  filter(Year == 2021) %>%
  filter(Sector != "Human") %>% 
  mutate(Sector = if_else(Sector == "Human", "Human - Total", Sector)) %>%
  mutate(Species = if_else(grepl("Veterinary", Sector), "Veterinary", "Human"))


(category_graph1<-    category_AMC_2021 %>%
    mutate(volume = if_else(Species == "Veterinary", -Value, Value)) %>%
    ggplot() +
    geom_bar(aes(x = Measurement, y = volume, fill = Sector, colour = Species), stat = "identity",
             position = "stack")+
    geom_vline(xintercept = 0, colour = "black")+
    scale_fill_manual(values = BELMAP_intersect_AMC_colourscheme)+
    scale_colour_manual(values = BELMAP_intersect_AMC_colourscheme, guide = 'none')+
    coord_flip()+
    scale_y_continuous(limits = c(-85,85),breaks = c(-80,40,0,40,80),labels = rep("",5))+
    geom_text(aes(y = -80, x = -0.1, label = "80"), size = 3.5,  angle = 90, colour = "black", hjust = 1.5)+
    geom_text(aes(y = -40, x = -0.1, label = "40"), size = 3.5,  angle = 90, colour = "black", hjust = 1.5)+
    geom_text(aes(y = 0, x = -0.1, label = "0"), size = 3.5,  angle = 90, colour = "black", hjust = 1.5)+
    geom_text(aes(y = 40, x = -0.1, label = "40"), size = 3.5,  angle = 90, colour = "black", hjust = 1.5)+
    geom_text(aes(y = 80, x = -0.1, label = "80"), size = 3.5,  angle = 90, colour = "black", hjust = 1.5)+
    labs(x = "Antimicrobial group", y = "Consumption (mg/kg)")+
    moiras_graph_theme())



