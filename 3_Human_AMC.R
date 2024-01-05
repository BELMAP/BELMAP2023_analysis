#-------------2. Human AMC ---------------------------------

#This script outlines the analysis for the BELMAP 2023 Chapter 3 - AMC in the Humanerinary Sector

#Contents-----
#1. Load libraries
#2. Load data
#3. Community DID analysis
#4. Community Ratio broad: narrow 
#5. Hospital Value perpatient days
#6. Hospital Ratio broad:narrow
#7. Hospital proportional consumption ATC groups




#1. Load libraries and themes----
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



#2. Load data--------
Human_AMC_Data <- read_csv2("Data/1.AMC/Human_AMC/Human_AMC_2023.csv",
                            na = c("", "NA", "-"),
                            trim_ws = TRUE)%>%
  mutate(Year_simple= Year - 2011)


unique(Human_AMC_Data$Sector)
unique(Human_AMC_Data$Unit)


#3. Community DID analysis-----------------------------------------------------
Community_DID <- Human_AMC_Data %>%
  filter(Sector == "Ambulant",
         Unit == "DDD/1000 inhabitants/day",
         `ATC group` == "J01")


# correlation analysis - Belgium
#make dataset
Community_DID_Belgium <- Community_DID %>%
  filter(Region == "Belgium")

#Normality check
qqnorm(Community_DID_Belgium$Year_simple, pch =1, frame = FALSE)
qqline(Community_DID_Belgium$Year_simple, col = "steelblue", lwd = 2)

qqnorm(Community_DID_Belgium$Value, pch =1, frame = FALSE)
qqline(Community_DID_Belgium$Value, col = "steelblue", lwd = 2)

library(car)
qqPlot(Community_DID_Belgium$Value)   # one outlier

#Shapiro test: if the  p-values are greater than the significance level 0.05 => the distribution of the data are not significantly different from normal distribution
# Shapiro-Wilk normality test for Year_simple
shapiro.test(Community_DID_Belgium$Year_simple) #p = 0,9 --> can treat as normally distributed

# Shapiro-Wilk normality test for Value
shapiro.test(Community_DID_Belgium$Value) # p = 0,08 --> can treat as normally distributed

#Pearson (parametric test, assumes linearity and normality) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <cor> and a p-value of <p-value>.
PcorBelgium <- cor(Community_DID_Belgium$Year_simple, Community_DID_Belgium$Value, method = "pearson")
PcortestBelgium <- cor.test(Community_DID_Belgium$Year_simple, Community_DID_Belgium$Value, method = "pearson")
PcorBelgium
PcortestBelgium

# cor 
# -0.8619936    p-value = 0.0003113

#Kendall (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <tau> and a p-value of <p-value>.
KcorBelgium <- cor(Community_DID_Belgium$Year_simple, Community_DID_Belgium$Value, method = "kendall")
KcortestBelgium <- cor.test(Community_DID_Belgium$Year_simple, Community_DID_Belgium$Value, method = "kendall")
KcorBelgium
KcortestBelgium

# tau 
# -0.7878788    p-value = 0.0001074

#Spearman (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <rho> and a p-value of <p-value>.
ScorBelgium <- cor(Community_DID_Belgium$Year_simple, Community_DID_Belgium$Value, method = "spearman")
ScortestBelgium <- cor.test(Community_DID_Belgium$Year_simple, Community_DID_Belgium$Value, method = "spearman")
ScorBelgium
ScortestBelgium

# rho 
# -0.9300699         p-value < 2.2e-16


# correlation analysis - Europe

# make dataset
Community_DID_Europe <- Community_DID %>%
  filter(Region == "Europe")

#Normality check
qqnorm(Community_DID_Europe$Year_simple, pch =1, frame = FALSE)
qqline(Community_DID_Europe$Year_simple, col = "steelblue", lwd = 2)

qqnorm(Community_DID_Europe$Value, pch =1, frame = FALSE)
qqline(Community_DID_Europe$Value, col = "steelblue", lwd = 2)


qqPlot(Community_DID_Europe$Value)   # two outliers

#Shapiro test: if the  p-values are greater than the significance level 0.05 => the distribution of the data are not significantly different from normal distribution
# Shapiro-Wilk normality test for Year_simple
shapiro.test(Community_DID_Europe$Year_simple) #p = 0,9 --> can treat as normally distributed

# Shapiro-Wilk normality test for Value
shapiro.test(Community_DID_Europe$Value) # p = 0.0006368 --> can't treat as normally distributed

# #Pearson (parametric test, assumes linearity and normality - does not apply here) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <cor> and a p-value of <p-value>.
# PcorEurope <- cor(Community_DID_Europe$Year_simple, Community_DID_Europe$Value, method = "pearson")
# PcortestEurope <- cor.test(Community_DID_Europe$Year_simple, Community_DID_Europe$Value, method = "pearson")
# PcorEurope
# PcortestEurope



#Kendall (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <tau> and a p-value of <p-value>.
KcorEurope <- cor(Community_DID_Europe$Year_simple, Community_DID_Europe$Value, method = "kendall")
KcortestEurope <- cor.test(Community_DID_Europe$Year_simple, Community_DID_Europe$Value, method = "kendall")
KcorEurope
KcortestEurope

# tau 
# -0.32062     p-value = 0.1489

#Spearman (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <rho> and a p-value of <p-value>.
ScorEurope <- cor(Community_DID_Europe$Year_simple, Community_DID_Europe$Value, method = "spearman")
ScortestEurope <- cor.test(Community_DID_Europe$Year_simple, Community_DID_Europe$Value, method = "spearman")
ScorEurope
ScortestEurope

# rho 
# -0.5253949 p-value = 0.0794


# no signif trend change over time --> is it oscillating or equals (for sign on figure)
(max(Community_DID_Europe$Value)-min(Community_DID_Europe$Value))/mean(Community_DID_Europe$Value) > 0.25
#[1] FALSE ---> not oscillating


# Community DID graph ------------
(Community_DID_graph <- Community_DID %>%
   mutate(label = case_when(
     Region == "Belgium" ~ paste(downward_arrow,"***",sep=" "),
     Region == "Europe" ~ equals)) %>%
   mutate(NAP_Target = 11.86)%>%
   mutate(ycoord = case_when(
     Region == "Belgium" ~ 20,
     Region == "Europe" ~ 15)) %>%
   ggplot()+
   geom_bar(aes(x= Year, y = Value, fill = Region),position = "dodge", stat = "identity") +
   geom_hline(aes(yintercept=NAP_Target), colour = "red", linetype = "dashed")+
   scale_fill_manual(values = BELMAP_location_colourscheme)+
   scale_colour_manual(values = BELMAP_location_colourscheme)+
   scale_x_continuous(limits = c(2010,2024), breaks = seq(2011,2022,1))+
   geom_richtext(aes(x = 2023.5, y = ycoord, 
                     label = label, 
                     fill = Region), stat = "unique", 
                 colour= "white", show.legend = FALSE)+
   labs(y= "DDD/1000 inhabitants/day")&
   moiras_graph_theme()
)




#4. Community Ratio broad: narrow -----------------------------------------------------

Community_ratio_broad_narrow <- Human_AMC_Data %>%
  filter(Sector == "Ambulant",
         Measurement == "Ratio of broad spectrum antibiotics")


# correlation analysis - Belgium
#make dataset
Community_ratio_broad_narrow_Belgium <- Community_ratio_broad_narrow %>%
  filter(Region == "Belgium")

#Normality check
qqnorm(Community_ratio_broad_narrow_Belgium$Year_simple, pch =1, frame = FALSE)
qqline(Community_ratio_broad_narrow_Belgium$Year_simple, col = "steelblue", lwd = 2)

qqnorm(Community_ratio_broad_narrow_Belgium$Value, pch =1, frame = FALSE)
qqline(Community_ratio_broad_narrow_Belgium$Value, col = "steelblue", lwd = 2)

qqPlot(Community_ratio_broad_narrow_Belgium$Value)   

#Shapiro test: if the  p-values are greater than the significance level 0.05 => the distribution of the data are not significantly different from normal distribution
# Shapiro-Wilk normality test for Year_simple
shapiro.test(Community_ratio_broad_narrow_Belgium$Year_simple) #p = 0,9 --> can treat as normally distributed

# Shapiro-Wilk normality test for Value
shapiro.test(Community_ratio_broad_narrow_Belgium$Value) #  p-value = 0.7627 --> can treat as normally distributed

#Pearson (parametric test, assumes linearity and normality) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <cor> and a p-value of <p-value>.
PcorBelgium <- cor(Community_ratio_broad_narrow_Belgium$Year_simple, Community_ratio_broad_narrow_Belgium$Value, method = "pearson")
PcortestBelgium <- cor.test(Community_ratio_broad_narrow_Belgium$Year_simple, Community_ratio_broad_narrow_Belgium$Value, method = "pearson")
PcorBelgium
PcortestBelgium

# cor 
# -0.8297688    p-value = 0.0008396

#Kendall (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <tau> and a p-value of <p-value>.
KcorBelgium <- cor(Community_ratio_broad_narrow_Belgium$Year_simple, Community_ratio_broad_narrow_Belgium$Value, method = "kendall")
KcortestBelgium <- cor.test(Community_ratio_broad_narrow_Belgium$Year_simple, Community_ratio_broad_narrow_Belgium$Value, method = "kendall")
KcorBelgium
KcortestBelgium

# tau 
# -0.6870429    p-value = 0.001981

#Spearman (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <rho> and a p-value of <p-value>.
ScorBelgium <- cor(Community_ratio_broad_narrow_Belgium$Year_simple, Community_ratio_broad_narrow_Belgium$Value, method = "spearman")
ScortestBelgium <- cor.test(Community_ratio_broad_narrow_Belgium$Year_simple, Community_ratio_broad_narrow_Belgium$Value, method = "spearman")
ScorBelgium
ScortestBelgium

# rho 
# -0.7986002         p-value = 0.00184


# correlation analysis - Europe

# make dataset
Community_ratio_broad_narrow_Europe <- Community_ratio_broad_narrow %>%
  filter(Region == "Europe")

#Normality check
qqnorm(Community_ratio_broad_narrow_Europe$Year_simple, pch =1, frame = FALSE)
qqline(Community_ratio_broad_narrow_Europe$Year_simple, col = "steelblue", lwd = 2)

qqnorm(Community_ratio_broad_narrow_Europe$Value, pch =1, frame = FALSE)
qqline(Community_ratio_broad_narrow_Europe$Value, col = "steelblue", lwd = 2)


qqPlot(Community_ratio_broad_narrow_Europe$Value)   # one outlier

#Shapiro test: if the  p-values are greater than the significance level 0.05 => the distribution of the data are not significantly different from normal distribution
# Shapiro-Wilk normality test for Year_simple
shapiro.test(Community_ratio_broad_narrow_Europe$Year_simple) #p = 0,9 --> can treat as normally distributed

# Shapiro-Wilk normality test for Value
shapiro.test(Community_ratio_broad_narrow_Europe$Value) # p = 0.006153 --> can't treat as normally distributed

# #Pearson (parametric test, assumes linearity and normality - does not apply here) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <cor> and a p-value of <p-value>.
# PcorEurope <- cor(Community_ratio_broad_narrow_Europe$Year_simple, Community_ratio_broad_narrow_Europe$Value, method = "pearson")
# PcortestEurope <- cor.test(Community_ratio_broad_narrow_Europe$Year_simple, Community_ratio_broad_narrow_Europe$Value, method = "pearson")
# PcorEurope
# PcortestEurope



#Kendall (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <tau> and a p-value of <p-value>.
KcorEurope <- cor(Community_ratio_broad_narrow_Europe$Year_simple, Community_ratio_broad_narrow_Europe$Value, method = "kendall")
KcortestEurope <- cor.test(Community_ratio_broad_narrow_Europe$Year_simple, Community_ratio_broad_narrow_Europe$Value, method = "kendall")
KcorEurope
KcortestEurope

# tau 
# 0.7575758     p-value = 0.00024

#Spearman (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <rho> and a p-value of <p-value>.
ScorEurope <- cor(Community_ratio_broad_narrow_Europe$Year_simple, Community_ratio_broad_narrow_Europe$Value, method = "spearman")
ScortestEurope <- cor.test(Community_ratio_broad_narrow_Europe$Year_simple, Community_ratio_broad_narrow_Europe$Value, method = "spearman")
ScorEurope
ScortestEurope

# rho 
# 0.8741259 p-value = 0.0003089


# Community broad to narrow ratio graph ------------
(Community_ratio_broad_narrow_graph <- Community_ratio_broad_narrow %>%
   mutate(label = case_when(
     Region == "Belgium" ~ paste(downward_arrow,"***",sep=" "),
     Region == "Europe" ~ paste(upward_arrow,"***",sep=" "))) %>%
   #mutate(NAP_Target = 11.86)%>%
   mutate(ycoord = case_when(
     Region == "Belgium" ~ 2,
     Region == "Europe" ~ 2.5)) %>%
   ggplot()+
   geom_bar(aes(x= Year, y = Value, fill = Region),position = "dodge", stat = "identity") +
   # geom_hline(aes(yintercept=NAP_Target), colour = "red", linetype = "dashed")+
   scale_fill_manual(values = BELMAP_location_colourscheme)+
   scale_colour_manual(values = BELMAP_location_colourscheme)+
   scale_x_continuous(limits = c(2010,2024), breaks = seq(2011,2022,1))+
   geom_richtext(aes(x = 2023.5, y = ycoord, 
                     label = label, 
                     fill = Region), stat = "unique", 
                 colour= "white", show.legend = FALSE)+
   labs(y= "DDD/1000 inhabitants/day")&
   moiras_graph_theme()
)


#5. Hospital DDD per patient days-----------------------------------------------------


Hospital_DDD_patient_days <- Human_AMC_Data %>%
  filter(grepl("Hospitals",Sector),
         Unit == "DDD/1000 patients/day")


# correlation analysis - Acute
#make dataset
Hospital_DDD_patient_days_Acute <- Hospital_DDD_patient_days %>%
  filter(Sector == "Acute Hospitals")

#Normality check
qqnorm(Hospital_DDD_patient_days_Acute$Year_simple, pch =1, frame = FALSE)
qqline(Hospital_DDD_patient_days_Acute$Year_simple, col = "steelblue", lwd = 2)

qqnorm(Hospital_DDD_patient_days_Acute$Value, pch =1, frame = FALSE)
qqline(Hospital_DDD_patient_days_Acute$Value, col = "steelblue", lwd = 2)

qqPlot(Hospital_DDD_patient_days_Acute$Value)   

#Shapiro test: if the  p-values are greater than the significance level 0.05 => the distribution of the data are not significantly different from normal distribution
# Shapiro-Wilk normality test for Year_simple
shapiro.test(Hospital_DDD_patient_days_Acute$Year_simple) #p = 0,9 --> can treat as normally distributed

# Shapiro-Wilk normality test for Value
shapiro.test(Hospital_DDD_patient_days_Acute$Value) #  p-value = 0.8206 --> can treat as normally distributed

#Pearson (parametric test, assumes linearity and normality) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <cor> and a p-value of <p-value>.
PcorAcute <- cor(Hospital_DDD_patient_days_Acute$Year_simple, Hospital_DDD_patient_days_Acute$Value, method = "pearson")
PcortestAcute <- cor.test(Hospital_DDD_patient_days_Acute$Year_simple, Hospital_DDD_patient_days_Acute$Value, method = "pearson")
PcorAcute
PcortestAcute

# cor 
# 0.8911427    -value = 0.0002313

#Kendall (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <tau> and a p-value of <p-value>.
KcorAcute <- cor(Hospital_DDD_patient_days_Acute$Year_simple, Hospital_DDD_patient_days_Acute$Value, method = "kendall")
KcortestAcute <- cor.test(Hospital_DDD_patient_days_Acute$Year_simple, Hospital_DDD_patient_days_Acute$Value, method = "kendall")
KcorAcute
KcortestAcute

# tau 
# 0.7454545   p-value = 0.0007595

#Spearman (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <rho> and a p-value of <p-value>.
ScorAcute <- cor(Hospital_DDD_patient_days_Acute$Year_simple, Hospital_DDD_patient_days_Acute$Value, method = "spearman")
ScortestAcute <- cor.test(Hospital_DDD_patient_days_Acute$Year_simple, Hospital_DDD_patient_days_Acute$Value, method = "spearman")
ScorAcute
ScortestAcute

# rho 
# 0.9090909          p-value = 5.554e-05



# correlation analysis - Chronic
#make dataset
Hospital_DDD_patient_days_Chronic <- Hospital_DDD_patient_days %>%
  filter(Sector == "Chronic Hospitals")

#Normality check
qqnorm(Hospital_DDD_patient_days_Chronic$Year_simple, pch =1, frame = FALSE)
qqline(Hospital_DDD_patient_days_Chronic$Year_simple, col = "steelblue", lwd = 2)

qqnorm(Hospital_DDD_patient_days_Chronic$Value, pch =1, frame = FALSE)
qqline(Hospital_DDD_patient_days_Chronic$Value, col = "steelblue", lwd = 2)

qqPlot(Hospital_DDD_patient_days_Chronic$Value)   

#Shapiro test: if the  p-values are greater than the significance level 0.05 => the distribution of the data are not significantly different from normal distribution
# Shapiro-Wilk normality test for Year_simple
shapiro.test(Hospital_DDD_patient_days_Chronic$Year_simple) #p = 0,9 --> can treat as normally distributed

# Shapiro-Wilk normality test for Value
shapiro.test(Hospital_DDD_patient_days_Chronic$Value) # p-value = 0.9647 --> can treat as normally distributed

#Pearson (parametric test, assumes linearity and normality) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <cor> and a p-value of <p-value>.
PcorChronic <- cor(Hospital_DDD_patient_days_Chronic$Year_simple, Hospital_DDD_patient_days_Chronic$Value, method = "pearson")
PcortestChronic <- cor.test(Hospital_DDD_patient_days_Chronic$Year_simple, Hospital_DDD_patient_days_Chronic$Value, method = "pearson")
PcorChronic
PcortestChronic

# cor 
# -0.2922573     p-value = 0.4824

#Kendall (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <tau> and a p-value of <p-value>.
KcorChronic <- cor(Hospital_DDD_patient_days_Chronic$Year_simple, Hospital_DDD_patient_days_Chronic$Value, method = "kendall")
KcortestChronic <- cor.test(Hospital_DDD_patient_days_Chronic$Year_simple, Hospital_DDD_patient_days_Chronic$Value, method = "kendall")
KcorChronic
KcortestChronic

# tau 
# -0.2142857    p-value = 0.5484

#Spearman (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <rho> and a p-value of <p-value>.
ScorChronic <- cor(Hospital_DDD_patient_days_Chronic$Year_simple, Hospital_DDD_patient_days_Chronic$Value, method = "spearman", na.rm=TRUE)
ScortestChronic <- cor.test(Hospital_DDD_patient_days_Chronic$Year_simple, Hospital_DDD_patient_days_Chronic$Value, method = "spearman", na.rm=TRUE)
ScorChronic
ScortestChronic

# rho 
# -0.2619048         p-value = 0.5364


# no significant trend - test oscillating vs equals icon:
(max(Hospital_DDD_patient_days_Chronic$Value,na.rm=TRUE)-min(Hospital_DDD_patient_days_Chronic$Value,na.rm=TRUE))/mean(Hospital_DDD_patient_days_Chronic$Value,na.rm=TRUE) > 0.25
#[1] FALSE

# correlation analysis - Psychiatric
#make dataset
Hospital_DDD_patient_days_Psychiatric <- Hospital_DDD_patient_days %>%
  filter(Sector == "Psychiatric Hospitals")

#Normality check
qqnorm(Hospital_DDD_patient_days_Psychiatric$Year_simple, pch =1, frame = FALSE)
qqline(Hospital_DDD_patient_days_Psychiatric$Year_simple, col = "steelblue", lwd = 2)

qqnorm(Hospital_DDD_patient_days_Psychiatric$Value, pch =1, frame = FALSE)
qqline(Hospital_DDD_patient_days_Psychiatric$Value, col = "steelblue", lwd = 2)

qqPlot(Hospital_DDD_patient_days_Psychiatric$Value)   

#Shapiro test: if the  p-values are greater than the significance level 0.05 => the distribution of the data are not significantly different from normal distribution
# Shapiro-Wilk normality test for Year_simple
shapiro.test(Hospital_DDD_patient_days_Psychiatric$Year_simple) #p = 0,9 --> can treat as normally distributed

# Shapiro-Wilk normality test for Value
shapiro.test(Hospital_DDD_patient_days_Psychiatric$Value) #  p-value = 0.7244 --> can treat as normally distributed

#Pearson (parametric test, assumes linearity and normality) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <cor> and a p-value of <p-value>.
PcorPsychiatric <- cor(Hospital_DDD_patient_days_Psychiatric$Year_simple, Hospital_DDD_patient_days_Psychiatric$Value, method = "pearson")
PcortestPsychiatric <- cor.test(Hospital_DDD_patient_days_Psychiatric$Year_simple, Hospital_DDD_patient_days_Psychiatric$Value, method = "pearson")
PcorPsychiatric
PcortestPsychiatric

# cor 
# -0.3031941      p-value = 0.3648

#Kendall (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <tau> and a p-value of <p-value>.
KcorPsychiatric <- cor(Hospital_DDD_patient_days_Psychiatric$Year_simple, Hospital_DDD_patient_days_Psychiatric$Value, method = "kendall")
KcortestPsychiatric <- cor.test(Hospital_DDD_patient_days_Psychiatric$Year_simple, Hospital_DDD_patient_days_Psychiatric$Value, method = "kendall")
KcorPsychiatric
KcortestPsychiatric

# tau 
# -0.3090909     p-value = 0.2183

#Spearman (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <rho> and a p-value of <p-value>.
ScorPsychiatric <- cor(Hospital_DDD_patient_days_Psychiatric$Year_simple, Hospital_DDD_patient_days_Psychiatric$Value, method = "spearman")
ScortestPsychiatric <- cor.test(Hospital_DDD_patient_days_Psychiatric$Year_simple, Hospital_DDD_patient_days_Psychiatric$Value, method = "spearman")
ScorPsychiatric
ScortestPsychiatric

# rho 
# -0.2909091          p-value = 0.3864

# no significant trend - test oscillating vs equals icon:
(max(Hospital_DDD_patient_days_Psychiatric$Value,na.rm=TRUE)-min(Hospital_DDD_patient_days_Psychiatric$Value,na.rm=TRUE))/mean(Hospital_DDD_patient_days_Psychiatric$Value,na.rm=TRUE) > 0.25
#[1] TRUE





# Hospital DDD per patient days ratio graph ------------
(Hospital_DDD_patient_days_graph <- Hospital_DDD_patient_days %>%
   mutate(label = case_when(
     Sector == "Acute Hospitals" ~ paste(upward_arrow,"***",sep=" "),
     Sector == "Chronic Hospitals" ~ equals,
     Sector == "Psychiatric Hospitals" ~ oscillate
   )) %>%
   #mutate(NAP_Target = 11.86)%>%
   mutate(ycoord = case_when(
     Sector == "Acute Hospitals" ~ 450,
     Sector == "Chronic Hospitals" ~ 150,
     Sector == "Psychiatric Hospitals" ~ 50)) %>%
   ggplot()+
   geom_bar(aes(x= Year, y = Value, fill = Sector),position = "dodge", stat = "identity") +
   # geom_hline(aes(yintercept=NAP_Target), colour = "red", linetype = "dashed")+
   scale_fill_manual(values = BELMAP_colourscheme)+
   scale_colour_manual(values = BELMAP_colourscheme)+
   scale_x_continuous(limits = c(2010,2024), breaks = seq(2011,2022,1))+
   geom_richtext(aes(x = 2023.5, y = ycoord, 
                     label = label, 
                     fill = Sector), stat = "unique", 
                 colour= "white", show.legend = FALSE)+
   labs(y= "DDD/1000 inhabitants/day")&
   moiras_graph_theme()
)



#6. Hospital Ratio broad:narrow-----------------------------------------------------



Hospital_ratio_broad_narrow <- Human_AMC_Data %>%
  filter(Sector == "Hospital",
         grepl("proportion of glycopeptides, third- and fourth-generation cephalosporins",Measurement))

# correlation analysis - Belgium
#make dataset
Hospital_ratio_broad_narrow_Belgium <- Hospital_ratio_broad_narrow %>%
  filter(Region == "Belgium")

#Normality check
qqnorm(Hospital_ratio_broad_narrow_Belgium$Year_simple, pch =1, frame = FALSE)
qqline(Hospital_ratio_broad_narrow_Belgium$Year_simple, col = "steelblue", lwd = 2)

qqnorm(Hospital_ratio_broad_narrow_Belgium$Value, pch =1, frame = FALSE)
qqline(Hospital_ratio_broad_narrow_Belgium$Value, col = "steelblue", lwd = 2)

qqPlot(Hospital_ratio_broad_narrow_Belgium$Value)   

#Shapiro test: if the  p-values are greater than the significance level 0.05 => the distribution of the data are not significantly different from normal distribution
# Shapiro-Wilk normality test for Year_simple
shapiro.test(Hospital_ratio_broad_narrow_Belgium$Year_simple) #p = 0,9 --> can treat as normally distributed

# Shapiro-Wilk normality test for Value
shapiro.test(Hospital_ratio_broad_narrow_Belgium$Value) #   p-value = 0.2093 --> can treat as normally distributed

#Pearson (parametric test, assumes linearity and normality) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <cor> and a p-value of <p-value>.
PcorBelgium <- cor(Hospital_ratio_broad_narrow_Belgium$Year_simple, Hospital_ratio_broad_narrow_Belgium$Value, method = "pearson")
PcortestBelgium <- cor.test(Hospital_ratio_broad_narrow_Belgium$Year_simple, Hospital_ratio_broad_narrow_Belgium$Value, method = "pearson")
PcorBelgium
PcortestBelgium

# cor 
# -0.2521143     p-value = 0.4292

#Kendall (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <tau> and a p-value of <p-value>.
KcorBelgium <- cor(Hospital_ratio_broad_narrow_Belgium$Year_simple, Hospital_ratio_broad_narrow_Belgium$Value, method = "kendall")
KcortestBelgium <- cor.test(Hospital_ratio_broad_narrow_Belgium$Year_simple, Hospital_ratio_broad_narrow_Belgium$Value, method = "kendall")
KcorBelgium
KcortestBelgium

# tau 
# -0.1984791     p-value = 0.3716 

#Spearman (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <rho> and a p-value of <p-value>.
ScorBelgium <- cor(Hospital_ratio_broad_narrow_Belgium$Year_simple, Hospital_ratio_broad_narrow_Belgium$Value, method = "spearman")
ScortestBelgium <- cor.test(Hospital_ratio_broad_narrow_Belgium$Year_simple, Hospital_ratio_broad_narrow_Belgium$Value, method = "spearman")
ScorBelgium
ScortestBelgium

# rho 
# -0.2697027          p-value = 0.3966



# no signif trend change over time --> is it oscillating or equals (for sign on figure)
(max(Hospital_ratio_broad_narrow_Belgium$Value)-min(Hospital_ratio_broad_narrow_Belgium$Value))/mean(Hospital_ratio_broad_narrow_Belgium$Value) > 0.25
#[1] FALSE ---> not oscillating


# correlation analysis - Europe

# make dataset
Hospital_ratio_broad_narrow_Europe <- Hospital_ratio_broad_narrow %>%
  filter(Region == "Europe")

#Normality check
qqnorm(Hospital_ratio_broad_narrow_Europe$Year_simple, pch =1, frame = FALSE)
qqline(Hospital_ratio_broad_narrow_Europe$Year_simple, col = "steelblue", lwd = 2)

qqnorm(Hospital_ratio_broad_narrow_Europe$Value, pch =1, frame = FALSE)
qqline(Hospital_ratio_broad_narrow_Europe$Value, col = "steelblue", lwd = 2)


qqPlot(Hospital_ratio_broad_narrow_Europe$Value) 

#Shapiro test: if the  p-values are greater than the significance level 0.05 => the distribution of the data are not significantly different from normal distribution
# Shapiro-Wilk normality test for Year_simple
shapiro.test(Hospital_ratio_broad_narrow_Europe$Year_simple) #p = 0,9 --> can treat as normally distributed

# Shapiro-Wilk normality test for Value
shapiro.test(Hospital_ratio_broad_narrow_Europe$Value) # p = 0.09441 --> can treat as normally distributed

#Pearson (parametric test, assumes linearity and normality - does not apply here) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <cor> and a p-value of <p-value>.
PcorEurope <- cor(Hospital_ratio_broad_narrow_Europe$Year_simple, Hospital_ratio_broad_narrow_Europe$Value, method = "pearson")
PcortestEurope <- cor.test(Hospital_ratio_broad_narrow_Europe$Year_simple, Hospital_ratio_broad_narrow_Europe$Value, method = "pearson")
PcorEurope
PcortestEurope

#cor
#0.705443  p-value = 0.01038


#Kendall (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <tau> and a p-value of <p-value>.
KcorEurope <- cor(Hospital_ratio_broad_narrow_Europe$Year_simple, Hospital_ratio_broad_narrow_Europe$Value, method = "kendall")
KcortestEurope <- cor.test(Hospital_ratio_broad_narrow_Europe$Year_simple, Hospital_ratio_broad_narrow_Europe$Value, method = "kendall")
KcorEurope
KcortestEurope

# tau 
# 0.4732962      p-value = 0.03311

#Spearman (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <rho> and a p-value of <p-value>.
ScorEurope <- cor(Hospital_ratio_broad_narrow_Europe$Year_simple, Hospital_ratio_broad_narrow_Europe$Value, method = "spearman")
ScortestEurope <- cor.test(Hospital_ratio_broad_narrow_Europe$Year_simple, Hospital_ratio_broad_narrow_Europe$Value, method = "spearman")
ScorEurope
ScortestEurope

# rho 
# 0.5954475  p-value = 0.04107



# Hospital broad to narrow ratio graph ------------
(Hospital_ratio_broad_narrow_graph <- Hospital_ratio_broad_narrow %>%
   mutate(label = case_when(
     Region == "Belgium" ~ equals,
     Region == "Europe" ~ paste(upward_arrow,"*",sep=" "))) %>%
   #mutate(NAP_Target = 11.86)%>%
   mutate(ycoord = case_when(
     Region == "Belgium" ~ 0.3,
     Region == "Europe" ~ 0.4)) %>%
   ggplot()+
   geom_bar(aes(x= Year, y = Value, fill = Region),position = "dodge", stat = "identity") +
   # geom_hline(aes(yintercept=NAP_Target), colour = "red", linetype = "dashed")+
   scale_fill_manual(values = BELMAP_location_colourscheme)+
   scale_colour_manual(values = BELMAP_location_colourscheme)+
   scale_x_continuous(limits = c(2010,2024), breaks = seq(2011,2022,1))+
   geom_richtext(aes(x = 2023.5, y = ycoord,
                     label = label,
                     fill = Region), stat = "unique",
                 colour= "white", show.legend = FALSE)+
   labs(y= "DDD/1000 inhabitants/day")&
   moiras_graph_theme()
)





#7. Hospital proportional consumption ATC groups-----------------------------------------------------

# 
# # make dataset
# Hospital_proportion_consumption <- Human_AMC_Data %>%
#   filter(grepl("^Proportion",Measurement))
# 
# 
# start_data_frame <- tibble(
#   "Year" = "",
#   "Measurement" = "",
#   "ATC_Group" = "",
#   "Value"= "",
#   "Unit"= "",
#   "Sector" = "",
#   "Region" = "",
#   "Source" = "",
#   "Year_simple"= "",
#   "icon" = "",
#   "signif" = ""
# )
# 
# 
# 
# write_csv(start_data_frame,file = "3_Human_AMC_analysis_outcomes.csv")
# 
# 
# 
# # list proportions indicators
# 
# proportions_categories<- unique(Hospital_proportion_consumption$Measurement)
# 
# #mkae output test results file
# capture.output(print("Analysis of Human AMC proportions data"),file = "3_Human_AMC_correlations/proportions_data_analysis.txt")  
# 
# 
# #make list for overdispersion plots
# qq_plot_list<- list()
# #open for loop proportions ---------------------------
# #i =  "Total proportions \\(mg\\/kg biomass\\)"
# for(i in proportions_categories){
#   
#   # 1. check assumptions for Pearsons
#   dataset_proportions<- Hospital_proportion_consumption %>%
#     filter(grepl(i,Measurement,  fixed = TRUE))
#   
#   capture.output(print(paste("Analysis for",i,sep = " : ")),file = "3_Human_AMC_correlations/proportions_data_analysis.txt", append = TRUE)
#   
#   #save qq plot
#   
#   #Normality check
#   # qqnorm(dataset_proportions$Year_simple, pch =1, frame = FALSE)
#   # qqline(dataset_proportions$Year_simple, col = "steelblue", lwd = 2)
#   #
#   # qqnorm(dataset_proportions$Value, pch =1, frame = FALSE)
#   # qqline(dataset_proportions$Value, col = "steelblue", lwd = 2)
#   
#   
#   qqPlot(dataset_proportions$Value)
#   qqPlot_value<- recordPlot()
#   
#   qqPlot(dataset_proportions$Year)
#   qqPlot_year<- recordPlot()
#   
#   #assign to local environment
#   assign(paste(i,"qqplot_value"),print(qqPlot_value))
#   assign(paste(i,"qqplot_year"),print(qqPlot_year))
#   
#   #Shapiro test: if the  p-values are greater than the significance level 0.05 => the distribution of the data are not significantly different from normal distribution
#   # Shapiro-Wilk normality test for Year_simple
#   
#   Shapiro_year <- shapiro.test(dataset_proportions$Year_simple) #p-value = 0.8698
#   capture.output(print("Testing Normality of year"),file = "3_Human_AMC_correlations/proportions_data_analysis.txt", append = TRUE)
#   capture.output(print(Shapiro_year),file = "3_Human_AMC_correlations/proportions_data_analysis.txt", append = TRUE)
#   
#   
#   # Shapiro-Wilk normality test for Value
#   Shapiro_value <- shapiro.test(dataset_proportions$Value) #p-value = 0.4784
#   
#   capture.output(print("Testing Normality of Values"),file = "3_Human_AMC_correlations/proportions_data_analysis.txt", append = TRUE)
#   capture.output(print(Shapiro_value),file = "3_Human_AMC_correlations/proportions_data_analysis.txt", append = TRUE)
#   
#   
#   # 2. calculate correlation tests
#   
#   #Pearson (parametric test, assumes linearity and normality) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <cor> and a p-value of <p-value>.
#   Pcor_Value <- cor(dataset_proportions$Year_simple, dataset_proportions$Value, method = "pearson")
#   Pcortest_Value <- cor.test(dataset_proportions$Year_simple, dataset_proportions$Value, method = "pearson")
#   Pcor_Value
#   #        cor -0.9856587; p-value = 2.905e-08
#   Pcortest_Value
#   #        cor -0.9856587; p-value = 2.905e-08
#   
#   capture.output(print("#Pearson (parametric test, assumes linearity and normality) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <cor> and a p-value of <p-value>.
# "),file = "3_Human_AMC_correlations/proportions_data_analysis.txt", append = TRUE)
#   capture.output(print(Pcortest_Value),file = "3_Human_AMC_correlations/proportions_data_analysis.txt", append = TRUE)
#   
#   
#   #Kendall (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <tau> and a p-value of <p-value>.
#   Kcor_Value <- cor(dataset_proportions$Year_simple, dataset_proportions$Value, method = "kendall")
#   Kcortest_Value <- cor.test(dataset_proportions$Year_simple, dataset_proportions$Value, method = "kendall")
#   Kcor_Value
#   # [1] -0.9272727
#   
#   Kcortest_Value
#   #  tau -0.9272727;  p-value = 3.257e-06
#   
#   capture.output(print("#Kendall (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <tau> and a p-value of <p-value>.
# "),file = "3_Human_AMC_correlations/proportions_data_analysis.txt",  append = TRUE)
#   capture.output(print(Kcortest_Value),file = "3_Human_AMC_correlations/proportions_data_analysis.txt", append = TRUE)
#   
#   
#   
#   #Spearman (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <rho> and a p-value of <p-value>.
#   Scor_Value <- cor(dataset_proportions$Year_simple, dataset_proportions$Value, method = "spearman")
#   Scortest_Value <- cor.test(dataset_proportions$Year_simple, dataset_proportions$Value, method = "spearman")
#   Scor_Value
#   #[1] -0.9818182
#   
#   Scortest_Value
#   
#   capture.output(print("#Spearman (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <rho> and a p-value of <p-value>.
# "),file = "3_Human_AMC_correlations/proportions_data_analysis.txt",  append = TRUE)
#   capture.output(print(Scortest_Value),file = "3_Human_AMC_correlations/proportions_data_analysis.txt", append = TRUE)
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
#     if((max(dataset_proportions$Value)-min(dataset_proportions$Value))/mean(dataset_proportions$Value) > 0.25){
#       model_icon <- "oscilate"
#       signif = ""
#     }else if((max(dataset_proportions$Value)-min(dataset_proportions$Value))/mean(dataset_proportions$Value) < 0.25){
#       model_icon <- "equals"
#       signif <- ""
#     }
#     
#   }
#   
#   # add icon and signif
#   
#   
#   dataset_proportions_analysis <- dataset_proportions %>%
#     mutate(icon = model_icon) %>%
#     mutate(signif = signif)
#   
#    # Save output to dataframe for graphs ---------------------
#   
#   write_csv(dataset_proportions_analysis,file = "3_Human_AMC_analysis_outcomes.csv", append = TRUE)

#close for loop proportions ---------------------------

}

Hospital_proportions_graph_data <-read.csv("3_Human_AMC_analysis_outcomes.csv", header = TRUE, sep = ";")

# Hospital proportions ratio graph ------------
(Hospital_proportions_graph <- Hospital_proportions_graph_data %>%
   mutate(ycoord = 0.1,
          Antibiotic = case_when(
            grepl("glycopeptides",Measurement) ~ "Glycopeptides",
            grepl("CSP",Measurement) ~ "3rd and 4th generation\ncephalosporins",
            grepl("FQ",Measurement) ~ "Fluoroquinolones",
            grepl("Polymixins",Measurement) ~ "Polymixins",
            grepl("piperacillin",Measurement) ~ "Piperacillin + INH",
            grepl("Linezolid",Measurement) ~ "Linezolid",
            grepl("carbapenems",Measurement) ~ "Carbapenems"
          )) %>%
   mutate(label_icon = case_when(
     icon == "upward_arrow" ~ upward_arrow,
     icon == "downward_arrow" ~ downward_arrow,
     icon == "equals" ~ equals,
     icon == "oscilate" ~ oscillate,
     .default = ""
   )) %>%
   mutate(label = if_else((is.na(signif) | signif == ""), label_icon, paste(label_icon,signif,sep=" "))) %>%
   ggplot()+
   geom_bar(aes(x= Year, y = Value),position = "dodge", stat = "identity", fill = "#3288bd") +
   # geom_hline(aes(yintercept=NAP_Target), colour = "red", linetype = "dashed")+
   facet_wrap(~Antibiotic, nrow=2) +
   #  scale_fill_manual(values = BELMAP_Human_AntiB_colourscheme )+
   #  scale_colour_manual(values = BELMAP_Human_AntiB_colourscheme )+
   scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.25))+
   scale_x_continuous(limits = c(2010,2024), breaks = seq(2011,2022,1))+
   geom_richtext(aes(x = 2023.5, y = ycoord,
                     label = label), stat = "unique", fill = "#3288bd",
                 colour= "white", show.legend = FALSE)+
   labs(y= "Proportion of hospital antimicrobial consumption")&
   moiras_graph_theme()
)












