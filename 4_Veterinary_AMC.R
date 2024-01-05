#-------------3. Veterinary AMC ------------------------------------

#This script outlines the analysis for the BELMAP 2023 Chapter 4 - AMC in the Veterinary Sector

#Contents-----
#1. Load libraries
#2. Load data
#3. Sales analysis
#4. Sales data graphs 
#5. BD100 Use analysis
#6. BD100 Use graphs


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

#2. Load data--------
Vet_AMC_Data <- read.csv2("Data/1.AMC/Veterinary_AMC/Veterinary_AMC_2023.csv",
                          sep = ";", header = TRUE)%>%
  mutate(Year_simple= Year - 2011) %>%
  dplyr::select(-X)

Vet_AMC_sales <-  Vet_AMC_Data %>%
  filter(grepl("mg/kg biomass", Unit))

Vet_AMC_use <-  Vet_AMC_Data %>%
  filter(grepl("Median BD100", Unit))

#3. Sales analysis--------------------------------------------------
#Correlation analysis - sales data
#---------------------


start_data_frame <- tibble(
  "Year" = "",
  "Measurement" = "",
  "Value"= "",
  "Unit"= "",
  "Year_simple"= "",
  "Target" = "",
  "icon" = "",
  "signif" = ""
)

#write_csv(start_data_frame,file = "4_vet_AMC_analysis_outcomes.csv")

# list targets -----------------------------------------

NAP_targets <- tibble(
  Measurement = c("Total sales (mg/kg biomass)", 
                  "Sales fluoro(quinolones)+ 3rd- and 4th-generation cephalosporins (mg/kg biomass)",
                  "Sales premixes (mg/kg biomass)",
                  "Sales polymyxins (mg/kg biomass)"),
  Target = c(51.28, 0.5025,14.05,1)
  
)

# total_target= 51.28 #35% of 2011 value              # 2020 target was: 73.25 #50% 2011 value  #
# premix_target= 14.05 #50% 2011 value
# red_target = 0.5025 # 25% 2011 value - quinolones and 3rd/4th ceph = red category
# colistin (polymixin = 1)


# list sales indicators

Sales_categories<- unique(Vet_AMC_sales$Measurement)

#mkae output test results file
capture.output(print("Analysis of Veterinary Sales data"),file = "4_Vet_AMC_correlations/Sales_data_analysis.txt")  


#make list for overdispersion plots
qq_plot_list<- list()
#open for loop sales ---------------------------
#i =  "Total sales \\(mg\\/kg biomass\\)"
for(i in Sales_categories){
  
  # 1. check assumptions for Pearsons
  dataset_sales<- Vet_AMC_sales %>%
    filter(grepl(i,Measurement,  fixed = TRUE))
  
  capture.output(print(paste("Analysis for",i,sep = " : ")),file = "4_Vet_AMC_correlations/Sales_data_analysis.txt", append = TRUE)
  
  #save qq plot
  
  #Normality check
  # qqnorm(dataset_sales$Year_simple, pch =1, frame = FALSE)
  # qqline(dataset_sales$Year_simple, col = "steelblue", lwd = 2)
  #
  # qqnorm(dataset_sales$Value, pch =1, frame = FALSE)
  # qqline(dataset_sales$Value, col = "steelblue", lwd = 2)
  
  
  qqPlot(dataset_sales$Value)
  qqPlot_value<- recordPlot()
  
  qqPlot(dataset_sales$Year)
  qqPlot_year<- recordPlot()
  
  #assign to local environment
  assign(paste(i,"qqplot_value"),print(qqPlot_value))
  assign(paste(i,"qqplot_year"),print(qqPlot_year))
  
  #Shapiro test: if the  p-values are greater than the significance level 0.05 => the distribution of the data are not significantly different from normal distribution
  # Shapiro-Wilk normality test for Year_simple
  
  Shapiro_year <- shapiro.test(dataset_sales$Year_simple) #p-value = 0.8698
  capture.output(print("Testing Normality of year"),file = "4_Vet_AMC_correlations/Sales_data_analysis.txt", append = TRUE)
  capture.output(print(Shapiro_year),file = "4_Vet_AMC_correlations/Sales_data_analysis.txt", append = TRUE)
  
  
  # Shapiro-Wilk normality test for Value
  Shapiro_value <- shapiro.test(dataset_sales$Value) #p-value = 0.4784
  
  capture.output(print("Testing Normality of Values"),file = "4_Vet_AMC_correlations/Sales_data_analysis.txt", append = TRUE)
  capture.output(print(Shapiro_value),file = "4_Vet_AMC_correlations/Sales_data_analysis.txt", append = TRUE)
  
  
  # 2. calculate correlation tests
  
  #Pearson (parametric test, assumes linearity and normality) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <cor> and a p-value of <p-value>.
  Pcor_Value <- cor(dataset_sales$Year_simple, dataset_sales$Value, method = "pearson")
  Pcortest_Value <- cor.test(dataset_sales$Year_simple, dataset_sales$Value, method = "pearson")
  Pcor_Value
  #        cor -0.9856587; p-value = 2.905e-08
  Pcortest_Value
  #        cor -0.9856587; p-value = 2.905e-08
  
  capture.output(print("#Pearson (parametric test, assumes linearity and normality) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <cor> and a p-value of <p-value>.
"),file = "4_Vet_AMC_correlations/Sales_data_analysis.txt", append = TRUE)
  capture.output(print(Pcortest_Value),file = "4_Vet_AMC_correlations/Sales_data_analysis.txt", append = TRUE)
  
  
  #Kendall (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <tau> and a p-value of <p-value>.
  Kcor_Value <- cor(dataset_sales$Year_simple, dataset_sales$Value, method = "kendall")
  Kcortest_Value <- cor.test(dataset_sales$Year_simple, dataset_sales$Value, method = "kendall")
  Kcor_Value
  # [1] -0.9272727
  
  Kcortest_Value
  #  tau -0.9272727;  p-value = 3.257e-06
  
  capture.output(print("#Kendall (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <tau> and a p-value of <p-value>.
"),file = "4_Vet_AMC_correlations/Sales_data_analysis.txt",  append = TRUE)
  capture.output(print(Kcortest_Value),file = "4_Vet_AMC_correlations/Sales_data_analysis.txt", append = TRUE)
  
  
  
  #Spearman (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <rho> and a p-value of <p-value>.
  Scor_Value <- cor(dataset_sales$Year_simple, dataset_sales$Value, method = "spearman")
  Scortest_Value <- cor.test(dataset_sales$Year_simple, dataset_sales$Value, method = "spearman")
  Scor_Value
  #[1] -0.9818182
  
  Scortest_Value
  
  capture.output(print("#Spearman (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <rho> and a p-value of <p-value>.
"),file = "4_Vet_AMC_correlations/Sales_data_analysis.txt",  append = TRUE)
  capture.output(print(Scortest_Value),file = "4_Vet_AMC_correlations/Sales_data_analysis.txt", append = TRUE)
  
  
  # manually check that results of correlation tests are aligned
  
  
  # 3.1 select correlation and save icon: sales-------------
  #select model
  if(Shapiro_value$p.value >  0.05 ){   # if Shapiro > 0.05 - can be treated as normally distributed
    select_correlation_value <-  Pcortest_Value$estimate
    select_correlation_pvalue <-  Pcortest_Value$p.value
  }else if(Shapiro_value$p.value <  0.05){ # if Shapiro < 0.05 - not normally distributed
    select_correlation_value <-  Scortest_Value$estimate
    select_correlation_pvalue <-  Scortest_Value$p.value
  }else{
    print(paste("ERROR no model selected for",i))
  }
  
  #if signif - extract direction of model
  if(select_correlation_pvalue<0.05){
    if(select_correlation_value > 0){
      model_icon <- "upward_arrow"
      signif = case_when(
        select_correlation_pvalue > 0.01 && select_correlation_pvalue < 0.05 ~ "*",
        select_correlation_pvalue > 0.001 && select_correlation_pvalue < 0.01 ~ "**",
        select_correlation_pvalue < 0.001 ~ "***"
      )
    }else if(select_correlation_value < 0){
      model_icon = "downward_arrow"
      signif = case_when(
        select_correlation_pvalue > 0.01 && select_correlation_pvalue < 0.05 ~ "*",
        select_correlation_pvalue > 0.001 && select_correlation_pvalue < 0.01 ~ "**",
        select_correlation_pvalue < 0.001 ~ "***"
      )
    }
  } else if(select_correlation_pvalue>0.05){
    if((max(dataset_sales$Value)-min(dataset_sales$Value))/mean(dataset_sales$Value) > 0.25){
      model_icon <- "oscilate"
      signif = ""
    }else if((max(dataset_sales$Value)-min(dataset_sales$Value))/mean(dataset_sales$Value) < 0.25){
      model_icon <- "equals"
      signif <- ""
    }
  }
  
  
  # Add NAP-AMR target-------------------------
  dataset_sales_analysis <- left_join(dataset_sales,NAP_targets, by = c("Measurement")) %>%
    mutate(icon = model_icon) %>%
    mutate(signif = signif)
  
  # Save output to dataframe for graphs ---------------------
  
  write_csv(dataset_sales_analysis,file = "4_vet_AMC_analysis_outcomes.csv", append = TRUE)
  
  #close for loop sales ---------------------------
  
}
# #3.2. Sales data graphs --------------------------------------------------

# adding icons to figures ------------------------------------

font_add('fa-solid', '../font_awesome_font_files/fontawesome-free-6.4.0-desktop/otfs/Font Awesome 6 Free-Solid-900.otf')

upward_arrow <- "<span style='font-family:fa-solid'>&#xf062;</span>"
downward_arrow <- "<span style='font-family:fa-solid'>&#xf063;</span>" 
equals <- "<span style='font-family:fa-solid'>&#xf52c;</span>"   #" = " 
oscillate <- "<span style='font-family:fa-solid'>&#xf83e;</span>"   #" ~ "



#load graph data --------------------------------------------------

graph_data_complete <- read.csv("4_vet_AMC_analysis_outcomes.csv", header = TRUE,
                                sep = ",") %>%
  mutate(label_icon = case_when(
    icon == "upward_arrow" ~ upward_arrow,
    icon == "downward_arrow" ~ downward_arrow,
    icon == "equals" ~ equals,
    icon == "oscilate" ~ oscillate
  )) %>%
  mutate(label = if_else(signif == "", label_icon, paste(label_icon,signif,sep=" ")))

Sales_categories

(Sales_graph <- graph_data_complete %>%
    filter(Measurement!="")%>%
    group_by(Measurement) %>%
    mutate(ycoord = min(Value)) %>%
    ungroup()%>%
    mutate(Measurement = fct_relevel(Measurement, c("Total sales (mg/kg biomass)",
                                                    "Sales fluoro(quinolones)+ 3rd- and 4th-generation cephalosporins (mg/kg biomass)",
                                                    "Sales premixes (mg/kg biomass)",                                                  
                                                    "Sales polymyxins (mg/kg biomass)" )))%>%
    ggplot()+
    geom_bar(aes(x = Year, y = Value,fill = Measurement), stat="identity", position = "dodge")+
    geom_hline(aes(yintercept=Target), colour = "red", linetype = "dashed")+
    scale_fill_manual(values = BELMAP_colourscheme)+
    scale_colour_manual(values = BELMAP_colourscheme)+
    facet_wrap( ~ Measurement, nrow = 2, scales = "free_y", labeller = label_wrap_gen(width=55))+
    # geom_richtext( size = 16, hjust = 0, label.colour = NA) +
    scale_x_continuous(limits = c(2010,2024), breaks = seq(2011,2022,1))+
    geom_richtext(aes(x = 2023.5, y = ycoord, 
                      label = label, 
                      fill = Measurement), stat = "unique", 
                  colour= "white", show.legend = FALSE)+
    labs(y= "Sales (mg/kg biomass)")&
    moiras_graph_theme()+
    theme(legend.position = "none",
          strip.text = element_text(size = 14))) 






#4. BD100 Use analysis--------------------------------------------------
start_data_frame <- tibble(
  "Year" = "",
  "Measurement" = "",
  "Value"= "",
  "Unit"= "",
  "Year_simple"= "",
  "icon" = "",
  "signif" = ""
)

# write_csv(start_data_frame,file = "4_vet_AMC_analysis_outcomes_use.csv")
# 
# # total_target= 73.25 #50% 2011 value
# # premix_target= 14.05 #50% 2011 value
# # red_target = 0.5025 # 25% 2011 value - quinolones and 3rd/4th ceph = red category
# # colistin (polymixin = 1)
# 
# 
# # list use indicators
# 
# Use_categories<-Vet_AMC_use$Measurement
# 
# #mkae output test results file
#  capture.output(print("Analysis of Veterinary Use data"),file = "4_Vet_AMC_correlations/use_data_analysis.txt")  
# 
# #make list for overdispersion plots
# qq_plot_list<- list()
# #open for loop use ---------------------------
# #i =  "Total sales \\(mg\\/kg biomass\\)"
# for(i in Use_categories){
# 
#   # 1. check assumptions for Pearsons
#   dataset_use<- Vet_AMC_use %>%
#     filter(grepl(i,Measurement,  fixed = TRUE))
# 
#   capture.output(print(paste("Analysis for",i,sep = " : ")),file = "4_Vet_AMC_correlations/use_data_analysis.txt", append = TRUE)
# 
#   #save qq plot
# 
#   #Normality check
#   # qqnorm(dataset_sales$Year_simple, pch =1, frame = FALSE)
#   # qqline(dataset_sales$Year_simple, col = "steelblue", lwd = 2)
#   #
#   # qqnorm(dataset_sales$Value, pch =1, frame = FALSE)
#   # qqline(dataset_sales$Value, col = "steelblue", lwd = 2)
# 
# 
#   qqPlot(dataset_use$Value)
#   qqPlot_value<- recordPlot()
# 
#   qqPlot(dataset_use$Year)
#   qqPlot_year<- recordPlot()
# 
#   #assign to local environment
#   assign(paste(i,"qqplot_value"),print(qqPlot_value))
#   assign(paste(i,"qqplot_year"),print(qqPlot_year))
# 
#   #Shapiro test: if the  p-values are greater than the significance level 0.05 => the distribution of the data are not significantly different from normal distribution
#   # Shapiro-Wilk normality test for Year_simple
# 
#   Shapiro_year <- shapiro.test(dataset_use$Year_simple) #p-value = 0.8698
#   capture.output(print("Testing Normality of year"),file = "4_Vet_AMC_correlations/use_data_analysis.txt", append = TRUE)
#   capture.output(print(Shapiro_year),file = "4_Vet_AMC_correlations/use_data_analysis.txt", append = TRUE)
# 
# 
#   # Shapiro-Wilk normality test for Value
#   Shapiro_value <- shapiro.test(dataset_use$Value) #p-value = 0.4784
# 
#   capture.output(print("Testing Normality of Values"),file = "4_Vet_AMC_correlations/use_data_analysis.txt", append = TRUE)
#   capture.output(print(Shapiro_value),file = "4_Vet_AMC_correlations/use_data_analysis.txt", append = TRUE)
# 
# 
#   # 2. calculate correlation tests
# 
#   #Pearson (parametric test, assumes linearity and normality) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <cor> and a p-value of <p-value>.
#   Pcor_Value <- cor(dataset_use$Year_simple, dataset_use$Value, method = "pearson")
#   Pcortest_Value <- cor.test(dataset_use$Year_simple, dataset_use$Value, method = "pearson")
#   Pcor_Value
#   #        cor -0.9856587; p-value = 2.905e-08
#   Pcortest_Value
#   #        cor -0.9856587; p-value = 2.905e-08
# 
#   capture.output(print("#Pearson (parametric test, assumes linearity and normality) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <cor> and a p-value of <p-value>.
# "),file = "4_Vet_AMC_correlations/use_data_analysis.txt", append = TRUE)
#   capture.output(print(Pcortest_Value),file = "4_Vet_AMC_correlations/use_data_analysis.txt", append = TRUE)
# 
# 
#   #Kendall (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <tau> and a p-value of <p-value>.
#   Kcor_Value <- cor(dataset_use$Year_simple, dataset_use$Value, method = "kendall")
#   Kcortest_Value <- cor.test(dataset_use$Year_simple, dataset_use$Value, method = "kendall")
#   Kcor_Value
#   # [1] -0.9272727
# 
#   Kcortest_Value
#   #  tau -0.9272727;  p-value = 3.257e-06
# 
#   capture.output(print("#Kendall (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <tau> and a p-value of <p-value>.
# "),file = "4_Vet_AMC_correlations/use_data_analysis.txt",  append = TRUE)
#   capture.output(print(Kcortest_Value),file = "4_Vet_AMC_correlations/use_data_analysis.txt", append = TRUE)
# 
# 
# 
#   #Spearman (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <rho> and a p-value of <p-value>.
#   Scor_Value <- cor(dataset_use$Year_simple, dataset_use$Value, method = "spearman")
#   Scortest_Value <- cor.test(dataset_use$Year_simple, dataset_use$Value, method = "spearman")
#   Scor_Value
#   #[1] -0.9818182
# 
#   Scortest_Value
# 
#   capture.output(print("#Spearman (non-parametric test) if p-value is less than 0.05. We can conclude that Year_simple and Value are significantly correlated with a correlation coefficient of <rho> and a p-value of <p-value>.
# "),file = "4_Vet_AMC_correlations/use_data_analysis.txt",  append = TRUE)
#   capture.output(print(Scortest_Value),file = "4_Vet_AMC_correlations/use_data_analysis.txt", append = TRUE)
# 
# 
#   # manually check that results of correlation tests are aligned
# 
# 
#   # 4.1 select correlation and save icon: use-------------
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
#     if((max(dataset_use$Value)-min(dataset_use$Value))/mean(dataset_use$Value) > 0.25){
#       model_icon <- "oscilate"
#       signif = ""
#     }else if((max(dataset_use$Value)-min(dataset_use$Value))/mean(dataset_use$Value) < 0.25){
#       model_icon <- "equals"
#       signif <- ""
#     }
#   }
# 
# 
#   # Add icons to datasheet-------------------------
#   dataset_use_analysis <- dataset_use %>%
#     mutate(icon = model_icon) %>%
#     mutate(signif = signif)
# 
#   # Save output to dataframe for graphs ---------------------
# 
#   write_csv(dataset_use_analysis,file = "4_vet_AMC_analysis_outcomes_use.csv", append = TRUE)
# 
#   #close for loop use ---------------------------
# 
# }
# #4.2 use data graphs --------------------------------------------------

# adding icons to figures ------------------------------------
font_add("Gill Sans MT", "GIL_____.TTF")

font_add('fa-solid', '../font_awesome_font_files/fontawesome-free-6.4.0-desktop/otfs/Font Awesome 6 Free-Solid-900.otf')

upward_arrow <- "<span style='font-family:fa-solid'>&#xf062;</span>"
downward_arrow <- "<span style='font-family:fa-solid'>&#xf063;</span>" 
equals <- "<span style='font-family:fa-solid'>&#xf52c;</span>"   #" = " 
oscillate <- "<span style='font-family:fa-solid'>&#xf83e;</span>"   #" ~ "



#load graph data --------------------------------------------------

graph_data_complete_use <- read.csv("4_vet_AMC_analysis_outcomes_use.csv", header = TRUE,
                                    sep = ",") %>%
  mutate(label_icon = case_when(
    icon == "upward_arrow" ~ upward_arrow,
    icon == "downward_arrow" ~ downward_arrow,
    icon == "equals" ~ equals,
    icon == "oscilate" ~ oscillate
  )) %>%
  mutate(label = if_else(signif == "", label_icon, paste(label_icon,signif,sep=" ")))

unique(graph_data_complete_use$Measurement)

(use_graph <- graph_data_complete_use %>%
    group_by(Measurement) %>%
    filter(!Measurement=="") %>%
    mutate(ycoord = min(Value)) %>%
    ungroup()%>%
    # mutate(Measurement1 = case_when(
    #   Measurement == "Breeding pigs" ~ "A) Breeding pigs",
    #   Measurement == "Sucklers" ~ "B) Sucklers",
    #   Measurement == "Weaners" ~ "C) Weaners",
    #   Measurement == "Fatteners" ~ "D) Fatteners",
    #   Measurement == "Laying hens" ~ "E) Laying hens",
    #   Measurement == "Broilers" ~ "F) Broilers",
    #   Measurement == "Veal calves" ~ "G) Veal calves"
    # )) %>%
    mutate(Measurement = fct_relevel(Measurement, c("Breeding pigs", "Sucklers",
                                                    "Weaners","Fatteners",
                                                    "Laying hens",  "Broilers",                                               
                                                    "Veal calves" )))%>%
    
    ggplot()+
    geom_bar(aes(x = Year, y = Value,fill = Measurement), stat="identity", position = "dodge")+
    scale_fill_manual(values = BELMAP_colourscheme)+
    scale_colour_manual(values = BELMAP_colourscheme)+
    facet_wrap( ~ Measurement, nrow = 1)+   #, scales = "free_y"
    # geom_richtext( size = 16, hjust = 0, label.colour = NA) +
    scale_x_continuous(limits = c(2017,2024), breaks = seq(2011,2022,1))+
    geom_richtext(aes(x = 2023.5, y = ycoord, 
                      label = label, 
                      fill = Measurement), stat = "unique", 
                  colour= "white", show.legend = FALSE)+
    # fill = NA, label.color = NA, # remove background and outline
    # label.padding = grid::unit(rep(0, 4), "pt")) +
    #  geom_text(aes(x = Year, y = 2, label = Sample_size.x), na.rm = TRUE, size = 2.5)+  #make this additional for interactive report
    labs(y= "Antimicrobial Use - BD100")&
    moiras_graph_theme()+
    theme(legend.position = "none")) 




