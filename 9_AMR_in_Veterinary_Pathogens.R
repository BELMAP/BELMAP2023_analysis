#--------------9. AMR in Veterinary Pathogens-------------------------------------------------------

# This script outlines the data analysis for Chapter 9 of the 2023 BELMAP report, namely it:
# a)	Uploads and combines the data
# b)	Cleans names (e.g. some in all caps etc.)
# c)	Per pathogen /indicator:
#   i.	Fits glm models with poisson and neg binom distributions
#   ii.	Checks for overdispersion, model fit and zero inflation, 
#   iii.	Compares these test values – if there is no overdispersion and the poisson model fit is better (by AIC and chi-square anova) select poisson model, if there is overdispersion and NB model fit is better select NB model, if there is disparity then flag a warning but select the model with best fit. 
#   iv.	Output message from zero inflation test and print all model and test outputs to a log file per host/pathogen/indicator for manual inspection, particularly if warning is flagged. – I no longer automatically run the zero-inflated models as the majority of the data contains no zeros- but the script includes the script to run and assess the ZI models, commented out, that can be run for specific cases if a warning of possible ZI is identified.
#   v.	Generate a dataset of predicted values based on the best fit model, which is used to plot the best fit model
#   vi.	Add Confidence intervals to this predicted dataset – based on the link function of the model (extract link function fit and standard error from best fit model, make columns of CI (fitted value +/- 2*SE based on link function, https://fromthebottomoftheheap.net/2018/12/10/confidence-intervals-for-glms/),
#   vii.	Identify p-value and coefficient for time variable in best fit model: 
#             •	if model shows significant change over time then identify direction of coefficient – if positive then assign model label to be upward arrow plus number of * based on model p-value,  
#             •	if model shows significant change over time and direction of coefficient is negative then assign model label to be downward arrow plus number of * based on model p-value; 
#             •	if non-significant but the range (max-min value) is greater than 25% of mean then assign “oscillating” icon, 
#             •	if non-significant but range < 0.25% of mean then assign “equals” icon. 
#   viii.	Generate plot coordinates for model icons based on ranking and number of indicators per graph (to standardise location of model direction icons in figures)
#   ix.	Combine real dataset with predicted dataset and model labels
#Contents-----
#1. Load libraries and themes
#2. load and combine datasets
#3. list pathogens and indicators - make 
#4. Run Models per pathogen, save to database, not mastitis 
#5. Run Models per pathogen for mastitis, save to database

#-----------------------------------------------------------------------------------
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
library(ggpattern)
library(ggtext)
library(ggrepel)
library(showtext)
#load.fontawesome()
# install.packages("showtext")

BELMAP_colourscheme <- c("#1f78b4","#a6cee3","#b2df8a","#fb9a99","#33a02c","#e31a1c",
                         "#fdbf6f","#ff7f00","#cab2d6", "#6a3d9a","#ffff99") #"MDR"="#e31a1c","pan-S" = "#33a02c",



BELMAP_Human_AntiB_colourscheme <- c("Ampicillin" ="#1f78b4",
                                     "Methicillin" ="#1f78b4",
                                     "Penicillin"="#1f78b4",
                                     "Carbapenem" = "#a6cee3",
                                     "Azole" = "#1f78b4",
                                     "Fluconazole" = "#1f78b4",
                                     "Colistin" = "#b2df8a",
                                     "3GC" = "#fb9a99",
                                     "Ceftriaxone"	= "#fb9a99",
                                     "ESBL"= "#fb9a99",
                                     "vancomycin" = "#33a02c",
                                     "Vancomycin" = "#33a02c",
                                     "Anidulafungine" = "#a6cee3",
                                     "Monoresistant isoniazid" = "#1f78b4",
                                     "Rifampicin" = "#a6cee3",
                                     "Fidaxomicin" ="#e31a1c",
                                     "Ciprofloxacin" = "#fdbf6f",
                                     "Moxifloxacin" = "#fdbf6f",
                                     "Macrolide" = "#ff7f00",
                                     "Clarithomycin" = "#ff7f00",
                                     "Azithromycin" = "#ff7f00",
                                     "Macrolide and Fluoroquinolone" = "#1a1a1a",
                                     "AG, FQ and 3GC" = "#1a1a1a",
                                     "MDR" = "#1a1a1a")





moiras_graph_theme<- function(..., base_size = 12){theme(
  panel.background = element_blank(), #transparent panel bg
  plot.background = element_blank(), #transparent plot bg
  panel.grid.major.y =  element_line(colour = "#ccc6c4"), #grey y major gridlines
  panel.grid.major.x =  element_blank(), #remove major gridlines
  panel.grid.minor = element_blank(), #remove minor gridlines
  # legend.background = element_rect(fill='#001f3f'), #transparent legend bg
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 14),
  #  legend.box.background = element_rect(fill='#001f3f'),
  text = element_text( family = "Gill Sans MT"),
  axis.text = element_text(size = 14),
  #plot.background=element_blank(),#, , size = 10, family = "calibri"
  axis.title = element_text(size = 16),
  axis.text.x = element_text(angle = -90), #vjust = 1, hjust = -1
  strip.background = element_blank(),
  strip.text = element_text(size = 16),
  legend.position = "top",
  legend.direction = "horizontal"
)}

font_add("Gill Sans MT", "GIL_____.TTF")


# Important step to enable showtext font rendering!
showtext_auto()


#2. Load and combine datasets---------------------------------------------------


##load data
# ARISA DATA
Vet_path_ARSIA_raw <-  read_csv2("Data/2.AMR/Veterinary_AMR/Veterinary_AMR_2023_ARSIA.csv",
                       #col_types = cols(
                       #  "NumValue" = col_double()#,
                       # "Sample_size"= col_double(),
                       # "Percent_resistant"= col_double(),
                       # "Number_resistant"= col_double(),
                       # "Pathogen" = col_character(),
                       # "Host" = col_character(),
                       # "Indicator" = col_character()
                       #),
                       na = c("", "NA", "-"),
                       trim_ws = TRUE
) %>%
  mutate(Indication = Source) %>%
  dplyr::select(Year,Antibiotic,Percent_resistance,Sample_size,Number_resistant,Source,Bacteria, Indication)%>%
  mutate(Organisation = "ARSIA")

# DGZ DATA
Vet_path_DGZ_raw <-  read_csv2("Data/2.AMR/Veterinary_AMR/Veterinary_AMR_2023_DGZ_MCC.csv",
                                 #col_types = cols(
                                 #  "NumValue" = col_double()#,
                                 # "Sample_size"= col_double(),
                                 # "Percent_resistant"= col_double(),
                                 # "Number_resistant"= col_double(),
                                 # "Pathogen" = col_character(),
                                 # "Host" = col_character(),
                                 # "Indicator" = col_character()
                                 #),
                                 na = c("", "NA", "-"),
                                 trim_ws = TRUE
) %>%
  mutate(Indication = if_else(is.na(Indication), Source, Indication)) %>%
  dplyr::select(Year,Antibiotic,Percent_resistance,Sample_size,Number_resistant,Source,Bacteria,Indication) %>%
  mutate(Organisation = "DGZ_MCC")

# DATA ON MASTITIS SAMPLES
mastitis_data<- read_csv2("Data/2.AMR/Veterinary_AMR/Veterinary_AMR_2023.csv",
                          #col_types = cols(
                          #  "NumValue" = col_double()#,
                          # "Sample_size"= col_double(),
                          # "Percent_resistant"= col_double(),
                          # "Number_resistant"= col_double(),
                          # "Pathogen" = col_character(),
                          # "Host" = col_character(),
                          # "Indicator" = col_character()
                          #),
                          na = c("", "NA", "-"),
                          trim_ws = TRUE
) %>%
  mutate(Indication = Source) %>%
  filter(grepl("Cattle_udder",Source)) %>%
  dplyr::select(Year,Antibiotic,Percent_resistance,Sample_size,Number_resistant,Source,Bacteria,Indication) %>%
  mutate(Organisation = "MCC_2022")

# combine datasets

vet_path_data_raw <- rbind(Vet_path_DGZ_raw,Vet_path_ARSIA_raw)

vet_path_data_clean_analysis <- vet_path_data_raw %>%
  group_by(Year,Antibiotic,Source,Bacteria) %>%
  summarise(Sample_size = sum(Sample_size, na.rm=TRUE),
            Number_resistant = sum(Number_resistant, na.rm=TRUE),
            Percent_resistant = Number_resistant/Sample_size*100) %>%
  ungroup() %>%
  filter(!(Sample_size == 0)) %>%
  mutate(category = paste(Bacteria,Antibiotic,Source, sep = "_"))

mastitis_data_combined <- rbind(mastitis_data,vet_path_data_raw) %>%
  filter(grepl("Cattle_udder",Source)) %>%
  filter(!is.na(Sample_size)) %>%
  arrange(Bacteria, Antibiotic, Year) %>%
  mutate(Antibiotic1= str_to_title(Antibiotic)) %>%
  mutate(Antibiotic = case_when(
    grepl("Tmpsmx",Antibiotic1) ~ "TMP-SMX",
    grepl("Penicilline_g",Antibiotic1) ~ "Penicilline G",
    grepl("Tmp-Smx",Antibiotic1) ~ "TMP-SMX",
    grepl("Amx_clav",Antibiotic1) ~ "Amoxyclav",
    .default = Antibiotic1)) %>%
  mutate(category = paste(Indication,Organisation,sep ="_")) %>%
 # mutate(category_nr = paste(Indication,Organisation,"num_resist", sep ="_")) %>%
    dplyr::select(-c(Organisation,Indication, Percent_resistance)) %>%
  pivot_wider(names_from = category, values_from = c(Sample_size,Number_resistant)) %>%
mutate(Sample_size_clinical = case_when(
  Sample_size_Cattle_udder_MCC_2022 == `Sample_size_Clinical mastitis_DGZ_MCC` ~ Sample_size_Cattle_udder_MCC_2022,
  is.na(Sample_size_Cattle_udder_MCC_2022) ~ `Sample_size_Clinical mastitis_DGZ_MCC`,
  is.na(`Sample_size_Clinical mastitis_DGZ_MCC`) ~ Sample_size_Cattle_udder_MCC_2022
)) %>%
  mutate(Number_resistant_clinical = case_when(
    Number_resistant_Cattle_udder_MCC_2022 == `Number_resistant_Clinical mastitis_DGZ_MCC` ~ Number_resistant_Cattle_udder_MCC_2022,
    is.na(Number_resistant_Cattle_udder_MCC_2022) ~ `Number_resistant_Clinical mastitis_DGZ_MCC`,
    is.na(`Number_resistant_Clinical mastitis_DGZ_MCC`) ~ Number_resistant_Cattle_udder_MCC_2022
  )) %>%
  dplyr::select(Year, Antibiotic,Source,Bacteria, Number_resistant_clinical, Sample_size_clinical, 
                `Sample_size_Subclinical mastitis_DGZ_MCC`,`Number_resistant_Subclinical mastitis_DGZ_MCC`) %>%
  # pivot_longer(cols = c(Number_resistant_clinical, `Number_resistant_Subclinical mastitis_DGZ_MCC`),
  #               names_to = "Measurement", values_to = "Number_resistant") %>%
  # mutate(Measurement = str_remove(Measurement, "Number_resistant_")) %>%
  # pivot_longer(cols = c(Sample_size_clinical, `Sample_size_Subclinical mastitis_DGZ_MCC`),names_pattern = "Sample_size_(.*)",
  #              names_to = "Measurement1", values_to = "Sample_size")
  mutate(category = paste(Antibiotic,Bacteria, sep="_"))

#4. Run Models not mastitis -------------------------------------------------

#make list of categories

combined_categories<-unique(vet_path_data_clean_analysis$category)

#write_csv2(vet_path_data_clean_analysis, "vet_path_data_clean_analysis.csv")

#calc degree variance for categories----------

var_list = list()


for(i in combined_categories){
  print(i)
  dataset_analysis<- vet_path_data_clean_analysis %>%
    filter(category == i) 
  degree_var= (max(dataset_analysis$Percent_resistant, na.rm = TRUE)-min(dataset_analysis$Percent_resistant, na.rm = TRUE))/mean(as.numeric(dataset_analysis$Percent_resistant), na.rm = TRUE)
  var_list[[i]] <-degree_var
  print(degree_var)
}

#make list for overdispersion plots
overdispersion_plot_list<- list()


#make csv for graph data with header

start_data_frame <- tibble(
  "Year" = "",
  "Bacteria" = "",
  "Antibiotic"= "",
  "Percent_resistant"= "",
  "Percent_resistant_predict"= "",
  "Source" = "",
  "CI_upper"= "",
  "CI_lower"= "",
  "Sample_size.x" = "",
  "icon" = "",
  "signif_level" = ""
)



#write_csv(start_data_frame,file = "9_Vet_AMR_data_and_GLM_predictions.csv")

# make for loop to run analysis-----
for(i in combined_categories){
  
  #4.0.subset data-----
  
  dataset_analysis_raw<- vet_path_data_clean_analysis %>%
    filter(grepl(i,category)) %>%
    mutate(Year = as.numeric(Year))
  min_year = min(dataset_analysis_raw$Year)
  max_year = max(dataset_analysis_raw$Year)
  dataset_analysis <- dataset_analysis_raw %>%  
    mutate(Year_simple = Year - min_year) %>%
    distinct()
  
  #make a file name to print each analysis
  
  output_file_name <- paste("9_Veterinary_AMR_GLM_output_summaries/GLM_output_",i,".txt",sep="") 
  
  #print title of analysis
  capture.output(print(paste("GLM analysis of ",i,sep="")),file = output_file_name)
  
  if(max_year-min_year < 3){
    capture.output(print("Too few data points for model"),file = output_file_name)
    graph_data <- dataset_analysis %>%
      mutate(Percent_resistant_predict = "",
             CI_upper = "",
             CI_lower = "",
             Sample_size.x = Sample_size,
             icon = "",
             signif_level = "") %>%
      dplyr::select(Year,Bacteria,Antibiotic, Percent_resistant, Percent_resistant_predict,Source,
                    CI_upper,CI_lower,Sample_size.x,icon, signif_level)
  } else if(sum(as.numeric(dataset_analysis$Percent_resistant)) == 0){
    capture.output(print("No resistance- can't fit model"),file = output_file_name)
    graph_data <- dataset_analysis %>%
      mutate(Percent_resistant_predict = "",
             CI_upper = "",
             CI_lower = "",
             Sample_size.x = Sample_size,
             icon = "",
             signif_level = "") %>%
      dplyr::select(Year,Bacteria,Antibiotic, Percent_resistant, Percent_resistant_predict,Source,
                    CI_upper,CI_lower,Sample_size.x,icon, signif_level)
  } else {
    
    #4.1. fit models:-----
    #4.1a. fit poisson-----
    
    glmpoissonirr <- glm(Number_resistant~Year_simple + offset(log(Sample_size)), data = dataset_analysis, family = poisson(link = "log"))
    
    
    #4.1b. fit NB-----
    
    nb <- glm.nb(Number_resistant~Year_simple + offset(log(Sample_size )), data = dataset_analysis )
    
    
    #--------------------------------------
    
    # #4.2.check for overdispersion-----
    # #visualise overdispersion 
    # plot(log(fitted(glmpoissonirr)),log((dataset_analysis$Number_resistant-fitted(glmpoissonirr))^2),xlab=expression(hat(mu)),ylab=expression((y-hat(mu))^2),pch=20,col="blue")
    # abline(0,1) ## 'variance = mean' line
    # 
    # overdisperions_plot<- recordPlot()
    # 
    # #assign to local environment
    # assign(paste(i,"overdispersion_plot"),print(overdisperions_plot))
    # 
    # #print dispersal graph to pdf
    # dispersion_graph_file<-paste("9_Vet_AMR_Dispersion_graphs/Dispersion_graph_poisson_model",i,".pdf",sep="") 
    # pdf(dispersion_graph_file, width = 8, height = 6)
    # # Call plot
    # overdisperions_plot
    # # Closing the graphical device
    # dev.off() 
    # 
    # calculate dispersion parameter
    dp = sum(residuals(glmpoissonirr,type ="pearson")^2)/glmpoissonirr$df.residual
    
    #print over vs underdispersion
    
    
    underdispersed<- dispersiontest(glmpoissonirr,alternative = c("less") )
    
    
    #4.3 select model and make predict dataset----------------
    #MAKE COMPARISONS
    
    #dispersion
    
    overdispersed<- dispersiontest(glmpoissonirr,alternative = c("greater") )
    
    #anova
    pois_vs_NB<- lrtest(glmpoissonirr,nb)
    
    #aic
    # Poisson AIC    
    ##extract log-likelihood
    LL <- logLik(glmpoissonirr)[1]
    ##extract number of parameters
    K.mod <- coef(glmpoissonirr) + 1
    ##compute AICc with full likelihood
    AICcCustom(LL, K.mod, nobs = nrow(dataset_analysis))
    
    AIC_poisson<- AICcCustom(LL, K.mod, nobs = nrow(dataset_analysis))   
    AIC_poisson_value<- AICcCustom(LL, K.mod, nobs = nrow(dataset_analysis))[1]
    
    # NB model AIC 
    AIC_NB <- AIC(nb)
    
    
    #select model
    if(
      (pois_vs_NB$`Pr(>Chisq)`[2] < 0.05) ){   #& (AIC_poisson_value > AIC_NB) 
      selected_model <- nb
      if(overdispersed$p.value > 0.05){
        print(paste("CHECK OUTPUTS - NB selected but no overdispersal for ", i ))
      }
    }else if(
      (pois_vs_NB$`Pr(>Chisq)`[2] > 0.05) ){  #& (AIC_poisson_value < AIC_NB)
      selected_model <- glmpoissonirr
      if(overdispersed$p.value < 0.05){
        print(paste("CHECK OUTPUTS - poisson selected but poss overdispersal for ",i))
      }
    }else if(grepl("pigs|SOW",i)){  
      selected_model <- NA  # no model selection for MRSA in pigs as changes in methodology
    }else{
      print(paste("ERROR no model selected for",i))
    }
    
    #make predicted dataset - save it 
    length_years = max_year - min_year
    
    ndata<- tibble(
      Year = seq(min_year,max_year,by =1),
      Year_simple = seq(0,length_years,by = 1),
      Sample_size = rep(200,length_years+1),
      category = rep(i,length_years+1)
    )
    
    ndata <- add_column(ndata, fit = predict(selected_model, newdata = ndata, type = 'response'))
    
    #make confidence intervals - based on link function
    #find inverse of link function:
    ilink<-family(selected_model)$linkinv
    
    ndata <- bind_cols(ndata, setNames(as_tibble(predict(selected_model, ndata, se.fit = TRUE)[1:2]),
                                       c('fit_link','se_link')))
    
    ndata <- mutate(ndata,
                    fit_resp  = ilink(fit_link),
                    right_upr = ilink(fit_link + (2 * se_link)),
                    right_lwr = ilink(fit_link - (2 * se_link)))
    
    
    # extract pvalue of model
    selected_model_p_value <- coef(summary(selected_model))[,'Pr(>|z|)'][2]
    
    #if signif - extract direction of model
    if(selected_model_p_value<0.05){
      if(coef(summary(selected_model))[,'Estimate'][2] > 0){
        model_icon <- "upward_arrow"
        signif = case_when(
          selected_model_p_value > 0.01 && selected_model_p_value < 0.05 ~ "*",
          selected_model_p_value > 0.001 && selected_model_p_value < 0.01 ~ "**",
          selected_model_p_value < 0.001 ~ "***"
        )
      }else if(coef(summary(selected_model))[,'Estimate'][2] < 0){
        model_icon = "downward_arrow"
        signif = case_when(
          selected_model_p_value > 0.01 && selected_model_p_value < 0.05 ~ "*",
          selected_model_p_value > 0.001 && selected_model_p_value < 0.01 ~ "**",
          selected_model_p_value < 0.001 ~ "***"
        )
      }
    } else if(selected_model_p_value>0.05){
      if((max(dataset_analysis$Percent_resistant)-min(dataset_analysis$Percent_resistant))/mean(dataset_analysis$Percent_resistant) > 0.25){
        model_icon <- "oscilate"
        signif = ""
      }else if((max(dataset_analysis$Percent_resistant)-min(dataset_analysis$Percent_resistant))/mean(dataset_analysis$Percent_resistant) < 0.25){
        model_icon <- "equals"
        signif <- ""
      }
    }
    
    predict_data <- ndata %>%
      mutate(Percent_resistant_predict = fit/Sample_size*100,
             CI_upper = right_upr/Sample_size*100,
             CI_lower = right_lwr/Sample_size*100,
             icon = model_icon,
             signif_level = signif)
    
    
    # if non-signif - decide if equals or oscilating
    #calc if range > mean/4
    
    
    
    graph_data <- left_join(dataset_analysis,predict_data, by = c("Year","Year_simple")) %>%
      dplyr::select(Year,Bacteria,Antibiotic, Percent_resistant, Percent_resistant_predict,Source,
                      CI_upper,CI_lower,Sample_size.x,icon, signif_level)
    
    
    
    
    #save dataset
    
    predicted_dataset_filename <- paste("9_Vet_AMR_GLM_predictions/GLM_predictions_",i,".csv",sep="") 
    
    write_csv(predict_data,file = predicted_dataset_filename)  
    write_csv(graph_data,file = "9_Vet_AMR_data_and_GLM_predictions.csv", append = TRUE)
    
    #print selected model to GLM output
    capture.output(print("SELECTED MODEL----------------------------"),file = output_file_name, append = "TRUE")
    capture.output(print(selected_model),file = output_file_name, append = "TRUE")
    capture.output(print(summary(selected_model)),file = output_file_name, append = "TRUE")
    
    
    #4.4 save outputs to manually check ------------
    
    
    #print poisson p value to 
    
    capture.output(print("TEST FOR DISPERSION-------------------- "),file = output_file_name, append = "TRUE") 
    capture.output(print(paste("dispersal parameter= ",dp,"If >1 then overdispersed")),file = output_file_name, append = "TRUE")
    
    capture.output(print(overdispersed),file = output_file_name, append = "TRUE")
    
    capture.output(print(underdispersed),file = output_file_name, append = "TRUE")
    #  --------------------------------------------------------------------------------
    
    capture.output(print("ANOVA COMPARING MODELS-------------------- "),file = output_file_name, append = "TRUE")
    
    #3.6. compare model likelihood ratio - NB vs poisson, ZI NB vs NB 
    #3.6.a poisson vs NB
    
    pois_vs_NB<- lrtest(glmpoissonirr,nb)
    
    capture.output(print("Likelihood Ratio Test Poisson vs NB model   "),file = output_file_name, append = "TRUE")
    capture.output(print(pois_vs_NB),file = output_file_name, append = "TRUE")
    
    
    # Poisson AIC
    capture.output(print(paste("AIC of Poisson model:",AIC_poisson)),file = output_file_name, append = "TRUE")
    
    capture.output(print(paste("AIC of NB model:", AIC_NB)),file = output_file_name, append = "TRUE")
    
    #  --------------------------------------------------------------------------------
    capture.output(print("MODELS-------------------- "),file = output_file_name, append = "TRUE")
    
    capture.output(print("Output from Poisson model"),file = output_file_name, append = "TRUE")
    capture.output(print(summary(glmpoissonirr)),file = output_file_name, append = "TRUE")
    #  capture.output(print(anova(glmpoissonirr,test='Chi')),file = output_file_name, append = "TRUE")
    
    capture.output(print("Output from Neg Binomial model"),file = output_file_name, append = "TRUE")
    capture.output(print(summary(nb)),file = output_file_name, append = "TRUE")
    # capture.output(print(anova(nb,test='Chi')),file = output_file_name, append = "TRUE")
    
    
    
    #check poisson zero inflation - if found manually run zero inflated models
    print(i)
    pois_ZI <- check_zeroinflation(glmpoissonirr, tolerance = 0.05)
    
    capture.output(print("Poisson ZI test= "),file = output_file_name, append = "TRUE")
    
    capture.output(print(pois_ZI),file = output_file_name, append = "TRUE")
    #check NB zero inflation
    
    NB_ZI <- check_zeroinflation(nb, tolerance = 0.05)
    
    capture.output(print("NB ZI test= "),file = output_file_name, append = "TRUE")
    capture.output(print(NB_ZI),file = output_file_name, append = "TRUE")
  }
  
  
  #close for loop -----
}



#5. Run Models mastitis -------------------------------------------------

#make list of categories

combined_categories_mast<-unique(mastitis_data_combined$category)

#checking model for Penicilline G_S. aureus
combined_categories_mast<-"Penicilline G_S. aureus"

#calc degree variance for categories----------

var_list = list()

for(i in combined_categories_mast){
  print(i)
  dataset_analysis<- mastitis_data_combined %>%
    filter(category == i) %>%
    mutate(Percent_resistant_clin = Number_resistant_clinical/Sample_size_clinical*100)
  degree_var= (max(dataset_analysis$Percent_resistant_clin, na.rm = TRUE)-min(dataset_analysis$Percent_resistant_clin, na.rm = TRUE))/mean(as.numeric(dataset_analysis$Percent_resistant_clin), na.rm = TRUE)
  var_list[[i]] <-degree_var
  print(degree_var)
}

#make list for overdispersion plots
overdispersion_plot_list<- list()


#make csv for graph data with header

start_data_frame <- tibble(
  "Year" = "",
  "Bacteria" = "",
  "Antibiotic"= "",
  "Percent_resistant_clinical"= "",
  "Percent_resistant_clinical_predict"= "",
  "Percent_resistant_subclinical"= "",
  "CI_upper"= "",
  "CI_lower"= "",
  "Sample_size_clinical" = "",
  "Sample_size_Subclinical mastitis_DGZ_MCC"  = "",
  "icon" = "",
  "signif_level" = ""
)


#write_csv(start_data_frame,file = "9_Vet_AMR_data_and_GLM_predictions_mast.csv")



# make for loop to run analysis-----
for(i in combined_categories_mast){
  
  #5.0.subset data-----
  
  dataset_analysis_raw<- mastitis_data_combined %>%
    filter(category == i) %>%
    mutate(Year = as.numeric(Year)) %>%
    mutate(Percent_resistant_clinical = Number_resistant_clinical/Sample_size_clinical*100) %>%
    mutate(Percent_resistant_subclinical = `Number_resistant_Subclinical mastitis_DGZ_MCC`/`Sample_size_Subclinical mastitis_DGZ_MCC`*100)
  min_year = min(dataset_analysis_raw$Year)
  max_year = max(dataset_analysis_raw$Year)
  dataset_analysis <- dataset_analysis_raw %>%  
    mutate(Year_simple = Year - min_year) %>%
    distinct()
  
  #make a file name to print each analysis
  
  output_file_name <- paste("9_Veterinary_AMR_GLM_output_summaries/GLM_output_",i,"_mast.txt",sep="") 
  
  #print title of analysis
  capture.output(print(paste("GLM analysis of ",i,sep="")),file = output_file_name)
  
  if(max_year-min_year < 3){
    capture.output(print("Too few data points for model"),file = output_file_name)
    graph_data <- dataset_analysis %>%
      mutate(Percent_resistant_clinical_predict = "",
             CI_upper = "",
             CI_lower = "",
             icon = "",
             signif_level = "") %>%
      dplyr::select(Year,Bacteria,Antibiotic, Percent_resistant_clinical,
                    Percent_resistant_clinical_predict, Percent_resistant_subclinical,
                    CI_upper,CI_lower,Sample_size_clinical,
                    `Sample_size_Subclinical mastitis_DGZ_MCC`,icon, signif_level)
  } else if(sum(as.numeric(dataset_analysis$Percent_resistant_clinical), na.rm = TRUE) == 0){
    capture.output(print("No resistance- can't fit model"),file = output_file_name)
    graph_data <- dataset_analysis %>%
      mutate(Percent_resistant_clinical_predict = "",
             CI_upper = "",
             CI_lower = "",
             icon = "",
             signif_level = "") %>%
      dplyr::select(Year,Bacteria,Antibiotic, Percent_resistant_clinical,
                    Percent_resistant_clinical_predict, Percent_resistant_subclinical,
                    CI_upper,CI_lower,Sample_size_clinical,
                    `Sample_size_Subclinical mastitis_DGZ_MCC`,icon, signif_level)
    } else {

  
    #5.1. fit models:-----
    #5.1a. fit poisson-----
    
    glmpoissonirr <- glm(Number_resistant_clinical~Year_simple + offset(log(Sample_size_clinical)), data = dataset_analysis, family = poisson(link = "log"))
    
    
    #5.1b. fit NB-----
    
    nb <- glm.nb(Number_resistant_clinical~Year_simple + offset(log(Sample_size_clinical)), data = dataset_analysis,control=glm.control(maxit=50) )
    
    
    #--------------------------------------
    
    # #5.2.check for overdispersion-----
    # #visualise overdispersion 
    # plot(log(fitted(glmpoissonirr)),log((dataset_analysis$Number_resistant-fitted(glmpoissonirr))^2),xlab=expression(hat(mu)),ylab=expression((y-hat(mu))^2),pch=20,col="blue")
    # abline(0,1) ## 'variance = mean' line
    # 
    # overdisperions_plot<- recordPlot()
    # 
    # #assign to local environment
    # assign(paste(i,"overdispersion_plot"),print(overdisperions_plot))
    # 
    # #print dispersal graph to pdf
    # dispersion_graph_file<-paste("9_Vet_AMR_Dispersion_graphs/Dispersion_graph_poisson_model",i,".pdf",sep="") 
    # pdf(dispersion_graph_file, width = 8, height = 6)
    # # Call plot
    # overdisperions_plot
    # # Closing the graphical device
    # dev.off() 
    # 
    # calculate dispersion parameter
    dp = sum(residuals(glmpoissonirr,type ="pearson")^2)/glmpoissonirr$df.residual
    
    #print over vs underdispersion
    
    
    underdispersed<- dispersiontest(glmpoissonirr,alternative = c("less") )
    
    
    #5.3 select model and make predict dataset----------------
    #MAKE COMPARISONS
    
    #dispersion
    
    overdispersed<- dispersiontest(glmpoissonirr,alternative = c("greater") )
    
    #anova
    pois_vs_NB<- lrtest(glmpoissonirr,nb)
    
    #aic
    # Poisson AIC    
    ##extract log-likelihood
    LL <- logLik(glmpoissonirr)[1]
    ##extract number of parameters
    K.mod <- coef(glmpoissonirr) + 1
    ##compute AICc with full likelihood
    AICcCustom(LL, K.mod, nobs = nrow(dataset_analysis))
    
    AIC_poisson<- AICcCustom(LL, K.mod, nobs = nrow(dataset_analysis))   
    AIC_poisson_value<- AICcCustom(LL, K.mod, nobs = nrow(dataset_analysis))[1]
    
    # NB model AIC 
    AIC_NB <- AIC(nb)
    
    
    #select model
    if(
      (pois_vs_NB$`Pr(>Chisq)`[2] < 0.05) ){   #& (AIC_poisson_value > AIC_NB) 
      selected_model <- nb
      if(overdispersed$p.value > 0.05){
        print(paste("CHECK OUTPUTS - NB selected but no overdispersal for ", i ))
      }
    }else if(
      (pois_vs_NB$`Pr(>Chisq)`[2] > 0.05) ){  #& (AIC_poisson_value < AIC_NB)
      selected_model <- glmpoissonirr
      if(overdispersed$p.value < 0.05){
        print(paste("CHECK OUTPUTS - poisson selected but poss overdispersal for ",i))
      }
    }else{
      print(paste("ERROR no model selected for",i))
    }
    
    #make predicted dataset - save it 
    length_years = max_year - min_year
    
    ndata<- tibble(
      Year = seq(min_year,max_year,by =1),
      Year_simple = seq(0,length_years,by = 1),
      Sample_size_clinical = rep(200,length_years+1),
      category = rep(i,length_years+1)
    )
    
    ndata <- add_column(ndata, fit = predict(selected_model, newdata = ndata, type = 'response'))
    
    #make confidence intervals - based on link function
    #find inverse of link function:
    ilink<-family(selected_model)$linkinv
    
    ndata <- bind_cols(ndata, setNames(as_tibble(predict(selected_model, ndata, se.fit = TRUE)[1:2]),
                                       c('fit_link','se_link')))
    
    ndata <- mutate(ndata,
                    fit_resp  = ilink(fit_link),
                    right_upr = ilink(fit_link + (2 * se_link)),
                    right_lwr = ilink(fit_link - (2 * se_link)))
    
    
    # extract pvalue of model
    selected_model_p_value <- coef(summary(selected_model))[,'Pr(>|z|)'][2]
    
    #if signif - extract direction of model
    if(selected_model_p_value<0.05){
      if(coef(summary(selected_model))[,'Estimate'][2] > 0){
        model_icon <- "upward_arrow"
        signif = case_when(
          selected_model_p_value > 0.01 && selected_model_p_value < 0.05 ~ "*",
          selected_model_p_value > 0.001 && selected_model_p_value < 0.01 ~ "**",
          selected_model_p_value < 0.001 ~ "***"
        )
      }else if(coef(summary(selected_model))[,'Estimate'][2] < 0){
        model_icon = "downward_arrow"
        signif = case_when(
          selected_model_p_value > 0.01 && selected_model_p_value < 0.05 ~ "*",
          selected_model_p_value > 0.001 && selected_model_p_value < 0.01 ~ "**",
          selected_model_p_value < 0.001 ~ "***"
        )
      }
    } else if(selected_model_p_value>0.05){
      if((max(dataset_analysis$Percent_resistant_clinical, na.rm = TRUE)-min(dataset_analysis$Percent_resistant_clinical, na.rm = TRUE))/mean(dataset_analysis$Percent_resistant_clinical, na.rm = TRUE) > 0.25){
        model_icon <- "oscilate"
        signif = ""
      }else if((max(dataset_analysis$Percent_resistant_clinical, na.rm = TRUE)-min(dataset_analysis$Percent_resistant_clinical, na.rm = TRUE))/mean(dataset_analysis$Percent_resistant_clinical, na.rm = TRUE) < 0.25){
        model_icon <- "equals"
        signif <- ""
      }
    }
    
    predict_data <- ndata %>%
      mutate(Percent_resistant_clinical_predict = fit/Sample_size_clinical*100,
             CI_upper = right_upr/Sample_size_clinical*100,
             CI_lower = right_lwr/Sample_size_clinical*100,
             icon = model_icon,
             signif_level = signif)
    
    
    # if non-signif - decide if equals or oscilating
    #calc if range > mean/4
    
    
    
    graph_data <- left_join(dataset_analysis,predict_data, by = c("Year","Year_simple")) %>%
      dplyr::select(Year,Bacteria,Antibiotic, Percent_resistant_clinical, 
                    Percent_resistant_clinical_predict, Percent_resistant_subclinical,
                    CI_upper,CI_lower,Sample_size_clinical.x,
                    `Sample_size_Subclinical mastitis_DGZ_MCC`,icon, signif_level)
    
    
    
    
    #save dataset
    
    predicted_dataset_filename <- paste("9_Vet_AMR_GLM_predictions/GLM_predictions_",i,"_mast.csv",sep="") 
    
    write_csv(predict_data,file = predicted_dataset_filename)  
    write_csv(graph_data,file = "9_Vet_AMR_data_and_GLM_predictions_mast.csv", append = TRUE)
    
    #print selected model to GLM output
    capture.output(print("SELECTED MODEL----------------------------"),file = output_file_name, append = "TRUE")
    capture.output(print(selected_model),file = output_file_name, append = "TRUE")
    capture.output(print(summary(selected_model)),file = output_file_name, append = "TRUE")
    
    
    #3.4 save outputs to manually check ------------
    
    
    #print poisson p value to 
    
    capture.output(print("TEST FOR DISPERSION-------------------- "),file = output_file_name, append = "TRUE") 
    capture.output(print(paste("dispersal parameter= ",dp,"If >1 then overdispersed")),file = output_file_name, append = "TRUE")
    
    capture.output(print(overdispersed),file = output_file_name, append = "TRUE")
    
    capture.output(print(underdispersed),file = output_file_name, append = "TRUE")
    #  --------------------------------------------------------------------------------
    
    capture.output(print("ANOVA COMPARING MODELS-------------------- "),file = output_file_name, append = "TRUE")
    
    #3.6. compare model likelihood ratio - NB vs poisson, ZI NB vs NB 
    #3.6.a poisson vs NB
    
    pois_vs_NB<- lrtest(glmpoissonirr,nb)
    
    capture.output(print("Likelihood Ratio Test Poisson vs NB model   "),file = output_file_name, append = "TRUE")
    capture.output(print(pois_vs_NB),file = output_file_name, append = "TRUE")
    
    
    # Poisson AIC
    capture.output(print(paste("AIC of Poisson model:",AIC_poisson)),file = output_file_name, append = "TRUE")
    
    capture.output(print(paste("AIC of NB model:", AIC_NB)),file = output_file_name, append = "TRUE")
    
    #  --------------------------------------------------------------------------------
    capture.output(print("MODELS-------------------- "),file = output_file_name, append = "TRUE")
    
    capture.output(print("Output from Poisson model"),file = output_file_name, append = "TRUE")
    capture.output(print(summary(glmpoissonirr)),file = output_file_name, append = "TRUE")
    #  capture.output(print(anova(glmpoissonirr,test='Chi')),file = output_file_name, append = "TRUE")
    
    capture.output(print("Output from Neg Binomial model"),file = output_file_name, append = "TRUE")
    capture.output(print(summary(nb)),file = output_file_name, append = "TRUE")
    # capture.output(print(anova(nb,test='Chi')),file = output_file_name, append = "TRUE")
    
    
    
    #check poisson zero inflation - if found manually run zero inflated models
    print(i)
    pois_ZI <- check_zeroinflation(glmpoissonirr, tolerance = 0.05)
    
    capture.output(print("Poisson ZI test= "),file = output_file_name, append = "TRUE")
    
    capture.output(print(pois_ZI),file = output_file_name, append = "TRUE")
    #check NB zero inflation
    
    NB_ZI <- check_zeroinflation(nb, tolerance = 0.05)
    
    capture.output(print("NB ZI test= "),file = output_file_name, append = "TRUE")
    capture.output(print(NB_ZI),file = output_file_name, append = "TRUE")
  }
  
  
  #close for loop -----
}

