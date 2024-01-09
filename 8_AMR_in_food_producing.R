#-------------8. AMR in food producing animals ------------------------------------

# This script outlines the data analysis for Chapter 8 of the 2023 BELMAP report, namely it:
# a)	Uploads the data
# b)	Cleans names (e.g. some in all caps etc.)
# c)	Per Host species / pathogen /indicator:
#   i.	Fits glm models with poisson and neg binom distributions
#   ii.	Checks for overdispersion, model fit and zero inflation, 
#   iii.	Compares these test values – if there is no overdispersion and the poisson model fit is better (by AIC and chi-square anova) select poisson model, if there is overdispersion and NB model fit is better select NB model, if there is disparity then flag a warning but select the model with best fit. 
#   iv.	Output message from zero inflation test and print all model and test outputs to a log file per host/pathogen/indicator for manual inspection, particularly if warning is flagged. – I no longer automatically run the zero-inflated models as the majority of the data contains no zeros- but the script includes the script to run and assess the ZI models, commented out, that can be run for specific cases if a warning of possible ZI is identified.
#   v.	Generate a dataset of predicted values based on the best fit model, which is used to plot the best fit model
#   vi.	Add Confidence intervals to this predicted dataset – based on the link function of the model (extract link function fit and standard error from best fit model, make columns of CI (fitted value +/- 2*SE based on link function),
#   vii.	Identify p-value and coefficient for time variable in best fit model: 
#             •	if model shows significant change over time then identify direction of coefficient – if positive then assign model label to be upward arrow plus number of * based on model p-value,  
#             •	if model shows significant change over time and direction of coefficient is negative then assign model label to be downward arrow plus number of * based on model p-value; 
#             •	if non-significant but the range (max-min value) is greater than 25% of mean then assign “oscillating” icon, 
#             •	if non-significant but range < 0.25% of mean then assign “equals” icon. 
#   viii.	Generate plot coordinates for model icons based on ranking and number of indicators per graph (to standardise location of model direction icons in figures)
#   ix.	Combine real dataset with predicted dataset and model labels
# d)	Make figures per pathogen, facetted by host (+/- primary/secondary indicator level for indicators such as E.coli), with barchart of real data and line of predicted model with ribbon for CI predicted from model link functions, add icons for models using font awesome HTML and showtext() package


#Contents-----
#1. Load libraries and themes
#2. Load data
#3. Run Models
#4. Make Figures

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
library(ggtext)
library(ggrepel)
library(showtext)
#load.fontawesome()
# install.packages("showtext")

BELMAP_colourscheme <- c("#1f78b4","#a6cee3","#b2df8a","#fb9a99","#33a02c","#e31a1c",
                                    "#fdbf6f","#ff7f00","#cab2d6", "#6a3d9a","#ffff99") #"MDR"="#e31a1c","pan-S" = "#33a02c",

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
#-----------------------------------------------------------------------------------
#2. Load data-----

Indicator_vet_data <- read_csv2("Data/2.AMR/Food_producing/Indicator-VET_data_combined.csv",
                                col_types = cols(
                                  "Year" = col_double(),
                                  "Sample_size"= col_double(),
                                  "Percent_resistance"= col_double(),
                                  "Number_resistants"= col_double(),
                                  "Pathogen" = col_character(),
                                  "Host" = col_character(),
                                  "Indicator" = col_character()
                                ),
                                na = c("", "NA"),
                                trim_ws = TRUE
                                ) %>%
  mutate(Year_simple = Year - 2011) %>%
  filter(!is.na(Sample_size)) 

unique(Indicator_vet_data$Indicator)
BELMAP_Indicator_colourscheme <- c("MDR"="#e31a1c",
                                    "pan-S" ="#33a02c",# ,
                                    "MRSA" = "#1f78b4",
                                    "ESBL producers" = "#a6cee3",
                                    "Ciprofloxacin resistance" = "#b2df8a",
                                    "Colistin resistance" = "#1f78b4")
                                    #"#33a02c","#e31a1c","#fdbf6f","#ff7f00","#cab2d6", "#6a3d9a","#ffff99") #

# Clean data

#correct for changes in methodology:
#MRSA - no models in pigs, but make same graph


#E.coli : , Amikacin was added to the monitoring in 2021 and thus the number of antibiotics considered for multi-drug resistance prevalence 
#was changed - indicate in text/Figure!




#-----------------------------------------------------------------------------------
#3. Run Models ----------------


# make category for each analysis - here one analysis for each antibiotic and source:
Indicator_vet_data_GLM_analysis <- Indicator_vet_data %>%
  mutate(Host = str_replace(string = Host, pattern = " \\(... isolation method .*\\)",replacement = ""))%>%
  mutate(category = paste(Host,Pathogen,Indicator,sep = "."))


#make list of categories

combined_categories<-unique(Indicator_vet_data_GLM_analysis$category)

#calc degree variance for categories----------

var_list = list()

for(i in combined_categories){
  print(i)
  dataset_analysis<- Indicator_vet_data_GLM_analysis %>%
    filter(category == i)
  degree_var= (max(dataset_analysis$Percent_resistance, na.rm = TRUE)-min(dataset_analysis$Percent_resistance, na.rm = TRUE))/mean(dataset_analysis$Percent_resistance, na.rm = TRUE)
  var_list[[i]] <-degree_var
  print(degree_var)
}


#make list for overdispersion plots
overdispersion_plot_list<- list()


#make csv for graph data with header

start_data_frame <- tibble(
  "Year" = "",
  "Pathogen" = "",
  "Host"= "",
  "Indicator"= "",
  "Percent_resistance"= "",
  "Percent_resistance_predict"= "",
  "CI_upper"= "",
  "CI_lower"= "",
  "Sample_size.x" = "",
  "icon" = "", 
  "signif_level" = ""
)

write_csv(start_data_frame,file = "8_food_producing_data_and_GLM_predictions.csv")

i="Beef Cattle.E. coli.Colistin resistance"
# make for loop to run analysis-----
for(i in combined_categories){
  
  #3.0.subset data-----
  
  dataset_analysis_raw<- Indicator_vet_data_GLM_analysis %>%
    filter(category == i)
  if(grepl("E. coli.ESBL producers",i)){
    dataset_analysis <- dataset_analysis_raw %>%  
      filter(Year >=2017) %>%
      mutate(Year_simple = Year - 2017)
  }else if(grepl("E\\. f",i)){
    dataset_analysis <- dataset_analysis_raw %>%  
      filter(Year >=2019) %>%
      mutate(Year_simple = Year - 2019)
  }else{
    dataset_analysis <- dataset_analysis_raw %>%  
      filter(Year >=2011)  
  }
  
  #make a file name to print each analysis
  
  output_file_name <- paste("8_food_producing_GLM_output_summaries/GLM_output_",i,".txt",sep="") 
  
  #print title of analysis
  capture.output(print(paste("GLM analysis of ",i,sep="")),file = output_file_name)
  
  #3.1. fit models:-----
  #3.1a. fit poisson-----
  
  glmpoissonirr <- glm(Number_resistants~Year_simple + offset(log(Sample_size)), data = dataset_analysis, family = poisson(link = "log"))
  
  
  #3.1b. fit NB-----
  
  nb <- glm.nb(Number_resistants~Year_simple + offset(log(Sample_size )), data = dataset_analysis )
 

  #--------------------------------------
  
  #3.2.check for overdispersion-----
  #visualise overdispersion 
  plot(log(fitted(glmpoissonirr)),log((dataset_analysis$Number_resistants-fitted(glmpoissonirr))^2),xlab=expression(hat(mu)),ylab=expression((y-hat(mu))^2),pch=20,col="blue")
  abline(0,1) ## 'variance = mean' line
  
  overdisperions_plot<- recordPlot()
  
  #assign to local environment
  assign(paste(i,"overdispersion_plot"),print(overdisperions_plot))
  
  #print dispersal graph to pdf
  dispersion_graph_file<-paste("8_food_producing_Dispersion_graphs/Dispersion_graph_poisson_model",i,".pdf",sep="") 
  pdf(dispersion_graph_file, width = 8, height = 6)
  # Call plot
  overdisperions_plot
  # Closing the graphical device
  dev.off() 
  
    # calculate dispersion parameter
  dp = sum(residuals(glmpoissonirr,type ="pearson")^2)/glmpoissonirr$df.residual

  #print over vs underdispersion

  
  underdispersed<- dispersiontest(glmpoissonirr,alternative = c("less") )
  
 
  #3.3 select model and make predict dataset----------------
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
  ndata<- tibble(
    Year = seq(2011,2022,by =1),
    Year_simple = seq(0,11,by = 1),
    Sample_size = rep(200,12),
    category = rep(i,12)
  )
  
  if(grepl("E. coli.ESBL producers",i)){
    ndata <- ndata %>%  
      filter(Year >=2017) %>%
      mutate(Year_simple = Year - 2017)
  }else if(grepl("E\\. f",i)){
    ndata <- ndata %>%  
      filter(Year >=2019) %>%
      mutate(Year_simple = Year - 2019)
  }else{
    ndata <- ndata %>%  
      filter(Year >=2011)  
  }
  
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
    if((max(dataset_analysis$Percent_resistance)-min(dataset_analysis$Percent_resistance))/mean(dataset_analysis$Percent_resistance) > 0.25){
      model_icon <- "oscilate"
      signif = ""
    }else if((max(dataset_analysis$Percent_resistance)-min(dataset_analysis$Percent_resistance))/mean(dataset_analysis$Percent_resistance) < 0.25){
      model_icon <- "equals"
      signif <- ""
      }
  }

  predict_data <- ndata %>%
    mutate(Percent_resistance_predict = fit/Sample_size*100,
           CI_upper = right_upr/Sample_size*100,
           CI_lower = right_lwr/Sample_size*100,
           icon = model_icon,
           signif_level = signif)
  
  
  # if non-signif - decide if equals or oscilating
  #calc if range > mean/4
  
  
  
  graph_data <- left_join(dataset_analysis,predict_data, by = c("Year","Year_simple")) %>%
    dplyr::select(Year,Pathogen, Host, Indicator,Percent_resistance,Percent_resistance_predict,CI_upper,CI_lower,Sample_size.x,
                  icon, signif_level)
  

  
  
  #save dataset
  
  predicted_dataset_filename <- paste("8_food_producing_GLM_predictions/GLM_predictions_",i,".csv",sep="") 

  write_csv(predict_data,file = predicted_dataset_filename)  
  write_csv(graph_data,file = "8_food_producing_data_and_GLM_predictions.csv", append = TRUE)
  
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
  

   
  #close for loop -----
}


#ERRORS - CHECK:
#"Beef Cattle.E. coli.MDR"
#"Veal Calves.E. coli.pan-S"


#Warnings poss ZI:
#"Pigs.E. coli.ESBL producers" poisson
# Poultry.E. coli.Colistin resistance poisson


--------------------------------------#In case zero inflated ---------------

# #3.4c. fit ZI NB-----
# 
# tryCatch({        #try catch so for loop doesn't fail if no zeros so ZI/underdispersed model doesn't fit
#   # glmCOM can handle underdispersion -----
#   
#   glmCOMpoisson <- glm.cmp(Number_resistants~Year_simple + offset(log(Sample_size )), data = dataset_analysis)
#   glmCOMpoisson
#   summary(glmCOMpoisson)
#   AIC(glmCOMpoisson, k = 2)
#   print(glmCOMpoisson)
#   coef(glmCOMpoisson)
#   nu(glmCOMpoisson)[1]
#   sdev(glmCOMpoisson)
#   vcov(glmCOMpoisson)
#   
#   # Likelihood ratio test for dispersion parameter
#   # Test for H_0: dispersion equal to 1 vs. H_1: not equal to 1
#   # (i.e. Poisson vs. COM-Poisson regression models)
#   lrt = equitest(glmCOMpoisson)
#   lrt
#   
#   capture.output(print("underdisperal model summary"),file = output_file_name, append = "TRUE")
#   
#   capture.output(print(summary(glmCOMpoisson)),file = output_file_name, append = "TRUE")
#   
#   capture.output(print("compare underdisperal model"),file = output_file_name, append = "TRUE")
#   
#   capture.output(print(lrt),file = output_file_name, append = "TRUE")
#   
#   
#   #3.4.d ZI models------
#   
#   modelZINB<-zeroinfl(Number_resistants ~ Year_simple + offset(log(Sample_size )) | Year_simple, data = dataset_analysis,  link = "logit", dist = "negbin")
#   #modelZINB1<-zeroinfl(Number_resistants ~ Year_simple, offset = log(Sample_size ), data = dataset_analysis,  link = "logit", dist = "negbin")
#   summary(modelZINB)
#   
#   
#   capture.output(print("Output from ZI NB model"),file = output_file_name, append = "TRUE")
#   capture.output(print(summary(modelZINB)),file = output_file_name, append = "TRUE")
#   
#   
#   #3.4.e. fit ZI poisson -----
#   #zero-inflated Poisson model #logit? is for zero-inflated - two groups
#   #The zero-inflated Poisson model assumes that the population consists of two ‘latent’ (or unobserved) groups. The first group consists 
#   #of observational units who always have zero counts. The first group is referred to as the ‘excess’ zeroes group (inflation component). Second group is the count component.
#   #include the offset in the count component, but not in the zero inflation component
#   
#   modelZIP<-zeroinfl(Number_resistants~Year_simple + offset(log(Sample_size )) | Year_simple, data = dataset_analysis, link = "logit", dist = "poisson")
#   #modelZIP1<-zeroinfl(Number_resistants~Year_simple, offset = log(Sample_size ), data = dataset_analysis, link = "logit", dist = "poisson")
#   summary(modelZIP)
#   AIC(modelZIP)
#   
#   capture.output(print("Output from ZI Poisson model"),file = output_file_name, append = "TRUE")
#   capture.output(print(summary(modelZIP)),file = output_file_name, append = "TRUE")
#   
#   #3.4.f.compare ZI models AIC ------
#   #3.4.f. ZI NB model AIC -----
#   
#   AIC_ZINB <- AIC(modelZINB)
#   
#   capture.output(print(paste("AIC of ZI NB model", AIC_ZINB)),file = output_file_name, append = "TRUE")
#   
#   #3.4.f ZI poisson model AIC -----
#   AIC_ZI_poisson <- AIC(modelZIP)
#   
#   capture.output(print(paste("AIC of ZI Poisson model", AIC_ZI_poisson)),file = output_file_name, append = "TRUE")
#   
#   
#   #3.5.a check NB vs ZI NB---------
#   
#   #NB versus zero-inflated NB
#   vuong(nb,modelZINB)
#   
#   
#   NB_vs_ZINB<- vuong(nb,modelZINB)
#   
#   capture.output(print("Vuong Test NB vs ZI NB model   "),file = output_file_name, append = "TRUE")
#   capture.output(print(NB_vs_ZINB),file = output_file_name, append = "TRUE")  
#   
#   #3.5.b check poisson vs ZI poisson---------
#   vuong(glmpoissonirr,modelZIP)
#   
#   pois_vd_ZIpois<- vuong(glmpoissonirr,modelZIP)
#   
#   capture.output(print("Vuong Test Poisson vs ZI Poisson model   "),file = output_file_name, append = "TRUE")
#   capture.output(print(pois_vd_ZIpois),file = output_file_name, append = "TRUE")
#   
#   #3.5.c likelihood ratio test to compare Zero-inflated Poisson with Zero-inflated NB -----
#   
#   lrtest(modelZIP,modelZINB)
#   
#   ZIpois_vd_ZINB<- lrtest(modelZIP,modelZINB)
#   
#   capture.output(print("LR Test ZI Poisson vs NB model   "),file = output_file_name, append = "TRUE")
#   capture.output(print(ZIpois_vd_ZINB),file = output_file_name, append = "TRUE")
#   
#   
#   
#   #3.5.d likelihood ratio test to compare Zero-inflated NB with Zero-inflated Poisson -----
#   
#   lrtest(modelZINB,modelZIP)
#   
#   ZINB_vd_ZIPois<- lrtest(modelZINB,modelZIP)
#   
#   capture.output(print("LR Test ZI NB vs ZI Poisson model   "),file = output_file_name, append = "TRUE")
#   capture.output(print(ZINB_vd_ZIPois),file = output_file_name, append = "TRUE")
#   
#   
# }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
# 



#-----------------------------------------------------------------------------------
#manually inspect output from GLMs - identify best model and make table of best models

# Models_directions <- tibble(
#   #  Combined_category = combined_categories
#   Combined_category = c("Poultry.E. coli.pan-S","Poultry.E. coli.MDR" ,"Poultry.E. coli.ESBL producers", "Poultry.E. coli.Ciprofloxacin resistance","Poultry.E. coli.Colistin resistance")
# ) %>%
#   mutate( Model_direction = case_when(
#     Combined_category == "Poultry.E. coli.pan-S"~ "",
#     Combined_category == "Poultry.E. coli.MDR" ~ "",
#     Combined_category == "Poultry.E. coli.ESBL producers" ~ "",
#     Combined_category == "Poultry.E. coli.Ciprofloxacin resistance"~ "",
#     Combined_category == "Poultry.E. coli.Colistin resistance"~ "",
#     # Combined_category == "Beef Cattle.E. coli.pan-S" ~ "",
#     # Combined_category == "Beef Cattle.E. coli.MDR"~ "",
#     # Combined_category == "Beef Cattle.E. coli.ESBL producers" ~ "",
#     # Combined_category == "Beef Cattle.E. coli.Ciprofloxacin resistance"~ "",
#     # Combined_category == "Beef Cattle.E. coli.Colistin resistance" ~ "",
#     # Combined_category == "Pigs.E. coli.pan-S" ~ "",
#     # Combined_category == "Pigs.E. coli.MDR" ~ "",
#     # Combined_category == "Pigs.E. coli.ESBL producers" ~ "",
#     # Combined_category == "Pigs.E. coli.Ciprofloxacin resistance" ~ "",
#     # Combined_category == "Pigs.E. coli.Colistin resistance" ~ "",
#     # Combined_category == "Veal Calves.E. coli.pan-S" ~ "",
#     # Combined_category == "Veal Calves.E. coli.MDR"~ "",
#     # Combined_category == "Veal Calves.E. coli.ESBL producers" ~ "",
#     # Combined_category == "Veal Calves.E. coli.Ciprofloxacin resistance"~ "",
#     # Combined_category == "Veal Calves.E. coli.Colistin resistance"~ ""
#   ))
# 


#-----------------------------------------------------------------------------------

#4. Make Figures------------------------------

# adding icons to figures ------------------------------------

font_add('fa-solid', '../font_awesome_font_files/fontawesome-free-6.4.0-desktop/otfs/Font Awesome 6 Free-Solid-900.otf')

upward_arrow <- "<span style='font-family:fa-solid'>&#xf062;</span>"
downward_arrow <- "<span style='font-family:fa-solid'>&#xf063;</span>" 
equals <- "<span style='font-family:fa-solid'>&#xf52c;</span>"   #" = " 
oscillate <- "<span style='font-family:fa-solid'>&#xf83e;</span>"   #" ~ "



#load graph data --------------------------------------------------

graph_data_complete <- read.csv("8_food_producing_data_and_GLM_predictions.csv") %>%
  mutate(label_icon = case_when(
  icon == "upward_arrow" ~ upward_arrow,
  icon == "downward_arrow" ~ downward_arrow,
  icon == "equals" ~ equals,
  icon == "oscilate" ~ oscillate
)) %>%
  mutate(label = if_else(signif_level == "", label_icon, paste(label_icon,signif_level,sep=" ")))


y_coords_oth <- graph_data_complete %>%
  mutate(Level_indicator = if_else(grepl("MDR|pan-S",Indicator), "Primary","Secondary")) %>%
  filter(!grepl("MRSA",Pathogen))%>%
  filter(Year == 2022) %>%
  group_by(Host,Pathogen,Level_indicator)%>%
  mutate(rank = rank(Percent_resistance_predict))%>%
  ungroup() %>%
  mutate(y_coord = case_when(
    rank == 1 ~ 10,
    rank == 2 ~ 30,
    rank == 3 ~ 50
  )) %>%
  dplyr::select(Host,Pathogen,Indicator,y_coord)

y_coords_mr <- graph_data_complete %>%
  mutate(Level_indicator = if_else(grepl("MDR|pan-S",Indicator), "Primary","Secondary")) %>%
  filter(grepl("MRSA",Pathogen))%>%
  filter(grepl("2020|2021|2022",Year)) %>%
  filter(!is.na(Percent_resistance_predict))%>%
  group_by(Host,Pathogen,Level_indicator)%>%
  mutate(rank = rank(Percent_resistance_predict))%>%
  ungroup() %>%
  mutate(y_coord = case_when(
    rank == 1 ~ 10,
    rank == 2 ~ 30,
    rank == 3 ~ 50
  )) %>%
dplyr::select(Host,Pathogen,Indicator,y_coord)


y_coords<-rbind(y_coords_mr,y_coords_oth) 
  

graph_data_complete_labels <- left_join(graph_data_complete,y_coords, by = c("Host","Pathogen","Indicator"))



#4.1 E.coli -----------------------------------

E.coli_sample_size <- Indicator_vet_data %>%
    filter(Pathogen == "E_coli") %>%
    group_by(Year,Host,Pathogen) %>%
    summarise(sample_size = median(as.numeric(Sample_size), na.rm=TRUE)) %>%
    ungroup()

(E_coli_graph <- graph_data_complete_labels %>%
  filter(Pathogen == "E. coli") %>%
    mutate(Level_indicator = if_else(grepl("MDR|pan-S",Indicator), "Primary","Secondary")) %>%
    ggplot(aes(label = label_icon))+
    geom_bar(aes(x=Year, y = as.numeric(Percent_resistance), fill = Indicator), stat="identity", position = "dodge")+
    geom_ribbon(aes(x=Year, fill = Indicator, ymin = CI_lower, ymax = CI_upper, y = Percent_resistance_predict), alpha = 0.25)+
    geom_smooth(aes(x=Year, colour = Indicator, y = Percent_resistance_predict),stat = "identity")+
    scale_fill_manual(values = BELMAP_Indicator_colourscheme,breaks=c("MDR", "pan-S", "Ciprofloxacin resistance", "Colistin resistance", "ESBL producers"))+
    scale_colour_manual(values = BELMAP_Indicator_colourscheme,breaks=c("MDR", "pan-S", "Ciprofloxacin resistance", "Colistin resistance", "ESBL producers"))+
    facet_grid(Level_indicator ~ Host)+
   # geom_richtext( size = 16, hjust = 0, label.colour = NA) +
    scale_x_continuous(limits = c(2010,2024), breaks = seq(2011,2022,1))+
    scale_y_continuous(limits = c(-5,100), breaks = seq(0,100,50))+
    geom_richtext(aes(x = 2023.5, y = y_coord, 
                      label = label, 
                      fill = Indicator), stat = "unique", 
                      colour= "white", show.legend = FALSE)+
    # fill = NA, label.color = NA, # remove background and outline
    # label.padding = grid::unit(rep(0, 4), "pt")) +
  #  geom_text(aes(x = Year, y = 2, label = Sample_size.x), na.rm = TRUE, size = 2.5)+  #make this additional for interactive report
    labs(y= "% Resistance")&
    moiras_graph_theme())   


# 4.2 MRSA -------------------------------------------
graph_data_MRSA <- graph_data_complete_labels %>%
  filter(Pathogen == "MRSA") %>%
  mutate(label = if_else(grepl("pigs|Sow",Host),"",label))

(MRSA_graph <- graph_data_MRSA %>%
   ggplot()+
   geom_bar(aes(x=Year, y = as.numeric(Percent_resistance), fill = Indicator), stat="identity", position = "dodge")+
   geom_vline(data = filter(graph_data_MRSA, grepl("pigs|Sow",Host)),aes(xintercept =2021), colour = "red", linetype= "dashed")+
   geom_ribbon(data = filter(graph_data_MRSA, !grepl("pigs|Sow",Host)),aes(x=Year, fill = Indicator, ymin = CI_lower, ymax = CI_upper, y = Percent_resistance_predict), alpha = 0.1)+
   geom_smooth(data = filter(graph_data_MRSA, !grepl("pigs|Sow",Host)),aes(x=Year, colour = Indicator, y = Percent_resistance_predict),stat = "identity")+
   scale_fill_manual(values = BELMAP_colourscheme)+
   scale_colour_manual(values = BELMAP_colourscheme)+
    geom_richtext(aes(x = 2023.5, y = y_coord, 
                      label = label, 
                      fill = Indicator), stat = "unique", 
                  colour= "white", show.legend = FALSE)+
                  # fill = NA, label.color = NA, # remove background and outline
                  # label.padding = grid::unit(rep(0, 4), "pt")) +
    facet_wrap(~ Host, nrow = 2)+
   scale_x_continuous(limits = c(2010,2024), breaks = seq(2011,2022,1))+
   scale_y_continuous(limits = c(-5,100), breaks = seq(0,100,50))+
 #  geom_text(aes(x = Year, y = 2, label = Sample_size.x), na.rm = TRUE, size = 2.5)+
   labs(y= "% Resistance")+
   moiras_graph_theme()
)   


# 4.3 Enterococci ------------------------------------

# graph_data_enterococci <- graph_data_complete_labels %>%
#   filter(grepl("E. f",Pathogen))
# 
# (Enterococci_graph <- graph_data_enterococci %>%
#     ggplot()+
#     geom_bar(aes(x=Year, y = as.numeric(Percent_resistance), fill = Indicator), stat="identity", position = "dodge")+
#     geom_ribbon(aes(x=Year, fill = Indicator, ymin = CI_lower, ymax = CI_upper, y = Percent_resistance_predict), alpha = 0.1)+
#     geom_smooth(aes(x=Year, colour = Indicator, y = Percent_resistance_predict),stat = "identity")+
#     scale_fill_manual(values = BELMAP_colourscheme)+
#     scale_colour_manual(values = BELMAP_colourscheme)+
#     facet_wrap(Pathogen ~ Host, nrow = 2)+
#     geom_richtext(aes(x = 2023.5, y = y_coord, 
#                       label = label, 
#                       colour = Indicator),
#                   fill = NA, label.color = NA, # remove background and outline
#                   label.padding = grid::unit(rep(0, 4), "pt")) +
#     scale_x_continuous(limits = c(2018,2025), breaks = seq(2019,2022,1))+
#     scale_y_continuous(limits = c(-5,100), breaks = seq(0,100,50))+
#  #   geom_text(aes(x = Year, y = 2, label = Sample_size.x), na.rm = TRUE, size = 2.5)+
#     labs(y= "% Resistance")+
#     moiras_graph_theme()
# )   
# 



(E_faecalis_graph <- graph_data_complete_labels %>%
    filter(grepl("E. faecalis",Pathogen)) %>%
    ggplot()+
    geom_bar(aes(x=Year, y = as.numeric(Percent_resistance), fill = Indicator), stat="identity", position = "dodge")+
    geom_ribbon(aes(x=Year, fill = Indicator, ymin = CI_lower, ymax = CI_upper, y = Percent_resistance_predict), alpha = 0.1)+
    geom_smooth(aes(x=Year, colour = Indicator, y = Percent_resistance_predict),stat = "identity")+
    scale_fill_manual(values = BELMAP_colourscheme)+
    scale_colour_manual(values = BELMAP_colourscheme)+
    facet_wrap(~Host, nrow = 1)+
    scale_x_continuous(limits = c(2010,2024), breaks = seq(2011,2022,1))+
    scale_y_continuous(limits = c(-5,100), breaks = seq(0,100,50))+
   geom_richtext(aes(x = 2023.5, y = y_coord, 
                     label = label, 
                     fill = Indicator), stat = "unique", 
                 colour= "white", show.legend = FALSE)+
 #   geom_text(aes(x = Year, y = 2, label = Sample_size.x), na.rm = TRUE, size = 2.5)+
    labs(y= "% Resistance")+
    moiras_graph_theme()
) 

(E_faecium_graph <- graph_data_complete %>%
    filter(grepl("E. faecium",Pathogen)) %>%
    ggplot()+
    geom_bar(aes(x=Year, y = as.numeric(Percent_resistance), fill = Indicator), stat="identity", position = "dodge")+
    geom_ribbon(aes(x=Year, fill = Indicator, ymin = CI_lower, ymax = CI_upper, y = Percent_resistance_predict), alpha = 0.1)+
    geom_smooth(aes(x=Year, colour = Indicator, y = Percent_resistance_predict),stat = "identity")+
    scale_fill_manual(values = BELMAP_colourscheme)+
    scale_colour_manual(values = BELMAP_colourscheme)+
    geom_richtext(aes(x = 2023.5, y = y_coord, 
                      label = label, 
                      fill = Indicator), stat = "unique", 
                  colour= "white", show.legend = FALSE)+
    facet_wrap(~Host, nrow = 1)+
    scale_x_continuous(limits = c(2018,2024), breaks = seq(2019,2022,1))+
    scale_y_continuous(limits = c(-5,100), breaks = seq(0,100,50))+
  #  geom_text(aes(x = Year, y = 2, label = Sample_size.x), na.rm = TRUE, size = 2.5)+
    labs(y= "% Resistance")+
    moiras_graph_theme()
)   


#NOTES and UNUSED CODE

font_paths()
font_files()

