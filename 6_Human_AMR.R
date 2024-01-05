#--------------6. AMR in Human Pathogens-------------------------------------------------------

# This script outlines the data analysis for Chapter 6 of the 2023 BELMAP report, namely it:
# a)	Uploads the data
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
# d)	Make figures per pathogen, facetted by host (+/- primary/secondary indicator level for indicators such as E.coli), with barchart of real data and line of predicted model with ribbon for CI predicted from model link functions, add icons for models using font awesome HTML and showtext() package
#Contents-----
#1. Load libraries and themes
#2. make ECDC datasets
#3. list pathogens and indicators - make 
#4. per pathogen: Load data and combine with ECDC
#5. Run Models
#6. Make Figures

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


#2. make ECDC datasets-------------------------------------------------


##load ecdc data
ECDC_raw <-  read_csv2("Data/2.AMR/Human_AMR/ECDC_surveillance_data_Antimicrobial_resistance_num1.csv",
                       #col_types = cols(
                       #  "NumValue" = col_double()#,
                       # "Sample_size"= col_double(),
                       # "Percent_resistant"= col_double(),
                       # "Number_resistants"= col_double(),
                       # "Pathogen" = col_character(),
                       # "Host" = col_character(),
                       # "Indicator" = col_character()
                       #),
                       na = c("", "NA", "-"),
                       trim_ws = TRUE
)

ECDC_data <- distinct(ECDC_raw) %>%
  separate(col = Population, into = c("Pathogen", "Antibiotic"), sep = "\\|", remove = TRUE) %>%
  dplyr::select(-c(Unit,TxtValue,RegionCode)) %>%
  pivot_wider(names_from = Indicator, values_from = NumValue) %>%
  mutate(
    Year = Time,
    Sample_size = `Total tested isolates`,
    Percent_resistant = `R - resistant isolates, percentage`,
    Number_resistants = `R - resistant isolates`,
    Region = RegionName
  ) %>%
  dplyr::select(-c(`I - 'susceptible, increased exposure' isolates`,`S - susceptible isolates`,
                   `Completeness age`,`Completeness gender`,`Penicillin non-wild-type isolates, percentage`,
                   Time, `R - resistant isolates`, `R - resistant isolates, percentage`,
                   `Total tested isolates`, RegionName, HealthTopic))

# make EU dataframe

ECDC_data_EU <- ECDC_data %>%
  group_by(Pathogen,Antibiotic,Year) %>%
  summarise(Sample_size = sum(as.numeric(Sample_size, na.rm = T)),
            Number_resistants = sum(Number_resistants, na.rm = T))%>%
  mutate(Percent_resistant = Number_resistants/Sample_size*100,
         Region = "Europe",
         Surveillance = "EARS",
         Matrix = "Blood_CSF")
# make neighbours dataframe

ECDC_data_neighbours <- ECDC_data %>%
  filter(grepl("Germany|France|Luxembourg|Netherlands",Region))%>%
  filter(!(Antibiotic == "Data Quality")) %>%
  group_by(Pathogen,Antibiotic,Year) %>%
  summarise(Sample_size = sum(as.numeric(Sample_size, na.rm = T)),
            Number_resistants = sum(Number_resistants, na.rm = T))%>%
  mutate(Percent_resistant = Number_resistants/Sample_size*100,
         Region = "Neighbours",
         Surveillance = "EARS",
         Matrix = "Blood_CSF")
# combine
ECDC_data_combined <- rbind(ECDC_data_EU,ECDC_data_neighbours)


#3. list pathogens and indicators --------------------------------------------------

#list pathogen folders

list.files("Data/2.AMR/Human_AMR/",
           #     pattern=".pdf", all.files=TRUE,
           full.names=TRUE)
# [1] "Acinetobacter"                                          "Aspergillus"
# [3] "C.difficile"                                         "Candida"
# [5] "E.coli"                                              "E.faecalis"
# [7] "E.faecium"                                           "ECDC_surveillance_data_Antimicrobial_resistance.csv"
# [9] "H.influenzae"                                        "H.pylori"
# [11] "K.pneumoniae"                                        "M.genitalum"
# [13] "M.tuberculosis"                                      "MRSA"
# [15] "N.gonorrhea"                                         "P.aeruginosa"
# [17] "S.pneumoniae"                                        "Salmonella_invasive"
# [19] "Shigella"

# pathogen_list<- grep(list.files(path="Data/2.AMR/Human_AMR/"),
#      pattern='ECDC', invert=TRUE, value=TRUE) # list pathogens excluding ecdc file
# 
# 
# #pathogen_list <-pathogen_list[4:18]
# 
# data_file_list <- list()
# for (i in pathogen_list) {
#   data_file_list[[i]] <-  list.files(paste("Data/2.AMR/Human_AMR/",i,"/",sep=""),
#                                pattern=".csv", all.files=TRUE,
#                                full.names=TRUE)
# }
# 
# 
# #4. per pathogen: Load data and combine with ECDC -------------------------------------------------
# 
# # load data per pathogen
# Human_AMR_dataframe<- tibble(
#   Year = "",
#   Pathogen = "",
#   Antibiotic = "",
#   Sample_size = "",
#   Percent_resistant = "",
#   Number_resistants= "",
#   Region = "",
#   Surveillance = "",
#   Matrix = ""
# )
# 
# for (i in pathogen_list) {
#   pathogen_data <- read_csv2(file = data_file_list[[i]],
#                              na = c("", "NA", "na"),
#                              trim_ws = TRUE
#   )  %>%
#     dplyr::select(Year, Pathogen,Antibiotic, Sample_size, Percent_resistant,
#            Number_resistants,     Region, Surveillance,Matrix) %>%
#       mutate(Percent_resistant = as.numeric(str_replace(Percent_resistant,",",".")),
#              Sample_size = as.numeric(str_replace(Sample_size,",",".")),
#              Number_resistants = as.numeric(str_replace(Number_resistants,",",".")),
#              Matrix = str_replace(Matrix,"/","_")) #get rid of / as messes up name of files
# 
#   Human_AMR_dataframe <- rbind(Human_AMR_dataframe, pathogen_data)
#     }
# 
# 
#  #write human pathogen dataframe
#  write.csv(Human_AMR_dataframe, "6_Human_AMR_dataframe_raw3.csv")


Human_AMR_dataframe<- read.csv("6_Human_AMR_dataframe_raw3.csv",dec = ",", header=TRUE ) %>%
  dplyr::select(-X)

## match pathogens our to ecdc database
# unique(Human_AMR_dataframe$Pathogen)
#
# unique(ECDC_data_combined$Pathogen)
#


ECDC_data_combined <- ECDC_data_combined %>%
  mutate(Pathogen1 = Pathogen) %>%
  mutate(Pathogen = case_when(
    Pathogen1 == "Acinetobacter spp." ~ "Acinetobacter",
    Pathogen1 == "Enterococcus faecalis" ~ "E. faecalis",
    Pathogen1 == "Enterococcus faecium" ~ "E. faecium",
    Pathogen1 == "Escherichia coli" ~ "E.coli",
    Pathogen1 == "Klebsiella pneumoniae" ~ "K.pneumoniae",
    Pathogen1 == "Pseudomonas aeruginosa" ~ "P. aeruginosa",
    Pathogen1 == "Staphylococcus aureus" ~ "S. aureus",
    Pathogen1 == "Streptococcus pneumoniae" ~ "S. pneumoniae"
  )) %>%
  dplyr::select(-Pathogen1)

# ## match antiBs our to ecdc database
# unique(Human_AMR_dataframe[Human_AMR_dataframe$Region == "Europe",]$Antibiotic)
# unique(Human_AMR_dataframe$Antibiotic)
#unique(ECDC_data_combined$Antibiotic)



Human_AMR_dataframe <- Human_AMR_dataframe %>%
  filter(!is.na(Antibiotic))%>%
  mutate(Antibiotic1 = Antibiotic) %>%
  mutate(Antibiotic = case_when(
    Antibiotic1 == "carbapenem" ~ "Carbapenem",
    Antibiotic1 == "Carbepenem" ~ "Carbapenem",
    Antibiotic1 == "carbepenem" ~ "Carbapenem",
    Antibiotic1 == "carbapenem(m\xe9ropenem)"~ "Carbapenem",
    Antibiotic1 == "carbapenem (méropenem) "~ "Carbapenem",
    Antibiotic1 == "carbapenem (méropenem)"~ "Carbapenem",
    Antibiotic1 == "colistin" ~ "Colistin",
    Antibiotic1 == "3CG" ~ "3GC",
    Antibiotic1 == "ciprofloxacin" ~ "Ciprofloxacin",
    Antibiotic1 == "Vancomycin" ~ "Vancomycin",
    Antibiotic1 == "Ag, Fq And 3gc" ~ "AG, FQ And 3GC",
    Antibiotic1 == "AG, FQ and 3GC" ~ "AG, FQ And 3GC",
    Antibiotic1 == "Methicillin" ~ "Methicillin",
    Antibiotic1 == "Macrolide" ~ "Macrolide",
    Antibiotic1 == "Penicillin" ~ "Penicillin",
    Antibiotic1 == "Mdr" ~ "MDR",
    .default = Antibiotic1
  )) %>%
  dplyr::select(-Antibiotic1)

#unique(ECDC_data_combined$Antibiotic)
ECDC_data_combined <- ECDC_data_combined %>%
  mutate(Antibiotic1 = Antibiotic) %>%
  mutate(Antibiotic = case_when(
    Antibiotic1 == "Carbapenems" ~ "Carbapenem",
    Antibiotic1 == "Third-generation cephalosporins" ~ "3GC",
    Antibiotic1 == "Fluoroquinolones" ~ "Ciprofloxacin",
    Antibiotic1 == "Vancomycin" ~ "vancomycin",
    Antibiotic1 == "Combined resistance (fluoroquinolones, aminoglycosides and carbapenems)" ~ "AG, FQ And 3GC",
    Antibiotic1 == "Combined resistance (third-generation cephalosporin, fluoroquinolones and aminoglycoside)"~ "AG, FQ And 3GC",
    Antibiotic1 == "Meticillin (MRSA)" ~ "Methicillin",
    Antibiotic1 == "Aminopenicillins" ~ "Aminopenicillins",
    Antibiotic1 == "Ceftazidime" ~ "3GC",
    Antibiotic1 == "Macrolides" ~ "Macrolide",
    Antibiotic1 == "Aminoglycosides" ~ "Aminoglycosides",
    Antibiotic1 == "High-level gentamicin" ~ "High-level gentamicin",
    Antibiotic1 == "Penicillins" ~ "Penicillin",
    Antibiotic1 == "PiperacillinTazobactam"~ "Piperacillin-Tazobactam",
    Antibiotic1 == "Combined resistance (at least three of piperac. and tazob., fluoroq., ceftaz., aminogl. and carbapenems)" ~ "MDR"
  )) %>%
  filter(!is.na(Antibiotic)) %>%
  dplyr::select(-c(Antibiotic1)) %>%
  filter(Year > 2010)
#combine datasheets

Human_AMR_Belgium_EU <- rbind(Human_AMR_dataframe, ECDC_data_combined) %>%
  filter(!is.na(Sample_size)) %>%
  filter(Sample_size !="") %>%
  filter(!grepl("^na$",Sample_size)) %>%
  mutate(Surveillance = dplyr::case_when(
    Region == "Europe" ~ "EARS",
    Region == "Neighbours" ~ "EARS",
    is.na(Surveillance) ~ "NRC",
    .default = Surveillance
  )) %>%
  mutate(Matrix = dplyr::case_when(
    Region == "Europe" ~ "Blood_CSF",
    Region == "Neighbours" ~ "Blood_CSF",
    (is.na(Matrix) & Surveillance == "NRC") ~ "NRC_samples",
    .default = Matrix
  )) %>%
  mutate(category = paste(Pathogen,Antibiotic, Region,Matrix, sep = "_"),
         Number_resistants = as.numeric(Number_resistants),
         Sample_size = as.numeric(Sample_size),
         Percent_resistant = as.numeric(Percent_resistant))

Human_AMR_Belgium_EU1 <- Human_AMR_Belgium_EU %>%
  mutate(#Percent_resistant = as.numeric(Percent_resistant),
    Sample_size  = as.numeric(Sample_size),
    Number_resistants = as.numeric(Number_resistants),
    Matrix = if_else(grepl("Europe|Neighbour", Region), "Blood_CSF",Matrix),
    Surveillance = if_else(grepl("Europe|Neighbour", Region), "EARS",Surveillance)
  )

#write.csv(Human_AMR_Belgium_EU, "6_Human_AMR_Belgium_EU_raw3.csv")


# load processed data:
?read_csv()
Human_AMR_Belgium_EU <- read_delim( "6_Human_AMR_Belgium_EU_raw3.csv", delim= ",")# %>%
# dplyr::select(-`...1`) %>%
mutate(Percent_resistant = (Number_resistants/Sample_size)*100)

#5. Run Models -------------------------------------------------

#make list of categories

combined_categories<-unique(Human_AMR_Belgium_EU$category)
combined_categories1 <- c("S. aureus_Methicillin_Belgium_All samples",
                          "E.coli_3GC_Belgium_All samples",
                          "E.coli_Carbapenem_Belgium_All samples")
#calc degree variance for categories----------

var_list = list()

# i ="E. faecium_vancomycin_Neighbours_Blood_CSF"
# combined_categories1<- c("Acinetobacter_Carbapenem_Belgium_Blood_CSF"#, "Acinetobacter_Carbapenem_Europe_Blood_CSF",
#                          #"Acinetobacter_Carbapenem_Neighbours_Blood_CSF"
#                         )


for(i in combined_categories1){
  print(i)
  dataset_analysis<- Human_AMR_Belgium_EU %>%
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
  "Pathogen" = "",
  "Antibiotic"= "",
  "Percent_resistant"= "",
  "Region" = "",
  "Percent_resistant_predict"= "",
  "Surveillance" = "",
  "Matrix"="",
  "CI_upper"= "",
  "CI_lower"= "",
  "Sample_size.x" = "",
  "icon" = "",
  "signif_level" = ""
)



#write_csv(start_data_frame,file = "6_Human_AMR_data_and_GLM_predictions.csv")

# make for loop to run analysis-----
for(i in combined_categories1){
  
  #3.0.subset data-----
  
  dataset_analysis_raw<- Human_AMR_Belgium_EU %>%
    filter(category == i) %>%
    mutate(Year = as.numeric(Year))
  min_year = min(dataset_analysis_raw$Year)
  max_year = max(dataset_analysis_raw$Year)
  dataset_analysis <- dataset_analysis_raw %>%  
    mutate(Year_simple = Year - min_year) %>%
    distinct()
  
  #make a file name to print each analysis
  
  output_file_name <- paste("6_Human_AMR_GLM_output_summaries/GLM_output_",i,".txt",sep="") 
  
  #print title of analysis
  capture.output(print(paste("GLM analysis of ",i,sep="")),file = output_file_name)
  
  if(max_year-min_year < 3){
    capture.output(print("Too few data points for model"),file = output_file_name)
    graph_data <- dataset_analysis %>%
      mutate(Percent_resistant_predict = "",
             CI_upper = "",
             CI_lower = "",
             Sample_size.x = "Sample_size",
             icon = "",
             signif_level = "") %>%
      dplyr::select(Year,Pathogen, Antibiotic,Percent_resistant,Region, Percent_resistant_predict,Surveillance,Matrix,CI_upper,CI_lower,Sample_size.x,
                    icon, signif_level)
  } else if(sum(as.numeric(dataset_analysis$Percent_resistant)) == 0){
    capture.output(print("No resistance- can't fit model"),file = output_file_name)
    graph_data <- dataset_analysis %>%
      mutate(Percent_resistant_predict = "",
             CI_upper = "",
             CI_lower = "",
             Sample_size.x = "Sample_size",
             icon = "",
             signif_level = "") %>%
      dplyr::select(Year,Pathogen, Antibiotic,Percent_resistant,Region, Percent_resistant_predict,Surveillance,Matrix,CI_upper,CI_lower,Sample_size.x,
                    icon, signif_level)
  } else {
    
    #3.1. fit models:-----
    #3.1a. fit poisson-----
    
    glmpoissonirr <- glm(Number_resistants~Year_simple + offset(log(Sample_size)), data = dataset_analysis, family = poisson(link = "log"))
    
    
    #3.1b. fit NB-----
    
    nb <- glm.nb(Number_resistants~Year_simple + offset(log(Sample_size )), data = dataset_analysis )
    
    
    #--------------------------------------
    
    # #3.2.check for overdispersion-----
    # #visualise overdispersion 
    # plot(log(fitted(glmpoissonirr)),log((dataset_analysis$Number_resistants-fitted(glmpoissonirr))^2),xlab=expression(hat(mu)),ylab=expression((y-hat(mu))^2),pch=20,col="blue")
    # abline(0,1) ## 'variance = mean' line
    # 
    # overdisperions_plot<- recordPlot()
    # 
    # #assign to local environment
    # assign(paste(i,"overdispersion_plot"),print(overdisperions_plot))
    # 
    # #print dispersal graph to pdf
    # dispersion_graph_file<-paste("6_Human_AMR_Dispersion_graphs/Dispersion_graph_poisson_model",i,".pdf",sep="") 
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
      dplyr::select(Year,Pathogen, Antibiotic,Percent_resistant,Region, Percent_resistant_predict,Surveillance,Matrix,CI_upper,CI_lower,Sample_size.x,
                    icon, signif_level)
    
    
    
    
    #save dataset
    
    predicted_dataset_filename <- paste("6_Human_AMR_GLM_predictions/GLM_predictions_",i,".csv",sep="") 
    
    write_csv(predict_data,file = predicted_dataset_filename)  
    write_csv(graph_data,file = "6_Human_AMR_data_and_GLM_predictions.csv", append = TRUE)
    
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




#   #6. Make Figures  -------------------------------------------------
#   
#   # adding icons to figures ------------------------------------
#   
#   font_add('fa-solid', '../font_awesome_font_files/fontawesome-free-6.4.0-desktop/otfs/Font Awesome 6 Free-Solid-900.otf')
#   
#   upward_arrow <- "<span style='font-family:fa-solid'>&#xf062;</span>"
#   downward_arrow <- "<span style='font-family:fa-solid'>&#xf063;</span>" 
#   equals <- "<span style='font-family:fa-solid'>&#xf52c;</span>"   #" = " 
#   oscillate <- "<span style='font-family:fa-solid'>&#xf83e;</span>"   #" ~ "
# 
#   
#   graph_data_complete <- graph_data %>%
#     mutate(label_icon = case_when(
#       icon == "upward_arrow" ~ upward_arrow,
#       icon == "downward_arrow" ~ downward_arrow,
#       icon == "equals" ~ equals,
#       icon == "oscilate" ~ oscillate
#     )) %>%
#     mutate(label = if_else(signif_level == "", label_icon, paste(label_icon,signif_level,sep=" ")))
#   
#   
#   # make y coordinates for icons
#   
#   
#   
#   y_coords_oth <- graph_data_complete %>%
#     mutate(Level_indicator = if_else(grepl("MDR|pan-S",Indicator), "Primary","Secondary")) %>%
#     filter(!grepl("MRSA",Pathogen))%>%
#     filter(Year == 2022) %>%
#     group_by(Host,Pathogen,Level_indicator)%>%
#     mutate(rank = rank(Percent_resistant_predict))%>%
#     ungroup() %>%
#     mutate(y_coord = case_when(
#       rank == 1 ~ 10,
#       rank == 2 ~ 30,
#       rank == 3 ~ 50
#     )) %>%
#     dplyr::select(Host,Pathogen,Indicator,y_coord)
#   
#   y_coords_mr <- graph_data_complete %>%
#     mutate(Level_indicator = if_else(grepl("MDR|pan-S",Indicator), "Primary","Secondary")) %>%
#     filter(grepl("MRSA",Pathogen))%>%
#     filter(grepl("2020|2021|2022",Year)) %>%
#     filter(!is.na(Percent_resistant_predict))%>%
#     group_by(Host,Pathogen,Level_indicator)%>%
#     mutate(rank = rank(Percent_resistant_predict))%>%
#     ungroup() %>%
#     mutate(y_coord = case_when(
#       rank == 1 ~ 10,
#       rank == 2 ~ 30,
#       rank == 3 ~ 50
#     )) %>%
#     dplyr::select(Host,Pathogen,Indicator,y_coord)
#   
#   
#   y_coords<-rbind(y_coords_mr,y_coords_oth) 
#   
#   
#   graph_data_complete_labels <- left_join(graph_data_complete,y_coords, by = c("Host","Pathogen","Indicator"))
#   
# # figure
#   
#   E.coli_sample_size <- Indicator_vet_data %>%
#     filter(Pathogen == "E_coli") %>%
#     group_by(Year,Host,Pathogen) %>%
#     summarise(sample_size = median(as.numeric(Sample_size), na.rm=TRUE)) %>%
#     ungroup()
#   
#   (E_coli_graph <- graph_data_complete_labels %>%
#       filter(Pathogen == "E. coli") %>%
#       mutate(Level_indicator = if_else(grepl("MDR|pan-S",Indicator), "Primary","Secondary")) %>%
#       ggplot(aes(label = label_icon))+
#       geom_bar(aes(x=Year, y = as.numeric(Percent_resistant), fill = Indicator), stat="identity", position = "dodge")+
#       geom_ribbon(aes(x=Year, fill = Indicator, ymin = CI_lower, ymax = CI_upper, y = Percent_resistant_predict), alpha = 0.25)+
#       geom_smooth(aes(x=Year, colour = Indicator, y = Percent_resistant_predict),stat = "identity")+
#       scale_fill_manual(values = BELMAP_Indicator_colourscheme,breaks=c("MDR", "pan-S", "Ciprofloxacin resistance", "Colistin resistance", "ESBL producers"))+
#       scale_colour_manual(values = BELMAP_Indicator_colourscheme,breaks=c("MDR", "pan-S", "Ciprofloxacin resistance", "Colistin resistance", "ESBL producers"))+
#       facet_grid(Level_indicator ~ Host)+
#       # geom_richtext( size = 16, hjust = 0, label.colour = NA) +
#       scale_x_continuous(limits = c(2010,2024), breaks = seq(2011,2022,1))+
#       scale_y_continuous(limits = c(-5,100), breaks = seq(0,100,50))+
#       geom_richtext(aes(x = 2023.5, y = y_coord, 
#                         label = label, 
#                         fill = Indicator), stat = "unique", 
#                     colour= "white", show.legend = FALSE)+
#       # fill = NA, label.color = NA, # remove background and outline
#       # label.padding = grid::unit(rep(0, 4), "pt")) +
#       #  geom_text(aes(x = Year, y = 2, label = Sample_size.x), na.rm = TRUE, size = 2.5)+  #make this additional for interactive report
#       labs(y= "% Resistance")&
#       moiras_graph_theme())   
#   
#   # save figure to tiff
#   figure_file_name <- paste("6_Human_AMR_Figures/",i,"_Figure.tiff",sep = "")
#   
#   




# which included in AMR?

graph_AMR <- graph_data_complete %>%
  filter(Surveillance == "AMR")

unique(graph_AMR$Pathogen)




# MODEL NSIH INCIDENCES ----------------------------------------------------------
# load processed data:

NSIH_incidence <- read_csv2( "Data/2.AMR/Human_AMR/NSIH_Incidence/Incidence.csv") %>%
  mutate(Category = paste(Pathogen,Antibiotic, sep = "_"))

#5. Run Models -------------------------------------------------

#make list of categories

combined_categories<-unique(NSIH_incidence$Category)

#calc degree variance for categories----------

var_list = list()


for(i in combined_categories){
  print(i)
  dataset_analysis<- NSIH_incidence %>%
    filter(Category == i) 
  degree_var= (max(dataset_analysis$Incidence, na.rm = TRUE)-min(dataset_analysis$Incidence, na.rm = TRUE))/mean(as.numeric(dataset_analysis$Incidence), na.rm = TRUE)
  var_list[[i]] <-degree_var
  print(degree_var)
}

#make list for overdispersion plots
overdispersion_plot_list<- list()


#make csv for graph data with header

start_data_frame <- tibble(
  "Year" = "",
  "Pathogen" = "",
  "Antibiotic"= "",
  "Number_resistant"= "",
  "Patient_days" = "",
  "Incidence"= "",
  "Incidence_predict"= "",
  "CI_upper"= "",
  "CI_lower"= "",
  "Patient_days.x" = "",
  "icon" = "",
  "signif_level" = ""
)


write_csv(start_data_frame,file = "6_Human_Incidence_GLM_predictions.csv")


# make for loop to run analysis-----
for(i in combined_categories){
  
  #3.0.subset data-----
  
  dataset_analysis_raw<- NSIH_incidence %>%
    filter(Category == i)
  min_year = min(dataset_analysis_raw$Year)
  max_year = max(dataset_analysis_raw$Year)
  dataset_analysis <- dataset_analysis_raw %>%  
    mutate(Year_simple = Year - min_year) %>%
    distinct()
  
  #make a file name to print each analysis
  
  output_file_name <- paste("6_Human_Incidence_GLM_output_summaries/GLM_output_",i,".txt",sep="") 
  
  #print title of analysis
  capture.output(print(paste("GLM analysis of ",i,sep="")),file = output_file_name)
  
  #3.1. fit models:-----
  #3.1a. fit poisson-----
  
  glmpoissonirr <- glm(Number_resistant~Year_simple + offset(log(Patient_days)), data = dataset_analysis, family = poisson(link = "log"))
  
  
  #3.1b. fit NB-----
  
  nb <- glm.nb(Number_resistant~Year_simple + offset(log(Patient_days )), data = dataset_analysis )
  
  
  #--------------------------------------
  
  # #3.2.check for overdispersion-----
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
  # dispersion_graph_file<-paste("6_Human_AMR_Dispersion_graphs/Dispersion_graph_poisson_model",i,".pdf",sep="") 
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
  length_years = max_year - min_year
  
  ndata<- tibble(
    Year = seq(min_year,max_year,by =1),
    Year_simple = seq(0,length_years,by = 1),
    Patient_days = rep(1000,length_years+1),
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
    if((max(dataset_analysis$Incidence)-min(dataset_analysis$Incidence))/mean(dataset_analysis$Incidence) > 0.25){
      model_icon <- "oscilate"
      signif = ""
    }else if((max(dataset_analysis$Incidence)-min(dataset_analysis$Incidence))/mean(dataset_analysis$Incidence) < 0.25){
      model_icon <- "equals"
      signif <- ""
    }
  }
  
  predict_data <- ndata %>%
    mutate(Incidence_predict = fit/Patient_days*1000,
           CI_upper = right_upr/Patient_days*1000,
           CI_lower = right_lwr/Patient_days*1000,
           icon = model_icon,
           signif_level = signif)
  
  
  # if non-signif - decide if equals or oscilating
  #calc if range > mean/4
  
  
  
  graph_data <- left_join(dataset_analysis,predict_data, by = c("Year","Year_simple")) %>%
    mutate(Patient_days = Patient_days.x) %>%
    dplyr::select(Year,Pathogen, Antibiotic,Number_resistant,Patient_days,Incidence,Incidence_predict,
                  CI_upper,CI_lower,Patient_days.x,icon, signif_level)
  
  
  #save dataset
  
  #   predicted_dataset_filename <- "6_Human_Incidence_GLM_predictions.csv" 
  
  #  write_csv(predict_data,file = predicted_dataset_filename)  
  write_csv(graph_data,file = "6_Human_Incidence_GLM_predictions.csv", append = TRUE)
  
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




##### Aspergillus--------------------------------------

# Due to significant changes in sample size after 2016, this dataset was reanalysed - 
#with data pre 2016 shown only as datapoints not included in trend analysis

# make for loop to run analysis-----
i = "Aspergillus_rerun"
Asperg_dataset <- read_csv2("Data/2.AMR/Human_AMR/Aspergillus/Human_AMR_Aspergillus.csv") %>%
  filter(Year>2015)

#3.0.subset data-----

dataset_analysis_raw<- Asperg_dataset %>%
  #  filter(category == i) %>%
  mutate(Year = as.numeric(Year))
min_year = min(dataset_analysis_raw$Year)
max_year = max(dataset_analysis_raw$Year)
dataset_analysis <- dataset_analysis_raw %>%  
  mutate(Year_simple = Year - min_year) %>%
  distinct()

#make a file name to print each analysis

output_file_name <- paste("6_Human_AMR_GLM_output_summaries/GLM_output_",i,".txt",sep="") 

#print title of analysis
capture.output(print(paste("GLM analysis of ",i,sep="")),file = output_file_name)

if(max_year-min_year < 3){
  capture.output(print("Too few data points for model"),file = output_file_name)
  graph_data <- dataset_analysis %>%
    mutate(Percent_resistant_predict = "",
           CI_upper = "",
           CI_lower = "",
           Sample_size.x = "Sample_size",
           icon = "",
           signif_level = "") %>%
    dplyr::select(Year,Pathogen, Antibiotic,Percent_resistant,Region, Percent_resistant_predict,Surveillance,Matrix,CI_upper,CI_lower,Sample_size.x,
                  icon, signif_level)
} else if(sum(as.numeric(dataset_analysis$Percent_resistant)) == 0){
  capture.output(print("No resistance- can't fit model"),file = output_file_name)
  graph_data <- dataset_analysis %>%
    mutate(Percent_resistant_predict = "",
           CI_upper = "",
           CI_lower = "",
           Sample_size.x = "Sample_size",
           icon = "",
           signif_level = "") %>%
    dplyr::select(Year,Pathogen, Antibiotic,Percent_resistant,Region, Percent_resistant_predict,Surveillance,Matrix,CI_upper,CI_lower,Sample_size.x,
                  icon, signif_level)
} else {
  
  #3.1. fit models:-----
  #3.1a. fit poisson-----
  
  glmpoissonirr <- glm(Number_resistants~Year_simple + offset(log(Sample_size)), data = dataset_analysis, family = poisson(link = "log"))
  
  
  #3.1b. fit NB-----
  
  nb <- glm.nb(Number_resistants~Year_simple + offset(log(Sample_size )), data = dataset_analysis )
  
  
  #--------------------------------------
  
  # #3.2.check for overdispersion-----
  # #visualise overdispersion 
  # plot(log(fitted(glmpoissonirr)),log((dataset_analysis$Number_resistants-fitted(glmpoissonirr))^2),xlab=expression(hat(mu)),ylab=expression((y-hat(mu))^2),pch=20,col="blue")
  # abline(0,1) ## 'variance = mean' line
  # 
  # overdisperions_plot<- recordPlot()
  # 
  # #assign to local environment
  # assign(paste(i,"overdispersion_plot"),print(overdisperions_plot))
  # 
  # #print dispersal graph to pdf
  # dispersion_graph_file<-paste("6_Human_AMR_Dispersion_graphs/Dispersion_graph_poisson_model",i,".pdf",sep="") 
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
    dplyr::select(Year,Pathogen, Antibiotic,Percent_resistant,Region, Percent_resistant_predict,Surveillance,Matrix,CI_upper,CI_lower,Sample_size.x,
                  icon, signif_level)
  
  
  
  
  #save dataset
  
  predicted_dataset_filename <- paste("6_Human_AMR_GLM_predictions/GLM_predictions_",i,".csv",sep="") 
  
  write_csv(predict_data,file = predicted_dataset_filename)  
  write_csv(graph_data,file = "6_Human_AMR_data_and_GLM_predictions.csv", append = TRUE)
  
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



