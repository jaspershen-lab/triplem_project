  rm(list = ls())
  setwd(r4projects::get_project_wd())
  source("1_code/100_tools.R")
  source("1_code/mantel_Procrustes_code.R")
  library(tidyverse)
  library(tidymass)
  library(readxl)
  
  
  ###load("data)
  load("3_data_analysis/oral_microbiome/data_preparation/object_cross_section")
  
  oral_object<-object_cross_section
  
  load("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section")
  
  metabolomics_object<-object_cross_section
  
  metabolomics_class<-read_excel("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx")
  
  dir.create("3_data_analysis/oral_microbiome/MMC",recursive = TRUE)
  
  
  
  ####only remain the genus level
  library(microbiomedataset)
  
  oral_object <-
    oral_object %>%
    summarize_variables(what = "sum_intensity", group_by = "Genus")
  
  ##only remain the genus at least 10% samples
  dim(oral_object)
  
  non_zero_per <-
    apply(oral_object, 1, function(x) {
      sum(x != 0) / ncol(oral_object)
    })
  
  idx <-
    which(non_zero_per > 0.1)
  
  oral_object <-
    oral_object[idx, ]
  
  
  oral_object <-
    oral_object %>%
    transform2relative_intensity()
  
  
  
  ##
  ##adjust BMI, sex, and IRIS, ethnicity
  library(tidyverse)
  library(ggpubr)
  library(rstatix)
  
  oral_expression_data <-
    extract_expression_data(oral_object) %>%
    apply(1, function(x) {
      (x - mean(x)) / sd(x)
    }) %>%
    t() %>%
    as.data.frame()
  
  library(plyr)
  
  oral_sample_info <-
    oral_object@sample_info
  
  #######adjust BMI, sex, and IRIS, ethnicity
  oral_expression_data <-
    lm_adjust(expression_data = oral_expression_data,
              sample_info = oral_sample_info,
              threads = 3)
  
  oral_temp_object <- oral_object
  oral_temp_object@expression_data <- oral_expression_data
  
  
  
  # 
  expression_data <-
    extract_expression_data(metabolomics_object) %>%
    `+`(1) %>%
    log(2) %>%
    apply(1, function(x) {
      (x - mean(x)) / sd(x)
    }) %>%
    t() %>%
    as.data.frame()
  
  library(plyr)
  
  sample_info <-
    metabolomics_object@sample_info
  
  #######adjust BMI, sex, and IRIS, ethnicity
  expression_data <-
    lm_adjust(expression_data = expression_data,
              sample_info = sample_info,
              threads = 3)
  
  metabolomics_temp_object <- metabolomics_object
  metabolomics_temp_object@expression_data <- expression_data
  
  
  ########
  
  microbiome_data<-data.frame(oral_temp_object@expression_data,check.names = FALSE)
  metabolomics_data<-data.frame(metabolomics_temp_object@expression_data,check.names = FALSE)
  
  metabolomics_data<-metabolomics_data[,colnames(microbiome_data)]
  
  metabolomics_Microbial<-metabolomics_class%>%filter(HMDB.Source.Microbial=="TRUE")
  metabolomics_data<-metabolomics_data[metabolomics_Microbial$variable_id,]
  
  library(ade4)
  
  # Translated comment.
  # Translated comment.
  # Translated comment.
  
  # Translated comment.
  # Translated comment.
  microbiome_scaled <- t(microbiome_data)
  metabolome_scaled <- t(metabolomics_data)
  
  # Translated comment.
  pca_microbiome <- dudi.pca(data.frame(microbiome_scaled), scannf = FALSE, nf = 5)
  pca_metabolome <- dudi.pca(data.frame(metabolome_scaled), scannf = FALSE, nf = 5)
  
  # Translated comment.
  coia <- coinertia(pca_microbiome, pca_metabolome, scannf = FALSE, nf = 2)
  
  # Translated comment.
  # Translated comment.
  coia$RV
  
  # Translated comment.
  distances <- sqrt(rowSums((coia$mX - coia$mY)^2))
  
  # Translated comment.
  # Translated comment.
  str(coia)
  
  # Translated comment.
  coia$eig
  
  # Translated comment.
  percent_var <- (coia$eig/sum(coia$eig))*100
  print(percent_var)
  
  # Translated comment.
  # Translated comment.
  micro_scores <- coia$li  # translated comment
  micro_loadings <- coia$c1 # translated comment
  
  # Translated comment.
  metab_scores <- coia$li # translated comment
  metab_loadings <- coia$l1 # translated comment
  
  # Translated comment.
  # Translated comment.
  head(micro_loadings)
  # Translated comment.
  head(metab_loadings)
  
  
  
  
  
  
  
  a<-cbind(oral_temp_object@sample_info,distances)
  
  
  
  
  ggplot(a, aes(x=adjusted_age, y= distances)) +
    geom_point(shape=21,size=4,fill="#A1D0C7",color="white") +
    geom_smooth(method="lm",colour = "grey50") +theme_light() +stat_cor(method = "pearson")+theme(legend.position="none", # translated comment
                                                                                                  axis.text.x=element_text(colour="black",size=14), # translated comment
                                                                                                  axis.text.y=element_text(size=14,face="plain"), # translated comment
                                                                                                  axis.title.y=element_text(size = 14,face="plain"), # translated comment
                                                                                                  axis.title.x=element_text(size = 14,face="plain"), # translated comment
                                                                                                  plot.title = element_text(size=15,face="bold",hjust = 0.5))+xlab("Age")