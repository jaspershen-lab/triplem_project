# load data
rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")
source("1_code/mantel_Procrustes_code.R")
library(tidyverse)
library(tidymass)
library(readxl)
library(ade4)
###load("data)
load("3_data_analysis/gut_microbiome/data_preparation/object_cross_section")

gut_object<-object_cross_section

load("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section")

metabolomics_object<-object_cross_section
metabolomics_class<-read_excel("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx")
dir.create("3_data_analysis/gut_microbiome/spearman/cross_section/",recursive = TRUE)

setwd("3_data_analysis/gut_microbiome/spearman/cross_section/")




####only remain the genus level
library(microbiomedataset)

gut_object <-
  gut_object %>%
  summarize_variables(what = "sum_intensity", group_by = "Genus")

##only remain the genus at least 10% samples
dim(gut_object)

non_zero_per <-
  apply(gut_object, 1, function(x) {
    sum(x != 0) / ncol(gut_object)
  })

idx <-
  which(non_zero_per > 0.1)

gut_object <-
  gut_object[idx, ]


gut_object <-
  gut_object %>%
  transform2relative_intensity()



##
##adjust BMI, sex, and IRIS, ethnicity
library(tidyverse)
library(ggpubr)
library(rstatix)

gut_expression_data <-
  extract_expression_data(gut_object) %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

library(plyr)

gut_sample_info <-
  gut_object@sample_info

#######adjust BMI, sex, and IRIS, ethnicity
gut_expression_data <-
  lm_adjust(expression_data = gut_expression_data,
            sample_info = gut_sample_info,
            threads = 3)

gut_temp_object <- gut_object
gut_temp_object@expression_data <- gut_expression_data



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







gut_data <- gut_temp_object@expression_data 
metabolome_data <- metabolomics_temp_object@expression_data

metabolomics_class<-metabolomics_class%>%filter(HMDB.Source.Microbial=="TRUE")
metabolome_data<-metabolome_data[metabolomics_class$variable_id,]
sample_info <- metabolomics_temp_object@sample_info
rownames(sample_info)<-sample_info$sample_id
# Translated comment.
selection_results <- select_significant_metabolites(
  gut_microbiome = gut_data,
  metabolome = metabolome_data,
  cor_threshold = 0.3,
  p_threshold = 0.05
)

# Translated comment.
metabolite_plot <- plot_metabolite_selection(selection_results,metabolomics_class = metabolomics_class)

gut_data<-gut_data

# Translated comment.
results <- analyze_associations(
  gut_microbiome = gut_data[,], 
  metabolome = metabolome_data,
  selected_metabolites = selection_results$significant_metabolites,
  metadata =  sample_info  # translated comment
)

# Translated comment.
plots <- plot_associations(results)

distance<-plots$distances$data
distance<-cbind(distance,gut_temp_object@sample_info[,15:22])
distance<-subset(distance,diabetes_class=="Prediabetic")
ggplot(distance, aes(x=distance$adjusted_age, y= coinertia_distance)) +
  geom_point(shape=21,size=4,fill="#A1D0C7",color="white") +
  geom_smooth(method="lm",colour = "grey50") +theme_light() +stat_cor(method = "pearson")+theme(legend.position="none", # translated comment
                                                                                                axis.text.x=element_text(colour="black",size=14), # translated comment
                                                                                                axis.text.y=element_text(size=14,face="plain"), # translated comment
                                                                                                axis.title.y=element_text(size = 14,face="plain"), # translated comment
                                                                                                axis.title.x=element_text(size = 14,face="plain"), # translated comment
                                                                                                plot.title = element_text(size=15,face="bold",hjust = 0.5))+xlim(c(33,70))