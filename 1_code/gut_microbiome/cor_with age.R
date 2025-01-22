rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")
source("1_code/mantel_Procrustes_code.R")
library(tidyverse)
library(tidymass)
library(readxl)


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





## 

gut_microbiome<-data.frame(t(gut_expression_data),check.names = FALSE)

colnames(gut_microbiome)<-gut_temp_object@variable_info$Genus
gut_microbiome$sample_id<-rownames(gut_microbiome)


plasma_metabolome<-data.frame(t(expression_data),check.names = FALSE)
plasma_metabolome$sample_id<-rownames(plasma_metabolome)



sample_info_mm<-sample_info%>%
  full_join(gut_microbiome, by = "sample_id") %>%
  full_join(plasma_metabolome, by = "sample_id") 


sample_info_mm<-subset(sample_info_mm,!IRIS=="IS")

ggplot(sample_info_mm, aes(x=adjusted_age, y= M263T353_NEG_HILIC)) +
  geom_point(shape=21,size=4,fill="#A1D0C7",color="white") +
  geom_smooth(method="lm",colour = "grey50") +theme_light() +stat_cor(method = "pearson")+theme(legend.position="none", #不需要图例
                                                                                                axis.text.x=element_text(colour="black",size=14), #设置x轴刻度标签的字体属性
                                                                                                axis.text.y=element_text(size=14,face="plain"), #设置x轴刻度标签的字体属性
                                                                                                axis.title.y=element_text(size = 14,face="plain"), #设置y轴的标题的字体属性
                                                                                                axis.title.x=element_text(size = 14,face="plain"), #设置x轴的标题的字体属性
                                                                                                plot.title = element_text(size=15,face="bold",hjust = 0.5))+
  labs(title = "Phenylacetylglutamine",
       y = "Standar Adundance",
       x = "Age")+ylim(c(-2,2.5))


