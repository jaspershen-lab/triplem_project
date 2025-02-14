# load data
rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")

library(tidyverse)
library(tidymass)
library(plyr)
library(microbiomedataset)
library(readxl)
###load("data)
load("3_data_analysis/gut_microbiome/data_preparation/object_cross_section")

gut_object<-object_cross_section

load("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section")

metabolomics_object<-object_cross_section
metabolite_annotation<-read_excel("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx")

setwd("1_code/gut_oral_microbiome/")

####only remain the genus level


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


microbiome_data<-gut_temp_object@expression_data
metabolite_data<-metabolomics_temp_object@expression_data
# 运行预测（自动使用并行计算）


results <- analyze_metabolite_ev(
  microbiome_data, 
  metabolite_data,
  do_feature_selection = TRUE,
  correlation_method = "spearman",
  p_threshold = 0.05,
  p_adjust_method = "none",
  rho_threshold = 0.1
)


gut_results<-results


###load("data)
load("3_data_analysis/oral_microbiome/data_preparation/object_cross_section")

oral_object<-object_cross_section


####only remain the genus level


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


microbiome_data<-oral_temp_object@expression_data
metabolite_data<-metabolomics_temp_object@expression_data

results <- analyze_metabolite_ev(
  microbiome_data, 
  metabolite_data,
  do_feature_selection = TRUE,
  correlation_method = "spearman",
  p_threshold = 0.05,
  p_adjust_method = "none",
  rho_threshold = 0.1
)


oral_results<-results



##### 合并 gut 和 oral 的GBDT数据

oral_results_summary<-oral_results$summary
gut_results_summary<-gut_results$summary

gut_oral_results_summary<-cbind(gut_results_summary[,c(1,2,6)],oral_results_summary[,c(2,6)])


colnames(gut_oral_results_summary)<-c("metabolite","gut_R2","gut_features","oral_R2","oral_features")


gut_oral_results_summary <- gut_oral_results_summary %>%
  mutate(
    gut_R2 = if_else(gut_R2 < 0.1, 0, gut_R2),
    oral_R2 = if_else(oral_R2 < 0.1, 0, oral_R2)
  )

gut_oral_results_summary <- gut_oral_results_summary %>%
  mutate(
    group = case_when(
      gut_R2 > 0 & oral_R2 == 0 ~ "gut",
      oral_R2 > 0 & gut_R2 == 0 ~ "oral",
      gut_R2 > 0 & oral_R2 > 0 ~ "co-influence",
      TRUE ~ "none"  # 这个处理其他情况，比如两个都是0的情况
    )
  )


gut_oral_results_summary

