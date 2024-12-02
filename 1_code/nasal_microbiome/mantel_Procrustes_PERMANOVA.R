# load data
rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")
source("1_code/mantel_Procrustes_code.R")
library(tidyverse)
library(tidymass)

###load("data)
load("3_data_analysis/nasal_microbiome/data_preparation/object_cross_section")

gut_object<-object_cross_section

load("3_data_analysis/nasal_metabolomics/data_preparation/metabolite/object_cross_section")

metabolomics_object<-object_cross_section

dir.create("3_data_analysis/nasal_microbiome/spearman/cross_section/",recursive = TRUE)

setwd("3_data_analysis/nasal_microbiome/spearman/cross_section/")




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
sample_info <- metabolomics_temp_object@sample_info
rownames(sample_info)<-sample_info$sample_id
# 1. 筛选显著相关的代谢物
selection_results <- select_significant_metabolites(
  gut_microbiome = gut_data,
  metabolome = metabolome_data,
  cor_threshold = 0.3,
  p_threshold = 0.05
)

# 2. 可视化代谢物筛选结果
metabolite_plot <- plot_metabolite_selection(selection_results)

# 3. 分析关联
results <- analyze_associations(
  gut_microbiome = gut_data, 
  metabolome = metabolome_data,
  selected_metabolites = selection_results$significant_metabolites,
  metadata =  sample_info  # 包含IRIS列的元数据
)

# 4. 绘制结果
plots <- plot_associations(results)