rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")

library(tidyverse)
library(tidymass)
library(readxl)
library(compositions)

## load data
### for gut microbiome

load(file = "2_data/gut_microbiome/data_preparation/expression_data")
load(file = "2_data/gut_microbiome/data_preparation/sample_info")
load(file = "2_data/gut_microbiome/data_preparation/variable_info")


gut_microbiome_table<-expression_data
gut_microbiome_metadata<-sample_info
gut_microbiome_tax<-variable_info

# 标准化微生物组数据

gut_microbiome_metadata <- gut_microbiome_metadata %>% filter(!is.na(Gender))
gut_microbiome_metadata <- gut_microbiome_metadata %>% filter(!is.na(IRIS))
gut_microbiome_metadata <- gut_microbiome_metadata %>% filter(!is.na(BMI))


gut_microbiome_table<-gut_microbiome_table[,gut_microbiome_metadata$sample_id]


gut_microbiome_table<-microbiome_genus_filt(gut_microbiome_table,gut_microbiome_tax,0.1)

gut_microbiome_table_clr <- t(data.frame(clr(gut_microbiome_table)))


gut_microbiome_table <-lm_adjust(expression_data = gut_microbiome_table_clr,
                                 sample_info = gut_microbiome_metadata,
                                 threads = 3)

gut_microbiome_table<-data.frame(t(gut_microbiome_table))


### for oral microbiome

load(file = "2_data/oral_microbiome/data_preparation/expression_data")
load(file = "2_data/oral_microbiome/data_preparation/sample_info")
load(file = "2_data/oral_microbiome/data_preparation/variable_info")


oral_microbiome_table<-expression_data
oral_microbiome_metadata<-sample_info
oral_microbiome_tax<-variable_info

# 标准化微生物组数据

oral_microbiome_metadata <- oral_microbiome_metadata %>% filter(!is.na(Gender))
oral_microbiome_metadata <- oral_microbiome_metadata %>% filter(!is.na(IRIS))
oral_microbiome_metadata <- oral_microbiome_metadata %>% filter(!is.na(BMI))


oral_microbiome_table<-oral_microbiome_table[,oral_microbiome_metadata$sample_id]


oral_microbiome_table<-microbiome_genus_filt(oral_microbiome_table,oral_microbiome_tax,0.1)

oral_microbiome_table_clr <- t(data.frame(clr(oral_microbiome_table)))


oral_microbiome_table <-lm_adjust(expression_data = oral_microbiome_table_clr,
                                  sample_info = oral_microbiome_metadata,
                                  threads = 3)

oral_microbiome_table<-data.frame(t(oral_microbiome_table))


### for skin microbiome

load(file = "2_data/skin_microbiome/data_preparation/expression_data")
load(file = "2_data/skin_microbiome/data_preparation/sample_info")
load(file = "2_data/skin_microbiome/data_preparation/variable_info")


skin_microbiome_table<-expression_data
skin_microbiome_metadata<-sample_info
skin_microbiome_tax<-variable_info

# 标准化微生物组数据

skin_microbiome_metadata <- skin_microbiome_metadata %>% filter(!is.na(Gender))
skin_microbiome_metadata <- skin_microbiome_metadata %>% filter(!is.na(IRIS))
skin_microbiome_metadata <- skin_microbiome_metadata %>% filter(!is.na(BMI))


skin_microbiome_table<-skin_microbiome_table[,skin_microbiome_metadata$sample_id]


skin_microbiome_table<-microbiome_genus_filt(skin_microbiome_table,skin_microbiome_tax,0.1)

skin_microbiome_table_clr <- t(data.frame(clr(skin_microbiome_table)))


skin_microbiome_table <-lm_adjust(expression_data = skin_microbiome_table_clr,
                                  sample_info = skin_microbiome_metadata,
                                  threads = 3)

skin_microbiome_table<-data.frame(t(skin_microbiome_table))


### for nasal microbiome

load(file = "2_data/nasal_microbiome/data_preparation/expression_data")
load(file = "2_data/nasal_microbiome/data_preparation/sample_info")
load(file = "2_data/nasal_microbiome/data_preparation/variable_info")


nasal_microbiome_table<-expression_data
nasal_microbiome_metadata<-sample_info
nasal_microbiome_tax<-variable_info


# 标准化微生物组数据

nasal_microbiome_metadata <- nasal_microbiome_metadata %>% filter(!is.na(Gender))
nasal_microbiome_metadata <- nasal_microbiome_metadata %>% filter(!is.na(IRIS))
nasal_microbiome_metadata <- nasal_microbiome_metadata %>% filter(!is.na(BMI))


nasal_microbiome_table<-nasal_microbiome_table[,nasal_microbiome_metadata$sample_id]


nasal_microbiome_table<-microbiome_genus_filt(nasal_microbiome_table,nasal_microbiome_tax,0.1)

nasal_microbiome_table_clr <- t(data.frame(clr(nasal_microbiome_table)))


nasal_microbiome_table <-lm_adjust(expression_data = nasal_microbiome_table_clr,
                                   sample_info = nasal_microbiome_metadata,
                                   threads = 3)

nasal_microbiome_table<-data.frame(t(nasal_microbiome_table))


## load data
## load data
library(readxl)

expression_data_metabolome <- read_excel("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/expression_data.xlsx")

sample_info_metabolome <- data.frame(read_excel("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/sample_info.xlsx"))

variable_info_metabolome <- data.frame(read_excel("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info.xlsx"))


row.names(expression_data_metabolome)<-variable_info_metabolome$variable_id
# data transform and 

sample_info_metabolome <- sample_info_metabolome %>% filter(!is.na(Gender))
sample_info_metabolome <- sample_info_metabolome %>% filter(!is.na(IRIS))
sample_info_metabolome <- sample_info_metabolome %>% filter(!is.na(BMI))


expression_data_metabolome<-expression_data_metabolome[,sample_info_metabolome$sample_id]


expression_data_metabolome<-log2(expression_data_metabolome)

expression_data_metabolome <-lm_adjust(expression_data = expression_data_metabolome,
                                       sample_info = sample_info_metabolome,
                                       threads = 3)

expression_data_metabolome<-data.frame(t(expression_data_metabolome))


colnames(expression_data_metabolome)<-variable_info_metabolome$variable_id


# 合并四个部位数据

Common_sample<-Reduce(intersect,list(row.names(expression_data_metabolome),row.names(gut_microbiome_table),row.names(oral_microbiome_table),row.names(nasal_microbiome_table),row.names(skin_microbiome_table)))

expression_data_metabolome<-expression_data_metabolome[Common_sample,]

gut_microbiome_table<-gut_microbiome_table[Common_sample,]
oral_microbiome_table<-oral_microbiome_table[Common_sample,]
nasal_microbiome_table<-nasal_microbiome_table[Common_sample,]
skin_microbiome_table<-skin_microbiome_table[Common_sample,]


colnames(gut_microbiome_table)<-paste0(colnames(gut_microbiome_table),"_gut")
colnames(oral_microbiome_table)<-paste0(colnames(oral_microbiome_table),"_oral")
colnames(nasal_microbiome_table)<-paste0(colnames(nasal_microbiome_table),"_nasal")
colnames(skin_microbiome_table)<-paste0(colnames(skin_microbiome_table),"_skin")


all_microbiome_data<-cbind(gut_microbiome_table,oral_microbiome_table,nasal_microbiome_table,skin_microbiome_table)


MM_adonis <- function(gut_microbiome_table, expression_data_metabolome) {
  library(compositions)
  library(vegan)
  # 与表达数据的交集
  
  shared_samples <- intersect(row.names(gut_microbiome_table), row.names(expression_data_metabolome))
  expression_data_metabolome <- expression_data_metabolome[shared_samples, ]
  gut_microbiome_table<-gut_microbiome_table[shared_samples,]
  
  # 计算距离并执行adonis分析
  metabolome_dist <- vegdist(expression_data_metabolome, method = "euclidean")
  results <- NULL
  for (i in colnames(gut_microbiome_table)) {
    adonis_result <- adonis2(metabolome_dist ~ gut_microbiome_table[, i], permutations = 1000,parallel=15)
    results <- rbind(results, adonis_result$`Pr(>F)`[1])
  }
  results<-data.frame(results)
  row.names(results) <- colnames(gut_microbiome_table)
  colnames(results)<-c("Pr(>F)")
  significant_results <- subset(results, `Pr(>F)` < 0.01)
  
  adonis_micro_sig <- gut_microbiome_table[, row.names(significant_results)]
  final_result <- adonis2(metabolome_dist ~ ., data = adonis_micro_sig, permutations = 1000,parallel=15)
  
  return(list(final_result = final_result, significant_results = significant_results))
}


all_microbiome_adonis<-MM_adonis(all_microbiome_data,expression_data_metabolome)


metabolite_annotation<-read_excel("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx")

metabolite_annotation<-subset(metabolite_annotation,metabolite_annotation$HMDB.Source.Microbial=="TRUE")
expression_data_metabolome_microbial<-expression_data_metabolome[,metabolite_annotation$variable_id]

# 计算微生物来源的代谢物的解释度

all_microbiome_data_micro<-MM_adonis(all_microbiome_data,expression_data_metabolome_microbial)