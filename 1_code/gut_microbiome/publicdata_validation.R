rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")


## read MetaCardis data
library(readxl)

MetaCardis_metagenomic<-data.frame(read_excel("2_data/MetaCardis cohort/European MetaCardis cohort _metagenomic.xlsx"))
row.names(MetaCardis_metagenomic)<-MetaCardis_metagenomic$ID

MetaCardis_metagenomic<-MetaCardis_metagenomic[,-1]


MetaCardis_tax<-data.frame(read_excel("2_data/MetaCardis cohort/taxonomy.xlsx"))

MetaCardis_metabolome<-data.frame(read_excel("2_data/MetaCardis cohort/European MetaCardis cohort_metabolome.xlsx"))

row.names(MetaCardis_metabolome)<-MetaCardis_metabolome$ID

MetaCardis_metabolome<-MetaCardis_metabolome[,-1]

MetaCardis_metabome_annotation<-data.frame(read_excel("2_data/MetaCardis cohort/metabolites_annotation.xlsx"))

MetaCardis_metadata<-data.frame(read_excel("2_data/MetaCardis cohort/metadata.xlsx"))


# 去除重复的人群

MetaCardis_metadata<-subset(MetaCardis_metadata,Status%in%c("MMC372","HC275","UMCC222","IHD372"))

MetaCardis_metagenomic<-MetaCardis_metagenomic[MetaCardis_metadata$ID,]

MetaCardis_metabolome<-MetaCardis_metabolome[MetaCardis_metadata$ID,]


# 筛选同时拥有代谢组和微生物组的人群
MetaCardis_metagenomic <- MetaCardis_metagenomic[!apply(MetaCardis_metagenomic == "NA", 1, all), ]
MetaCardis_metabolome <- MetaCardis_metabolome[!apply(MetaCardis_metabolome == "NA", 1, all), ]

com_sample<-intersect(row.names(MetaCardis_metagenomic),row.names(MetaCardis_metabolome))

MetaCardis_metagenomic<-MetaCardis_metagenomic[com_sample,]
MetaCardis_metabolome<-MetaCardis_metabolome[com_sample,]
MetaCardis_metadata<-subset(MetaCardis_metadata,ID%in%com_sample)



# 将MetaCardis_metagenomic合并为Genus levels




MetaCardis_metagenomic <- type.convert(MetaCardis_metagenomic, as.is = TRUE)

MetaCardis_metagenomic<-data.frame(t(MetaCardis_metagenomic))

MetaCardis_metagenomic<-cbind(MetaCardis_metagenomic,MetaCardis_tax[,8:4])


MetaCardis_metagenomic<-aggregate(MetaCardis_metagenomic[,1:1146],list(MetaCardis_metagenomic$genus),sum)

rownames(MetaCardis_metagenomic)<-  MetaCardis_metagenomic$Group.1

MetaCardis_metagenomic<-MetaCardis_metagenomic[,-1]

MetaCardis_metagenomic<-t(MetaCardis_metagenomic)
MetaCardis_metagenomic<-MetaCardis_metagenomic/apply(MetaCardis_metagenomic,1,sum)*100



# 计算相对丰度

# 读取iPOP数据


library(tidyverse)
library(tidymass)
library(plyr)
library(microbiomedataset)
###load("data)
load("3_data_analysis/gut_microbiome/data_preparation/object_cross_section")

gut_object<-object_cross_section

load("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section")

metabolomics_object<-object_cross_section

metabolite_annotation<-read_excel("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx")




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
  which(non_zero_per > 0.2)

gut_object <-
  gut_object[idx, ]


gut_object <-
  gut_object %>%
  transform2relative_intensity()

gut_expression_data <-
  extract_expression_data(gut_object) %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()



gut_sample_info <-
  gut_object@sample_info


gut_temp_object <- gut_object
gut_temp_object@expression_data <- gut_expression_data


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



##筛选共同的属
com_tax<-intersect(gut_temp_object@variable_info$Genus,colnames(MetaCardis_metagenomic))


microbiome_data_ipop<-gut_temp_object@expression_data
microbiome_tax_ipop<-gut_temp_object@variable_info

rownames(microbiome_data_ipop)<-microbiome_tax_ipop$Genus

microbiome_data_ipop<-microbiome_data_ipop[com_tax,]

MetaCardis_metagenomic<-data.frame(t(MetaCardis_metagenomic))

MetaCardis_metagenomic<-MetaCardis_metagenomic[com_tax,]

MetaCardis_metagenomic <-
  MetaCardis_metagenomic%>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()


### 挑选共同的代谢物

metabolome_data_ipop<-metabolomics_temp_object@expression_data
metabolome_ann_ipop<-data.frame(metabolite_annotation)

com_metabolites<-intersect(metabolome_ann_ipop$HMDB,MetaCardis_metabome_annotation$HMDB)
#去除NA
com_metabolites<-com_metabolites[-13]


# 

metabolome_ann_ipop<-subset(metabolome_ann_ipop,HMDB%in%com_metabolites)

metabolome_data_ipop<-metabolome_data_ipop[metabolome_ann_ipop$variable_id,]


MetaCardis_metabome_annotation<-subset(MetaCardis_metabome_annotation,HMDB%in%com_metabolites)

MetaCardis_metabolome<-MetaCardis_metabolome[,MetaCardis_metabome_annotation$Met_ID]

MetaCardis_metabolome <- type.convert(MetaCardis_metabolome, as.is = TRUE)



### 
microbiome_data_ipop<-data.frame(t(microbiome_data_ipop))

metabolome_data_ipop<-data.frame(t(metabolome_data_ipop))


results <- analyze_metabolite_ev(
  data.frame(t(microbiome_data_ipop)), 
  data.frame(t(metabolome_data_ipop)),
  do_feature_selection = TRUE,
  correlation_method = "spearman",
  p_threshold = 0.05,
  p_adjust_method = "none",
  rho_threshold = 0.1
)






new_predictions <- predict_metabolite(
  "metabolite_models/M182T470_2_POS_HILIC_model.rds", 
  data.frame(t(MetaCardis_metagenomic))
)


cor.test(new_predictions$predictions,MetaCardis_metabolome$Met00172,method = "pearson")



library(tidyverse)

predict_metabolites <- function(models_directory, new_data) {
  # 获取目录下所有RDS模型文件
  model_files <- list.files(models_directory, pattern = "\\.rds$", full.names = TRUE)
  
  # 初始化结果列表
  predictions_list <- list()
  
  # 遍历每个模型文件
  for (model_file in model_files) {
    # 加载模型信息
    model_info <- readRDS(model_file)
    
    # 确保输入数据包含所需的特征
    required_features <- model_info$selected_features
    missing_features <- setdiff(required_features, colnames(new_data))
    
    if (length(missing_features) > 0) {
      warning(paste("Skipping model", basename(model_file), "- missing required features:", 
                    paste(missing_features, collapse = ", ")))
      next
    }
    
    # 选择需要的特征并确保顺序正确
    filtered_data <- new_data[, required_features, drop = FALSE]
    
    # 进行预测
    predictions <- predict(model_info$model, filtered_data, n.trees = model_info$parameters$n.trees)
    
    # 保存预测结果
    predictions_list[[basename(model_file)]] <- list(
      predictions = predictions,
      model_performance = model_info$performance,
      feature_importance = model_info$feature_importance
    )
  }
  
  return(predictions_list)
}

# 调用函数并提供模型目录和数据
models_directory <- "metabolite_models"
new_data <- data.frame(t(MetaCardis_metagenomic))  # 确保输入数据格式正确

results <- predict_metabolites(models_directory, new_data)


predictions_df <- do.call(cbind, lapply(results, function(x) x$predictions))

# 设置列名为模型名称
colnames(predictions_df) <- names(results)

# 转换为数据框格式（如果需要）
predictions_df <- as.data.frame(predictions_df)

colnames(predictions_df)<-gsub("_model.rds","",colnames(predictions_df))


predictions_df
########

R2_df <- do.call(cbind, lapply(results, function(x) x$model_performance$mean_r2))

# 设置列名为模型名称
colnames(R2_df) <- names(results)


#######

predictions_df<-data.frame(t(predictions_df))

colnames(predictions_df)<-colnames(MetaCardis_metagenomic)
predictions_df$variable_id<-rownames(predictions_df)
predictions_df<-merge(predictions_df,metabolite_annotation,by="variable_id")

metabolites_name<-unique(predictions_df$HMDB)

for (i in metabolites_name){
  
  predictions_data<-subset(predictions_df,HMDB==i)[,2:1147]
  predictions_data<-colSums(predictions_data)
  
  obersered_sample_name<-subset(MetaCardis_metabome_annotation,HMDB==i)
  obersered_sample_name<-obersered_sample_name$Met_ID
  for (z in obersered_sample_name){
    
    obersered_data<-MetaCardis_metabolome[,z]
    
    merge_data<-rbind(predictions_data,obersered_data)
    merge_data<-data.frame(t(merge_data))
    colnames(merge_data)<-c("predictions","obersered")
    
    merge_data<-na.omit(merge_data)
    
    results_cor<-cor.test(merge_data$predictions,merge_data$obersered)
    
    results_cor<-data.frame(results_cor$p.value,results_cor$estimate)
    
    results_cor$HMDB_Name<-i
    
  }
  
  results_cor_name<-rbind(results_cor_name,results_cor)
  
}


