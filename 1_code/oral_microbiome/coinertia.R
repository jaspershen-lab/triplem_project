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
  
  # 假设你的数据是这样的:
  # microbiome_data: 行是样本,列是菌群
  # metabolome_data: 行是相同的样本,列是代谢物
  
  # 1. 数据预处理
  # 标准化数据,确保行名一致
  microbiome_scaled <- t(microbiome_data)
  metabolome_scaled <- t(metabolomics_data)
  
  # 2. 分别对两个数据集进行PCA
  pca_microbiome <- dudi.pca(data.frame(microbiome_scaled), scannf = FALSE, nf = 5)
  pca_metabolome <- dudi.pca(data.frame(metabolome_scaled), scannf = FALSE, nf = 5)
  
  # 3. 进行Coinertia分析
  coia <- coinertia(pca_microbiome, pca_metabolome, scannf = FALSE, nf = 2)
  
  # 4. 查看结果
  # RV系数(整体相关性)
  coia$RV
  
  # 计算每个样本的投影距离
  distances <- sqrt(rowSums((coia$mX - coia$mY)^2))
  
  # 为结果创建一个数据框
  # 查看coinertia对象的结构
  str(coia)
  
  # 查看特征值（eigenvalues）
  coia$eig
  
  # 查看各轴的解释百分比
  percent_var <- (coia$eig/sum(coia$eig))*100
  print(percent_var)
  
  # 获取两组数据在共同空间的坐标
  # 菌群数据的坐标
  micro_scores <- coia$li  # 样本在菌群空间的坐标
  micro_loadings <- coia$c1 # 菌群变量的贡献
  
  # 代谢物数据的坐标
  metab_scores <- coia$li # 样本在代谢物空间的坐标
  metab_loadings <- coia$l1 # 代谢物变量的贡献
  
  # 查看变量贡献
  # 菌群变量贡献
  head(micro_loadings)
  # 代谢物变量贡献
  head(metab_loadings)
  
  
  
  
  
  
  
  a<-cbind(oral_temp_object@sample_info,distances)
  
  
  
  
  ggplot(a, aes(x=adjusted_age, y= distances)) +
    geom_point(shape=21,size=4,fill="#A1D0C7",color="white") +
    geom_smooth(method="lm",colour = "grey50") +theme_light() +stat_cor(method = "pearson")+theme(legend.position="none", #不需要图例
                                                                                                  axis.text.x=element_text(colour="black",size=14), #设置x轴刻度标签的字体属性
                                                                                                  axis.text.y=element_text(size=14,face="plain"), #设置x轴刻度标签的字体属性
                                                                                                  axis.title.y=element_text(size = 14,face="plain"), #设置y轴的标题的字体属性
                                                                                                  axis.title.x=element_text(size = 14,face="plain"), #设置x轴的标题的字体属性
                                                                                                  plot.title = element_text(size=15,face="bold",hjust = 0.5))+xlab("Age")