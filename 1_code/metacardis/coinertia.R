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
metabolite_annotation<-read_excel("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx")
#metabolite_annotation<-subset(metabolite_annotation,metabolite_annotation$HMDB.Source.Microbial=="TRUE")

MetaCardis_metabome_annotation_microbial<-subset(MetaCardis_metabome_annotation,HMDB%in%metabolite_annotation$HMDB)

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
MetaCardis_metabolome <- type.convert(MetaCardis_metabolome, as.is = TRUE)

MetaCardis_metabolome<-na.omit(MetaCardis_metabolome)
# 将MetaCardis_metagenomic合并为Genus levels




MetaCardis_metagenomic <- type.convert(MetaCardis_metagenomic, as.is = TRUE)

MetaCardis_metagenomic<-data.frame(t(MetaCardis_metagenomic))

MetaCardis_metagenomic<-cbind(MetaCardis_metagenomic,MetaCardis_tax[,8:4])


MetaCardis_metagenomic<-aggregate(MetaCardis_metagenomic[,1:1146],list(MetaCardis_metagenomic$genus),sum)

rownames(MetaCardis_metagenomic)<-  MetaCardis_metagenomic$Group.1

MetaCardis_metagenomic<-MetaCardis_metagenomic[,-1]

MetaCardis_metagenomic<-t(MetaCardis_metagenomic)
MetaCardis_metagenomic<-MetaCardis_metagenomic/apply(MetaCardis_metagenomic,1,sum)*100

MetaCardis_metagenomic<-MetaCardis_metagenomic[rownames(MetaCardis_metabolome),]
MetaCardis_metadata<-MetaCardis_metadata%>%filter(ID%in%rownames(MetaCardis_metabolome))



library(ade4)




# 读取数据
metabolome <-MetaCardis_metagenomic
metadata <- MetaCardis_metadata

# 确保代谢物数据是数值型
metabolome <- as.data.frame(apply(metabolome, 2, as.numeric))

# 提取年龄数据
age <- metadata$Age..years.

# 创建一个空的数据框来存储相关性结果
correlation_results <- data.frame(
  Metabolite = colnames(metabolome),
  Correlation = NA,
  P_value = NA
)

# 计算每个代谢物与年龄的相关性
for(i in 1:ncol(metabolome)) {
  # 使用complete.cases确保只使用完整的数据对
  complete_data <- complete.cases(metabolome[,i], age)
  
  # 计算相关性
  cor_test <- cor.test(metabolome[complete_data,i], 
                       age[complete_data], 
                       method = "pearson",  # 使用Spearman相关系数
                       exact = FALSE)        # 对于大样本量使用近似计算
  
  # 存储结果
  correlation_results$Correlation[i] <- cor_test$estimate
  correlation_results$P_value[i] <- cor_test$p.value
}

# 添加FDR校正的P值
correlation_results$FDR <- p.adjust(correlation_results$P_value, method = "BH")

# 按相关系数绝对值排序
correlation_results <- correlation_results[order(abs(correlation_results$Correlation), decreasing = TRUE),]

# 添加显著性标记
correlation_results$Significance <- ifelse(correlation_results$FDR < 0.05, 
                                           "Significant", 
                                           "Not Significant")


correlation_results_MGS<-subset(correlation_results,Correlation>0&Significance=="Significant")



metabolome <-MetaCardis_metabolome
metadata <- MetaCardis_metadata

# 确保代谢物数据是数值型
metabolome <- as.data.frame(apply(metabolome, 2, as.numeric))

# 提取年龄数据
age <- metadata$Age..years.

# 创建一个空的数据框来存储相关性结果
correlation_results <- data.frame(
  Metabolite = colnames(metabolome),
  Correlation = NA,
  P_value = NA
)

# 计算每个代谢物与年龄的相关性
for(i in 1:ncol(metabolome)) {
  # 使用complete.cases确保只使用完整的数据对
  complete_data <- complete.cases(metabolome[,i], age)
  
  # 计算相关性
  cor_test <- cor.test(metabolome[complete_data,i], 
                       age[complete_data], 
                       method = "pearson",  # 使用Spearman相关系数
                       exact = FALSE)        # 对于大样本量使用近似计算
  
  # 存储结果
  correlation_results$Correlation[i] <- cor_test$estimate
  correlation_results$P_value[i] <- cor_test$p.value
}

# 添加FDR校正的P值
correlation_results$FDR <- p.adjust(correlation_results$P_value, method = "BH")

# 按相关系数绝对值排序
correlation_results <- correlation_results[order(abs(correlation_results$Correlation), decreasing = TRUE),]

# 添加显著性标记
correlation_results$Significance <- ifelse(correlation_results$FDR < 0.05, 
                                           "Significant", 
                                           "Not Significant")

correlation_results_Met<-subset(correlation_results,Correlation>0.1&Significance=="Significant")




# 2. 分别对两个数据集进行PCA
pca_microbiome <- dudi.pca(data.frame(MetaCardis_metagenomic), scannf = FALSE, nf = 5)
pca_metabolome <- dudi.pca(data.frame(MetaCardis_metabolome[,correlation_results_Met$Metabolite]), scannf = FALSE, nf = 5)

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





a<-cbind(MetaCardis_metadata,distances)
a<- type.convert(a, as.is = TRUE)

a<-subset(a,a$Status=="UMCC222")
ggplot(a, aes(x=Age..years., y= distances)) +
  geom_point(shape=21,size=4,fill="#A1D0C7",color="white") +
  geom_smooth(method="lm",colour = "grey50") +theme_light() +stat_cor(method = "pearson")+theme(legend.position="none", #不需要图例
                                                                                                axis.text.x=element_text(colour="black",size=14), #设置x轴刻度标签的字体属性
                                                                                                axis.text.y=element_text(size=14,face="plain"), #设置x轴刻度标签的字体属性
                                                                                                axis.title.y=element_text(size = 14,face="plain"), #设置y轴的标题的字体属性
                                                                                                axis.title.x=element_text(size = 14,face="plain"), #设置x轴的标题的字体属性
                                                                                                plot.title = element_text(size=15,face="bold",hjust = 0.5))+xlab("Age")+xlim(c(27,70))



a<-subset(a,Age..years.>27&Age..years.<70)

