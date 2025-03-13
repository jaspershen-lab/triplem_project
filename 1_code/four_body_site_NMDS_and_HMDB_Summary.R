rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")
setwd("1_code/4_site_merge/")
library(tidymass)
library(tidyverse)
library(readxl)
###load("data)
load("../../3_data_analysis/gut_microbiome/data_preparation/object_cross_section")

gut_object<-object_cross_section

load("../../3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section")

metabolomics_object<-object_cross_section


load("../../3_data_analysis/oral_microbiome/data_preparation/object_cross_section")

oral_object<-object_cross_section

load("../../3_data_analysis/skin_microbiome/data_preparation/object_cross_section")

skin_object<-object_cross_section

load("../../3_data_analysis/nasal_microbiome/data_preparation/object_cross_section")
nasal_object<-object_cross_section
metabolite_annotation<-read_excel("../../3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx")
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







####only remain the genus level
library(microbiomedataset)

skin_object <-
  skin_object %>%
  summarize_variables(what = "sum_intensity", group_by = "Genus")

##only remain the genus at least 10% samples
dim(skin_object)

non_zero_per <-
  apply(skin_object, 1, function(x) {
    sum(x != 0) / ncol(skin_object)
  })

idx <-
  which(non_zero_per > 0.1)

skin_object <-
  skin_object[idx, ]


skin_object <-
  skin_object %>%
  transform2relative_intensity()





####only remain the genus level
library(microbiomedataset)

nasal_object <-
  nasal_object %>%
  summarize_variables(what = "sum_intensity", group_by = "Genus")

##only remain the genus at least 10% samples
dim(nasal_object)

non_zero_per <-
  apply(nasal_object, 1, function(x) {
    sum(x != 0) / ncol(nasal_object)
  })

idx <-
  which(non_zero_per > 0.1)

nasal_object <-
  nasal_object[idx, ]


nasal_object <-
  nasal_object %>%
  transform2relative_intensity()










# 合并四个部位的微生物组数据绘制PCOA图

# 读取四个不同区域的微生物组数据
# 假设文件路径为当前工作目录
gut_genus<-gut_object@expression_data
rownames(gut_genus)<-gut_object@variable_info$Genus

oral_genus<-oral_object@expression_data
rownames(oral_genus)<-oral_object@variable_info$Genus

skin_genus<-skin_object@expression_data
rownames(skin_genus)<-skin_object@variable_info$Genus

nasal_genus<-nasal_object@expression_data
rownames(nasal_genus)<-nasal_object@variable_info$Genus

# 加载必要的R包
library(vegan)      # 用于计算生态距离和NMDS
library(ggplot2)    # 用于绘图
library(readr)      # 用于读取CSV文件
library(dplyr)      # 用于数据处理
library(tidyr)      # 用于数据整理


# 为每个数据集添加来源标签
gut_samples <- colnames(gut_genus)
oral_samples <- colnames(oral_genus)
skin_samples <- colnames(skin_genus)
nasal_samples <- colnames(nasal_genus)

# 整合所有数据
# 转置矩阵使行为样本，列为物种
gut_t <- t(gut_genus)
oral_t <- t(oral_genus)
skin_t <- t(skin_genus)
nasal_t <- t(nasal_genus)

# 修改样本名称以避免重复，同时保留原始信息
# 假设相同样本名表示来自同一个人的不同部位
rownames(gut_t) <- paste0(rownames(gut_t), "_gut")
rownames(oral_t) <- paste0(rownames(oral_t), "_oral")
rownames(skin_t) <- paste0(rownames(skin_t), "_skin")
rownames(nasal_t) <- paste0(rownames(nasal_t), "_nasal")

# 创建样本类型标记
gut_labels <- data.frame(Sample = rownames(gut_t), Site = "Gut", 
                         Subject = sub("_gut$", "", rownames(gut_t)))
oral_labels <- data.frame(Sample = rownames(oral_t), Site = "Oral", 
                          Subject = sub("_oral$", "", rownames(oral_t)))
skin_labels <- data.frame(Sample = rownames(skin_t), Site = "Skin", 
                          Subject = sub("_skin$", "", rownames(skin_t)))
nasal_labels <- data.frame(Sample = rownames(nasal_t), Site = "Nasal", 
                           Subject = sub("_nasal$", "", rownames(nasal_t)))

# 合并所有物种
# 首先确保所有表格有相同的物种列
all_species <- unique(c(colnames(gut_t), colnames(oral_t), colnames(skin_t), colnames(nasal_t)))

# 修改填充缺失物种的函数，避免索引错误
fill_missing_species <- function(df, all_species) {
  # 创建一个新的数据框，包含所有可能的物种
  result <- matrix(0, nrow = nrow(df), ncol = length(all_species))
  rownames(result) <- rownames(df)
  colnames(result) <- all_species
  
  # 填充现有数据
  common_species <- intersect(colnames(df), all_species)
  for (sp in common_species) {
    result[, sp] <- df[, sp]
  }
  
  # 转换为数据框并返回
  return(as.data.frame(result))
}

# 应用修改后的函数
gut_complete <- fill_missing_species(gut_t, all_species)
oral_complete <- fill_missing_species(oral_t, all_species)
skin_complete <- fill_missing_species(skin_t, all_species)
nasal_complete <- fill_missing_species(nasal_t, all_species)

# 合并所有样本数据
all_data <- rbind(gut_complete, oral_complete, skin_complete, nasal_complete)

# 合并样本标签
sample_metadata <- rbind(gut_labels, oral_labels, skin_labels, nasal_labels)
rownames(sample_metadata) <- sample_metadata$Sample

# 确保样本顺序匹配
sample_metadata <- sample_metadata[rownames(all_data), ]

# 计算Bray-Curtis距离
bray_dist <- vegdist(all_data, method = "bray")

# 进行NMDS分析
set.seed(123)  # 设置随机种子以确保结果可重复
nmds_result <- metaMDS(bray_dist, k = 2, trymax = 100, autotransform = FALSE)

# 检查NMDS分析是否收敛，并输出应力值(stress)
cat("NMDS Stress:", nmds_result$stress, "\n")
if(nmds_result$stress > 0.2) {
  warning("NMDS 应力值 > 0.2，表明排序质量较差")
} else if(nmds_result$stress > 0.1) {
  cat("NMDS 应力值在0.1~0.2之间，表明排序质量一般\n")
} else {
  cat("NMDS 应力值 < 0.1，表明排序质量良好\n")
}

# 提取NMDS坐标
nmds_df <- as.data.frame(nmds_result$points)
colnames(nmds_df) <- c("NMDS1", "NMDS2")

# 将样本信息添加到NMDS数据
nmds_df$Sample <- rownames(nmds_df)
nmds_df <- merge(nmds_df, sample_metadata, by = "Sample")

# 绘制NMDS图
nmds_plot <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Site, shape = Site)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = Site), level = 0.95, linetype = 2) +
  scale_color_manual(values = c("Gut" = "#FF5733", "Oral" = "#33A8FF", 
                                "Skin" = "#33FF57", "Nasal" = "#D433FF")) +
  labs(title = "NMDS of Microbiome Samples Based on Bray-Curtis Distance",
       subtitle = paste("Stress =", round(nmds_result$stress, 3)),
       x = "NMDS1",
       y = "NMDS2") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

# 为每个数据集添加来源标签
gut_samples <- colnames(gut_genus)
oral_samples <- colnames(oral_genus)
skin_samples <- colnames(skin_genus)
nasal_samples <- colnames(nasal_genus)

# 整合所有数据
# 转置矩阵使行为样本，列为物种
gut_t <- t(gut_genus)
oral_t <- t(oral_genus)
skin_t <- t(skin_genus)
nasal_t <- t(nasal_genus)

# 修改样本名称以避免重复，同时保留原始信息
# 假设相同样本名表示来自同一个人的不同部位
rownames(gut_t) <- paste0(rownames(gut_t), "_gut")
rownames(oral_t) <- paste0(rownames(oral_t), "_oral")
rownames(skin_t) <- paste0(rownames(skin_t), "_skin")
rownames(nasal_t) <- paste0(rownames(nasal_t), "_nasal")

# 创建样本类型标记
gut_labels <- data.frame(Sample = rownames(gut_t), Site = "Gut", 
                         Subject = sub("_gut$", "", rownames(gut_t)))
oral_labels <- data.frame(Sample = rownames(oral_t), Site = "Oral", 
                          Subject = sub("_oral$", "", rownames(oral_t)))
skin_labels <- data.frame(Sample = rownames(skin_t), Site = "Skin", 
                          Subject = sub("_skin$", "", rownames(skin_t)))
nasal_labels <- data.frame(Sample = rownames(nasal_t), Site = "Nasal", 
                           Subject = sub("_nasal$", "", rownames(nasal_t)))

# 合并所有物种
# 首先确保所有表格有相同的物种列
all_species <- unique(c(colnames(gut_t), colnames(oral_t), colnames(skin_t), colnames(nasal_t)))

# 修改填充缺失物种的函数，避免索引错误
fill_missing_species <- function(df, all_species) {
  # 创建一个新的数据框，包含所有可能的物种
  result <- matrix(0, nrow = nrow(df), ncol = length(all_species))
  rownames(result) <- rownames(df)
  colnames(result) <- all_species
  
  # 填充现有数据
  common_species <- intersect(colnames(df), all_species)
  for (sp in common_species) {
    result[, sp] <- df[, sp]
  }
  
  # 转换为数据框并返回
  return(as.data.frame(result))
}

# 应用修改后的函数
gut_complete <- fill_missing_species(gut_t, all_species)
oral_complete <- fill_missing_species(oral_t, all_species)
skin_complete <- fill_missing_species(skin_t, all_species)
nasal_complete <- fill_missing_species(nasal_t, all_species)

# 合并所有样本数据
all_data <- rbind(gut_complete, oral_complete, skin_complete, nasal_complete)

# 合并样本标签
sample_metadata <- rbind(gut_labels, oral_labels, skin_labels, nasal_labels)
rownames(sample_metadata) <- sample_metadata$Sample

# 确保样本顺序匹配
sample_metadata <- sample_metadata[rownames(all_data), ]

# 计算Bray-Curtis距离
bray_dist <- vegdist(all_data, method = "bray")

# 进行NMDS分析
set.seed(123)  # 设置随机种子以确保结果可重复
nmds_result <- metaMDS(bray_dist, k = 2, trymax = 100, autotransform = FALSE)

# 检查NMDS分析是否收敛，并输出应力值(stress)
cat("NMDS Stress:", nmds_result$stress, "\n")
if(nmds_result$stress > 0.2) {
  warning("NMDS 应力值 > 0.2，表明排序质量较差")
} else if(nmds_result$stress > 0.1) {
  cat("NMDS 应力值在0.1~0.2之间，表明排序质量一般\n")
} else {
  cat("NMDS 应力值 < 0.1，表明排序质量良好\n")
}

# 提取NMDS坐标
nmds_df <- as.data.frame(nmds_result$points)
colnames(nmds_df) <- c("NMDS1", "NMDS2")

# 将样本信息添加到NMDS数据
nmds_df$Sample <- rownames(nmds_df)
nmds_df <- merge(nmds_df, sample_metadata, by = "Sample")

# 绘制NMDS图
ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, fill = Site,color=Site)) +
  geom_point(size = 4, alpha = 0.8,colour = "white",shape = 21) +
  stat_ellipse(aes(group = Site), level = 0.95, linetype = 2) +
  scale_fill_manual(values = c("Gut" = "#edd064", "Oral" = "#a1d5b9", 
                               "Skin" = "#f2ccac", "Nasal" = "#a17db4")) +
  scale_color_manual(values = c("Gut" = "#edd064", "Oral" = "#a1d5b9", 
                                "Skin" = "#f2ccac", "Nasal" = "#a17db4"))+
  labs(title = "NMDS of Microbiome Samples",
       subtitle = paste("Stress =", round(nmds_result$stress, 3)),
       x = "NMDS1",
       y = "NMDS2") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 12),
    panel.grid.minor = element_blank()
  )


## 统计代谢物检测的种类


# 加载必要的包
library(ggplot2)
library(dplyr)

class_counts <- as.data.frame(table(metabolite_annotation$HMDB.Class,metabolite_annotation$HMDB.Source.Microbial))
colnames(class_counts) <- c("HMDB.Class","HMDB.Source.Microbial", "Count")


class_counts<-subset(class_counts,!(HMDB.Class=="NA"))
class_counts<-subset(class_counts,!(HMDB.Source.Microbial=="NA"))
class_counts<-subset(class_counts,Count>=5)

total_by_class <- aggregate(Count ~ HMDB.Class, data = class_counts, sum)
# 按总计数从小到大排序
class_order <- total_by_class$HMDB.Class[order(total_by_class$Count)]
# 将 HMDB.Class 转换为有序因子
class_counts$HMDB.Class <- factor(class_counts$HMDB.Class, levels = class_order)

# 然后绘图
ggplot(class_counts, aes(x = HMDB.Class, y = Count, fill=HMDB.Source.Microbial)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#597c8b", "#b36a6f"))+
  coord_flip() +
  labs(x = "HMDB.Class", y = "Counts") +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12)
  )
### 来源维恩图


metabolite_annotation_Source<-metabolite_annotation[,29:31]

metabolite_annotation_Source<-subset(metabolite_annotation_Source,!(HMDB.Source.Endogenous=="NA"))
metabolite_annotation_Source<-subset(metabolite_annotation_Source,!(HMDB.Source.Food=="NA"))
metabolite_annotation_Source<-subset(metabolite_annotation_Source,!(HMDB.Source.Microbial=="NA"))


library(ggVennDiagram)
library(ggplot2)

# Assuming your data is already loaded as metabolite_annotation
# If not, use the following code to read it (replace with your file path)
# metabolite_annotation <- read.csv("path_to_your_file.csv", stringsAsFactors = FALSE)

# Create lists of metabolites for each source
endogenous_metabolites <- rownames(metabolite_annotation)[metabolite_annotation$HMDB.Source.Endogenous == TRUE]
food_metabolites <- rownames(metabolite_annotation)[metabolite_annotation$HMDB.Source.Food == TRUE]
microbial_metabolites <- rownames(metabolite_annotation)[metabolite_annotation$HMDB.Source.Microbial == TRUE]

# Create a list for the Venn diagram
venn_list <- list(
  Endogenous = endogenous_metabolites,
  Food = food_metabolites,
  Microbial = microbial_metabolites
)

# Create the Venn diagram
venn_plot <- ggVennDiagram(venn_list, 
                           label = "count",
                           category.names = c("Endogenous", "Food", "Microbial"))

# Customize the plot
venn_plot <- venn_plot + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none") +
  labs(title = "Distribution of Metabolite Sources",
       caption = "Data source: metabolite_annotation")

# Display the plot
