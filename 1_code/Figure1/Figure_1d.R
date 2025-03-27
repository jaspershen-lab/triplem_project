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

# 加载必要的包
library(microbiome)
library(phyloseq)
library(dplyr)
library(ggplot2)
library(gridExtra) # 用于组合多个图形



# 创建OTU表
otu_table <- otu_table(as.matrix(t(all_data)), taxa_are_rows = TRUE)

# 创建样本数据表
sample_data <- sample_data(sample_metadata)

# 创建phyloseq对象
physeq <- phyloseq(otu_table, sample_data)

# 按位点分组
site_list <- c("Gut", "Oral", "Skin", "Nasal")
core_taxa_results <- list()

#------------------------ 第一部分：识别核心物种 ------------------------#

# 为每个位点确定核心物种
for (site in site_list) {
  # 筛选特定位点的样本
  site_samples <- subset_samples(physeq, Site == site)
  
  # 使用microbiome包的core_members函数找出核心物种
  # detection参数：相对丰度阈值
  # prevalence参数：在该组中出现的样本百分比阈值
  core_detection <- 0.001  # 相对丰度至少0.1%
  core_prevalence <- 0.5   # 在至少50%的样本中出现
  
  # 计算核心分类群
  core_taxa <- core_members(site_samples, 
                            detection = core_detection, 
                            prevalence = core_prevalence)
  
  # 如果核心物种超过10个，只取丰度最高的10个
  if (length(core_taxa) > 10) {
    # 计算每个分类群的平均丰度
    taxa_sums <- taxa_sums(site_samples)
    taxa_sums <- taxa_sums[names(taxa_sums) %in% core_taxa]
    taxa_sums <- sort(taxa_sums, decreasing = TRUE)
    core_taxa <- names(taxa_sums)[1:10]
  }
  
  # 储存结果
  core_taxa_results[[site]] <- core_taxa
  
  # 打印结果
  cat("\n", site, "核心物种 (", length(core_taxa), "):\n", sep="")
  print(core_taxa)
}

# 合并所有位点的核心物种
all_core_taxa <- unique(unlist(core_taxa_results))
cat("\n总共识别出", length(all_core_taxa), "个核心物种\n")

#------------------------ 第二部分：使用plot_core可视化 ------------------------#
library(RColorBrewer)
# 创建一个函数来绘制每个位点的核心物种图
plot_site_core <- function(physeq, site, core_taxa) {
  # 筛选特定位点的样本
  site_samples <- subset_samples(physeq, Site == site)
  
  # 转换为相对丰度
  site_samples_rel <- site_samples
  
  # 只保留该位点的核心物种
  site_samples_rel <- prune_taxa(all_core_taxa, site_samples_rel)
  
  
  prevalences <- seq(.05, 1, .02)
  
  detections <- round(10^seq(log10(0.0001), log10(.2), length = 30), 3)
  
  # Also define gray color palette
  gray <- rev(brewer.pal(5, "RdBu"))
  
  # 使用plot_core函数绘图
  # 这将创建一个热图显示核心物种在不同检测和流行阈值下的存在情况
  p <- plot_core(site_samples_rel,
                 plot.type = "heatmap", 
                 colours = gray,
                 prevalences = prevalences, 
                 taxa.order = rev(all_core_taxa),
                 detections = detections) +
    labs(x = "Detection Threshold\n(Relative Abundance (%))") +
    
    #Adjusts axis text size and legend bar height
    theme(axis.text.y= element_text(size=14, face="italic"),
          axis.text.x.bottom=element_text(size=8),
          axis.title = element_text(size=10),
          legend.text = element_text(size=8),
          legend.title = element_text(size=10))
  
  return(p)
}

# 为每个位点绘制核心物种图
plot_list <- list()
for (site in site_list) {
  plot_list[[site]] <- plot_site_core(physeq, site, core_taxa_results[[site]])
}

# 将四个图组合在一起
combined_plot <- grid.arrange(
  plot_list[["Gut"]],
  plot_list[["Oral"]],
  plot_list[["Skin"]],
  plot_list[["Nasal"]],
  nrow = 1
)

p_gut<-plot_list[["Gut"]]+theme(legend.position = "none",axis.text.y  = element_blank())
p_oral<-plot_list[["Oral"]]+theme(legend.position = "none",axis.text.y  = element_blank())
p_skin<-plot_list[["Skin"]]+theme(legend.position = "none",axis.text.y  = element_blank())
p_nasal<-plot_list[["Nasal"]]+theme(legend.position = "none",axis.text.y  = element_blank())


combined_plot <- grid.arrange(
  p_gut,
  p_oral,
  p_skin,
  p_nasal,
  nrow = 1
)

pdf("../../4_manuscript/Figures/Figure_1/figure_1d.pdf", width = 10, height = 6)
plot(combined_plot)
dev.off()
