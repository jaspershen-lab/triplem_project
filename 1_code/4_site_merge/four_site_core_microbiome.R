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
