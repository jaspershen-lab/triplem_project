rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")
library(tidyverse)
library(tidymass)
library(readxl)



demographic_data<-read.table("2_data/FBIP-main/metadata/FBIP_metadata.txt")

# 筛选时间点0
demographic_data<-subset(demographic_data,Time=="W0")


library(ComplexHeatmap)
library(circlize)
library(grid)


demographic_data$Gender<-gsub("F","Female",demographic_data$Gender)
demographic_data$Gender<-gsub("M","Male",demographic_data$Gender)
# 将不透明度设置为0.7（相当于30%的透明度）
alpha_value = 0.7

# 创建注释对象
ha = columnAnnotation(
  # 连续型变量使用单一颜色，添加透明度
  Age = anno_barplot(demographic_data$Age, 
                     gp = gpar(fill = scales::alpha("#E69F00", alpha_value)),
                     border = TRUE,
                     width = unit(2, "cm")),  # 使用橙色
  BMI = anno_barplot(demographic_data$BMI,
                     gp = gpar(fill = scales::alpha("#56B4E9", alpha_value)),
                     border = TRUE,
                     width = unit(2, "cm")),  # 使用蓝色
  TC = anno_barplot(demographic_data$TC,
                     gp = gpar(fill = scales::alpha("#F21451", alpha_value)),
                     border = TRUE,
                     width = unit(2, "cm")),  
  Insulin = anno_barplot(demographic_data$Insulin,
                    gp = gpar(fill = scales::alpha("#F45112", alpha_value)),
                    border = TRUE,
                    width = unit(2, "cm")),  
  
  border = TRUE,
  # 分类型变量
  Gender = demographic_data$Gender,

  # 设置分类变量的颜色，添加透明度
  col = list(
    Gender = mapply(function(x) scales::alpha(x, alpha_value), sex_color, USE.NAMES = TRUE)
  ),
  
  # 设置注释的样式
  annotation_name_gp = gpar(fontsize = 10),
  annotation_name_side = "left",
  simple_anno_size = unit(0.5, "cm"),
  
  # 添加图例
  show_legend = TRUE,
  annotation_legend_param = list(
    Gender = list(title = "Gender")
  ),
  
  # 设置列间距
  gap = unit(c(2, 2, 1, 1, 1), "mm")  # 在各列之间添加间距
)

# 创建一个空矩阵用于绘制热图
mat = matrix(0, nrow = 1, ncol = nrow(demographic_data))

# 创建并绘制热图，禁用聚类
ht = Heatmap(mat,
             top_annotation = ha,
             show_row_names = FALSE,
             show_column_names = FALSE,
             show_heatmap_legend = FALSE,
             cluster_rows = FALSE,    # 禁用行聚类
             cluster_columns = FALSE, # 禁用列聚类
             height = unit(0.2, "cm"))   # 设置整体宽度

# 绘制并设置图例位置



pdf("4_manuscript/Figures/Figure_1/figure_s1_FBIP_demographic.pdf", width = 10, height = 6)
ht_drawn <-draw(ht, annotation_legend_side = "bottom")  # 捕获绘制的对象
dev.off()