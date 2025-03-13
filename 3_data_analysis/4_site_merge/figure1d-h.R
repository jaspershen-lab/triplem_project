rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")
source("1_code/mantel_Procrustes_code.R")
library(tidyverse)
library(tidymass)
library(readxl)



load("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section")
metabolomics_object<-object_cross_section


setwd("3_data_analysis/4_site_merge/")

### 计算四个身体部位菌群的alpha多样性

metabolite_annotation<-read_excel("../../3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx")

load("../../3_data_analysis/gut_microbiome/data_preparation/object_cross_section")

gut_object<-object_cross_section


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







shannon_div <- diversity(t(gut_object@expression_data), index = "shannon")

# 创建结果数据框
results_gut <- data.frame(
  Sample = names(shannon_div),
  Shannon = shannon_div
)


load("../../3_data_analysis/oral_microbiome/data_preparation/object_cross_section")

oral_object<-object_cross_section


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







shannon_div <- diversity(t(oral_object@expression_data), index = "shannon")

# 创建结果数据框
results_oral <- data.frame(
  Sample = names(shannon_div),
  Shannon = shannon_div
)


load("../../3_data_analysis/skin_microbiome/data_preparation/object_cross_section")

skin_object<-object_cross_section


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







shannon_div <- diversity(t(skin_object@expression_data), index = "shannon")

# 创建结果数据框
results_skin <- data.frame(
  Sample = names(shannon_div),
  Shannon = shannon_div
)


load("../../3_data_analysis/nasal_microbiome/data_preparation/object_cross_section")

nasal_object<-object_cross_section


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







shannon_div <- diversity(t(nasal_object@expression_data), index = "shannon")

# 创建结果数据框
results_nasal <- data.frame(
  Sample = names(shannon_div),
  Shannon = shannon_div
)


# 合并alpha diversity

Sample_ID<-data.frame(metabolomics_object@sample_info)[,1:2]
colnames(Sample_ID)[1]<-"Sample"

alpha_diversity<-Sample_ID%>%
  full_join(results_gut, by = "Sample") %>%
  full_join(results_oral, by = "Sample") %>%
  full_join(results_skin, by = "Sample") %>%
  full_join(results_nasal, by = "Sample")

rownames(alpha_diversity)<-alpha_diversity$Sample

alpha_diversity<-alpha_diversity[,-1:-2]

colnames(alpha_diversity)<-c("gut","oral","skin","nasal")






## 提取metabolite数据

#######adjust BMI, sex, and IRIS, ethnicity
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

metabolome_data<-metabolomics_temp_object@expression_data


metabolome_data<-data.frame(t(metabolome_data))

#筛选共有的样本

metabolome_data<-metabolome_data[rownames(alpha_diversity),]


# 加载必要的R包
library(tidyverse)



# 确保样本ID是匹配的
common_samples <- intersect(rownames(alpha_diversity), rownames(metabolome_data))
if(length(common_samples) == 0) {
  stop("两个数据集没有共同的样本ID")
}
cat("共有", length(common_samples), "个样本可用于分析\n")

# 使用共同样本筛选数据
alpha_diversity_filtered <- alpha_diversity[common_samples, ]
metabolome_data_filtered <- metabolome_data[common_samples, ]

# 计算每个位点alpha多样性与代谢物的Spearman相关性
sites <- colnames(alpha_diversity_filtered)
metabolites <- colnames(metabolome_data_filtered)

# 准备存储结果的数据框
results <- data.frame(
  Site = character(),
  Metabolite = character(),
  Rho = numeric(),
  P_value = numeric(),
  Adjusted_P_value = numeric(),
  stringsAsFactors = FALSE
)

# 计算相关性
for(site in sites) {
  alpha_values <- alpha_diversity_filtered[[site]]
  
  for(metabolite in metabolites) {
    metabolite_values <- metabolome_data_filtered[[metabolite]]
    
    # 计算Spearman相关性
    cor_test <- cor.test(alpha_values, metabolite_values, 
                         method = "pearson", 
                         exact = FALSE)
    
    # 添加结果到数据框
    results <- rbind(results, data.frame(
      Site = site,
      Metabolite = metabolite,
      Rho = cor_test$estimate,
      P_value = cor_test$p.value,
      Adjusted_P_value = NA,  # 先设为NA，后面统一调整
      stringsAsFactors = FALSE
    ))
  }
}

# 对p值进行BH校正
results$Adjusted_P_value <- p.adjust(results$P_value, method = "BH")

results<-subset(results,P_value<0.05)

results<-merge(results,metabolite_annotation[,c("variable_id","HMDB.Name","HMDB.Class","HMDB.Source.Microbial")],by.x="Metabolite",by.y="variable_id")


results<-subset(results,!(HMDB.Class=="NA"))


# 定义保留的类别
keep_classes <- c("Benzene and substituted derivatives", 
                  "Carboxylic acids and derivatives",
                  "Fatty Acyls",
                  "Glycerophospholipids",
                  "Organic sulfuric acids and derivatives",
                  "Organooxygen compounds",
                  "Piperidines",
                  "Steroids and steroid derivatives")

# 将不在保留类别列表中的HMDB.Class值重新分类为"Others"
results$HMDB.Class <- ifelse(results$HMDB.Class %in% keep_classes, 
                             results$HMDB.Class, 
                             "Others")



# Updated visualization for metabolite correlations by HMDB Class using original Rho values
library(ggplot2)
library(dplyr)
library(scales)
library(ggrepel)
library(patchwork) # For combining plots

# Data preprocessing
plot_data <- results %>%
  # Use original Rho values
  mutate(
    # Transform p-values for better visualization
    neg_log_p = -log10(P_value),
    # Size based on p-value significance
    p_size = case_when(
      P_value < 0.001 ~ 4,
      P_value < 0.01 ~ 3,
      P_value < 0.05 ~ 2,
      TRUE ~ 1
    )
  ) %>%
  # Calculate the number of metabolites per class and filter
  group_by(HMDB.Class) %>%
  dplyr::mutate(
    class_size = n(),
    class_median_rho = median(Rho, na.rm = TRUE),
    class_pct_significant = mean(Adjusted_P_value < 0.05, na.rm = TRUE) * 100,
    # Rank within each class by absolute correlation (for labeling purposes)
    class_rank = rank(-abs(Rho))
  ) %>%
  ungroup() %>%
  # Filter classes with at least 10 members (adjust this threshold as needed)
  filter(class_size >= 10) %>%
  # Sort by median Rho within class
  arrange(HMDB.Class, desc(Rho)) %>%
  # Add significance markers and labels
  mutate(
    # Label top 3 metabolites in each class
    label = ifelse(class_rank <= 3 & Adjusted_P_value < 0.05, HMDB.Name, "")
  )

# Assign sequential numbers for x-axis
plot_data$metabolite_sort <- 1:nrow(plot_data)

# Calculate class ranges for background shading
class_summary <- plot_data %>%
  group_by(HMDB.Class) %>%
  dplyr::summarise(
    start = min(metabolite_sort),
    end = max(metabolite_sort),
    mid = mean(c(min(metabolite_sort), max(metabolite_sort))),
    median_rho = median(Rho, na.rm = TRUE),
    pct_significant = mean(Adjusted_P_value < 0.05, na.rm = TRUE) * 100,
    n = n()
  )

# Custom color palette for sites
site_colors <- c(
  "gut" = "#edd064",
  "oral" = "#a1d5b9",
  "skin" = "#f2ccac",
  "nasal" = "#a17db4"
)

# Main plot
p1 <- ggplot(plot_data, aes(x = metabolite_sort, y = Rho)) +
  # Add alternating backgrounds for classes
  geom_rect(data = class_summary,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = "gray95", alpha = 0.5, inherit.aes = FALSE) +
  # Add zero line to distinguish positive from negative correlations
  geom_hline(yintercept = 0, linetype = "solid", color = "gray50", size = 0.7) +
  # Add points with site colors
  geom_point(aes(color = Site, 
                 size = p_size), 
             alpha = 0.85) +
  # Add FDR significance thresholds (both for positive and negative correlations)
  geom_hline(yintercept = median(subset(plot_data, Adjusted_P_value == 0.05 & Rho > 0)$Rho, na.rm = TRUE), 
             linetype = "dashed", color = "darkred", size = 0.5) +
  geom_hline(yintercept = median(subset(plot_data, Adjusted_P_value == 0.05 & Rho < 0)$Rho, na.rm = TRUE), 
             linetype = "dashed", color = "darkblue", size = 0.5) +
  # Add labels for top metabolites
  geom_text_repel(
    data = subset(plot_data, label != ""),
    aes(label = label),
    size = 3.5,
    box.padding = 0.4,
    point.padding = 0.3,
    force = 8,
    max.overlaps = 30,
    segment.color = "grey50",
    segment.size = 0.2,
    min.segment.length = 0.1
  ) +
  # Customize scales
  scale_color_manual(values = site_colors, name = "Site") +
  scale_size_continuous(
    name = "P-value", 
    breaks = c(1, 2, 3, 4),
    labels = c("ns", "P < 0.05", "P < 0.01", "P < 0.001"),
    range = c(2, 5)
  ) +
  scale_y_continuous(
    breaks = seq(-1, 1, by = 0.1),
    minor_breaks = seq(-1, 1, by = 0.05),
    limits = c(min(plot_data$Rho) * 1.05, max(plot_data$Rho) * 1.05),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  # Customize theme
  theme_minimal() +
  theme(
    panel.grid.minor = element_line(color = "gray95"),
    panel.grid.major = element_line(color = "gray90"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10),
    plot.caption = element_text(hjust = 0, size = 8)
  ) +
  # Add class names as x-axis labels
  scale_x_continuous(
    breaks = class_summary$mid,
    labels = class_summary$HMDB.Class,
    expand = c(0.01, 0.01)
  ) +
  # Add informative labels
  labs(
    x = "HMDB Class", 
    y = "Correlation (Rho)",
    title = "Metabolite Correlations by HMDB Class",
    subtitle = paste0("Showing ", nrow(plot_data), " metabolites across ", 
                      length(unique(plot_data$HMDB.Class)), " HMDB classes")
  )

# Display the plot
p1



library(ggpubr)

# 绘制是否为微生物来源代谢物的Rho大小

plot_microbial_source_data<-plot_data

plot_microbial_source_data$HMDB.Source.Microbial<-gsub("NA","FALSE",plot_microbial_source_data$HMDB.Source.Microbial)

HMDB.Source.Microbial_color<-c("","")

ggplot(data=plot_microbial_source_data,aes(x=HMDB.Source.Microbial,y=abs(Rho),fill=HMDB.Source.Microbial))+ #”fill=“设置填充颜色
  stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+ #由于自带的箱形图没有胡须末端没有短横线，使用误差条的方式补上
  geom_boxplot(size=0.5,outlier.fill="white",outlier.color="white")+ #size设置箱线图的边框线和胡须的线宽度，fill设置填充颜色，outlier.fill和outlier.color设置异常点的属性
  geom_jitter(fill="black",width =0.2,shape = 21,size=2)+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
  scale_fill_manual(values = c( "#597c8b","#b36a6f"))+
  scale_color_manual(values=c("black","black","black"))+ #设置散点图的圆圈的颜色为黑色
  ggtitle("")+ #设置总的标题
  theme_bw()+ theme(legend.position="none", #不需要图例
                    axis.text.x=element_text(colour="black",size=14,"Helvetica"), #设置x轴刻度标签的字体属性
                    axis.text.y=element_text(size=14,family = "Helvetica"), #设置x轴刻度标签的字体属性
                    axis.title.y=element_text(size = 14,family = "Helvetica"), #设置y轴的标题的字体属性
                    axis.title.x=element_text(size = 14,family = "Helvetica"), #设置x轴的标题的字体属性
                    plot.title = element_text(size=15,face="bold","Helvetica",hjust = 0.5))+
  ylab("|Rho|")+xlab("Microbial Source")+ #设置x轴和y轴的标题
  stat_compare_means() 

## 绘制四个身体部位的rho值山峦图

# Ridge plot (mountain plot) of Rho values by microbiome site
library(ggplot2)
library(dplyr)
library(ggridges)

# Custom color palette for sites


# Data preparation
ridge_data <- results %>%
  # If needed, ensure Site is a factor with desired order
  mutate(
    Site = factor(Site, levels = c("gut", "oral", "skin", "nasal")),
    # Add significance flag
    is_significant = Adjusted_P_value < 0.05
  )



# Create ridge plot
p1 <- ggplot(ridge_data, aes(x = abs(Rho), y = Site, fill = Site)) +
  # Add density ridges
  geom_density_ridges(
    aes(height = ..density..),
    alpha = 0.8,
    scale = 3,
    rel_min_height = 0.01,
    quantile_lines = TRUE,
    quantiles = 1
  ) +
  # Customize colors
  scale_fill_manual(values = site_colors) +

  labs(
    title = "",
    subtitle = "",
    x = "Spearman Correlation Coefficient (Rho)",
    y = NULL
  ) +
  # Customize theme
  theme_ridges(center_axis_labels = TRUE) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.title.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    plot.caption = element_text(hjust = 0, size = 8)
  )




# 统计总的代谢物种类的柱状图


# 读取数据
# 假设你的数据在名为"metabolite_annotation.csv"的文件中
# 如果数据在不同位置，请修改文件路径
data <- metabolite_annotation



# 处理空值（NA或空字符串）
data<-subset(data,!(HMDB.Class=="NA"))

# 加载所需的库
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(forcats)

# 假设我们有两个数据框：
# 1. metabolite_annotation - 包含总体代谢物数据
# 2. plot_data - 包含各位点的代谢物数据

# 第一步：处理总体代谢物数据
total_data <- metabolite_annotation
total_data <- subset(total_data, !(HMDB.Class=="NA"))

# 统计每个类别的数量和百分比
class_summary <- total_data %>%
  dplyr::count(HMDB.Class) %>%  # 修正之前的dpcount为count
  mutate(
    total = sum(n),
    percentage = n / total * 100
  ) %>%
  arrange(desc(percentage))

# 只保留前7个最常见的类别，其余归为"Others"
top_n <- 7
if(nrow(class_summary) > top_n) {
  top_classes <- class_summary[1:top_n,]
  others <- data.frame(
    HMDB.Class = "Others",
    n = sum(class_summary$n[(top_n+1):nrow(class_summary)]),
    total = class_summary$total[1],
    percentage = sum(class_summary$percentage[(top_n+1):nrow(class_summary)])
  )
  class_summary <- rbind(top_classes, others)
}

# 添加一个"Site"列，标记为"Total"
class_summary <- class_summary %>%
  mutate(Site = "Total")

# 确保Site列已经添加成功
print(names(class_summary))  # 检查列名是否包含Site

# 第二步：处理各位点代谢物数据
site_data <- plot_data
site_data$Site <- tolower(site_data$Site)
expected_sites <- c("gut", "oral", "skin", "nasal")
site_data <- site_data %>% filter(Site %in% expected_sites)

# 提取总体数据中的类别，以便在位点数据中保持一致
important_classes <- class_summary$HMDB.Class
if("Others" %in% important_classes) {
  important_classes <- important_classes[important_classes != "Others"]
}

# 统计各位点每个类别的数量和百分比
site_class_summary <- site_data %>%
  group_by(Site, HMDB.Class) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  group_by(Site) %>%
  mutate(
    total = sum(count),
    percentage = count / total * 100
  )

# 处理各位点的数据，确保类别一致性
site_class_final <- site_class_summary %>%
  # 标记是否属于重要类别
  mutate(is_important = HMDB.Class %in% important_classes) %>%
  group_by(Site) %>%
  # 对于每个位点，将非重要类别归为"Others"
  mutate(
    HMDB.Class = if_else(is_important, HMDB.Class, "Others")
  ) %>%
  # 重新统计"Others"类别
  group_by(Site, HMDB.Class) %>%
  dplyr::summarise(
    percentage = sum(percentage),
    .groups = "drop"
  )

# 检查结果
print("site_class_final structure:")
print(str(site_class_final))
print("site_class_final first few rows:")
print(head(site_class_final))

# 如果结果不正确，尝试替代方法
if(ncol(site_class_final) < 3 || !"Site" %in% names(site_class_final)) {
  print("使用替代方法重新计算site_class_final")
  
  # 替代方法：逐步执行并检查中间结果
  site_data_with_important <- site_class_summary %>%
    mutate(is_important = HMDB.Class %in% important_classes)
  
  print("Intermediate data:")
  print(head(site_data_with_important))
  
  # 手动重分类
  site_data_reclassified <- site_data_with_important %>%
    mutate(HMDB.Class_new = ifelse(is_important, as.character(HMDB.Class), "Others"))
  
  # 按新分类汇总
  site_class_final <- site_data_reclassified %>%
    group_by(Site, HMDB.Class_new) %>%
    summarise(
      percentage = sum(percentage),
      .groups = "drop"
    ) %>%
    rename(HMDB.Class = HMDB.Class_new)
  
  print("New site_class_final:")
  print(head(site_class_final))
}

# 第三步：合并总体和各位点数据
# 检查两个数据框的结构
print("class_summary columns:")
print(names(class_summary))
print("site_class_final columns:")
print(names(site_class_final))

# 确保两个数据框有相同的列名
class_summary_selected <- class_summary %>% 
  dplyr::select(Site, HMDB.Class, percentage)
site_class_selected <- site_class_final %>% 
  dplyr::select(Site, HMDB.Class, percentage)

# 合并数据
combined_data <- rbind(
  class_summary_selected,
  site_class_selected
)

# 确保"Others"类别在图例中排在最后
combined_data$HMDB.Class <- fct_relevel(as.factor(combined_data$HMDB.Class), "Others", after = Inf)

# 设置Site的顺序，使"Total"在最左边
combined_data$Site <- factor(combined_data$Site, levels = c("Total", expected_sites))

# 生成足够的颜色
unique_classes <- levels(combined_data$HMDB.Class)
num_colors <- length(unique_classes)

# 使用更多样化的调色板
if(num_colors <= 8) {
  colors <- brewer.pal(max(3, num_colors), "Set1")
} else if(num_colors <= 12) {
  colors <- brewer.pal(max(3, num_colors), "Paired")
} else {
  # 组合多个调色板以获得更多颜色
  colors <- c(
    brewer.pal(8, "Set1"),
    brewer.pal(8, "Set2"),
    brewer.pal(8, "Set3")
  )
  # 如果还不够，则使用颜色渐变
  if(length(colors) < num_colors) {
    colors <- colorRampPalette(colors)(num_colors)
  } else {
    colors <- colors[1:num_colors]
  }
}

# 如果有第6个颜色且需要删除
if(length(colors) >= 6) {
  colors <- colors[-6]
}

# 打印颜色和类别的长度进行诊断
print(paste("Length of colors:", length(colors)))
print(paste("Length of unique_classes:", length(unique_classes)))

# 确保颜色和类别数量一致
if(length(colors) != length(unique_classes)) {
  # 如果颜色少于类别，添加额外的颜色
  if(length(colors) < length(unique_classes)) {
    additional_colors_needed <- length(unique_classes) - length(colors)
    additional_colors <- colorRampPalette(colors)(additional_colors_needed)
    colors <- c(colors, additional_colors)
  } 
  # 如果颜色多于类别，截断颜色向量
  else {
    colors <- colors[1:length(unique_classes)]
  }
}

# 再次检查长度
print(paste("Adjusted length of colors:", length(colors)))

# 确保"Others"类别使用灰色
if("Others" %in% unique_classes) {
  names(colors) <- unique_classes  # 现在长度应该匹配了
  colors["Others"] <- "#999999" # 灰色
}
combined_data$Site<-factor(combined_data$Site,levels = c("Total","gut","oral","skin","nasal"))
# 创建堆叠柱状图
p <- ggplot(combined_data, aes(x = Site, y = percentage, fill = HMDB.Class)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(
    title = "",
    x = "Site",
    y = "Percent (%)",
    fill = "HMDB.Class"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  ) +
  scale_y_continuous(breaks = seq(0, 100, 10))

# 为每个分段添加标签（当分段足够大时）
# 首先计算每个分段的位置
label_data <- combined_data %>%
  group_by(Site) %>%
  arrange(desc(HMDB.Class)) %>%
  mutate(
    ymin = lag(cumsum(percentage), default = 0),
    ymax = cumsum(percentage),
    pos = (ymin + ymax) / 2
  ) %>%
  ungroup()

# 只为较大的分段添加标签（例如>=5%）
label_threshold <- 5
p <- p +
  geom_text(
    data = subset(label_data, percentage >= label_threshold),
    aes(
      x = Site,
      y = pos,
      label = sprintf("%.1f%%", percentage)
    ),
    size = 3,
    color = "white"
  )

# 最终设置主题和字体
p + theme(
  axis.text.x = element_text(colour = "black", size = 14, family = "Helvetica"), 
  axis.text.y = element_text(size = 14, family = "Helvetica"), 
  axis.title.y = element_text(size = 14, family = "Helvetica"), 
  axis.title.x = element_text(size = 14, family = "Helvetica"), 
  plot.title = element_text(size = 15, face = "bold", family = "Helvetica", hjust = 0.5)
)

# 加载所需的库
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)

# 使用site_class_final进行分析
# 首先，确保Total位点也包含在数据中
# 如果不包含，您需要添加它

# 检查数据中包含哪些位点
print(unique(site_class_final$Site))

# 1. 准备卡方检验的数据
# 创建一个函数进行两两位点间的卡方检验
chi_square_test <- function(site1, site2, data) {
  # 提取两个位点的数据
  site_data <- data %>% 
    filter(Site %in% c(site1, site2))
  
  # 将数据转换为宽格式，以便进行卡方检验
  # 使用实际的count列进行检验
  contingency_table <- site_data %>%
    select(Site, HMDB.Class, count) %>%
    spread(key = Site, value = count, fill = 0)  # 使用0填充缺失值
  
  # 打印检查转换后的频数表
  print(paste("Contingency table for", site1, "vs", site2))
  print(contingency_table)
  
  # 删除HMDB.Class列，只保留频数用于卡方检验
  freq_table <- contingency_table %>% 
    select(-HMDB.Class) %>%
    as.matrix()
  
  # 执行卡方检验
  chi_test <- tryCatch({
    chisq.test(freq_table)
  }, error = function(e) {
    # 如果出错（例如频数过低），返回NA值
    warning(paste("Error in chi-square test for", site1, "vs", site2, ":", e$message))
    return(list(
      p.value = NA,
      parameter = NA,
      statistic = NA
    ))
  })
  
  # 返回p值和自由度
  return(list(
    site1 = site1,
    site2 = site2,
    p_value = chi_test$p.value,
    df = chi_test$parameter,
    statistic = chi_test$statistic
  ))
}

# 2. 获取所有位点
sites <- unique(site_class_final$Site)
site_pairs <- expand.grid(site1 = sites, site2 = sites, stringsAsFactors = FALSE)
# 删除自己与自己比较的行
site_pairs <- site_pairs %>% filter(site1 != site2)

# 3. 执行所有对的卡方检验
results <- list()
for(i in 1:nrow(site_pairs)) {
  pair <- site_pairs[i, ]
  cat("Testing", pair$site1, "vs", pair$site2, "\n")
  result <- chi_square_test(pair$site1, pair$site2, site_class_final)
  results[[i]] <- result
}

# 4. 将结果转换为数据框
result_df <- do.call(rbind, lapply(results, function(x) {
  data.frame(
    site1 = x$site1,
    site2 = x$site2,
    p_value = x$p_value,
    statistic = x$statistic,
    df = x$df,
    stringsAsFactors = FALSE
  )
}))

# 打印结果表格
print(result_df)

# 5. 创建p值矩阵用于热图显示
p_value_matrix <- matrix(NA, length(sites), length(sites))
rownames(p_value_matrix) <- sites
colnames(p_value_matrix) <- sites

# 填充矩阵
for(i in 1:nrow(result_df)) {
  row_idx <- which(sites == result_df$site1[i])
  col_idx <- which(sites == result_df$site2[i])
  p_value_matrix[row_idx, col_idx] <- result_df$p_value[i]
}

# 对角线设为1（自己与自己比较）
diag(p_value_matrix) <- 1

# 6. 创建热图
# 转换矩阵为长格式用于ggplot
p_value_long <- melt(p_value_matrix)
names(p_value_long) <- c("Site1", "Site2", "p_value")

# 添加显著性标记
p_value_long$significance <- "NS"
p_value_long$significance[p_value_long$p_value < 0.05] <- "*"
p_value_long$significance[p_value_long$p_value < 0.01] <- "**"
p_value_long$significance[p_value_long$p_value < 0.001] <- "***"
# NA值表示无法计算
p_value_long$significance[is.na(p_value_long$p_value)] <- "NA"

# 创建热图，-log10转换p值
p_value_long$neg_log_p <- -log10(p_value_long$p_value)
# 把NA和Inf替换为0
p_value_long$neg_log_p[is.na(p_value_long$neg_log_p) | is.infinite(p_value_long$neg_log_p)] <- 0

# 创建热图
p <- ggplot(p_value_long, aes(x = Site2, y = Site1, fill = neg_log_p)) +
  geom_tile(color = "white") +
  geom_text(aes(label = significance), color = "black", size = 5) +
  scale_fill_gradient2(
    low = "white", 
    high = "red", 
    mid = "pink",
    midpoint = 1,
    limit = c(0, max(p_value_long$neg_log_p, na.rm = TRUE)),
    name = "-log10(p-value)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  labs(
    title = "各位点代谢物分布差异的卡方检验结果",
    x = "",
    y = ""
  ) +
  coord_fixed()

# 7. 显示热图
print(p)

# 8. 创建显著性摘要表
significance_summary <- result_df %>%
  mutate(
    significance = case_when(
      is.na(p_value) ~ "NA",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "NS"
    ),
    comparison = paste(site1, "vs", site2)
  ) %>%
  select(comparison, p_value, statistic, df, significance) %>%
  arrange(p_value)

# 打印摘要表
print(significance_summary)