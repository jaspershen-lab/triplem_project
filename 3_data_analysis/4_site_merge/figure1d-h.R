rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")
source("1_code/mantel_Procrustes_code.R")
library(tidyverse)
library(tidymass)
library(readxl)


metabolite_annotation<-read_excel("../../3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx")
load("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section")
metabolomics_object<-object_cross_section


setwd("3_data_analysis/4_site_merge/")

### 计算四个身体部位菌群的alpha多样性



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



# Updated visualization for metabolite correlations by HMDB Class
library(ggplot2)
library(dplyr)
library(scales)
library(ggrepel)
library(patchwork) # For combining plots

# Data preprocessing
plot_data <- results %>%
  # Calculate absolute Rho values
  mutate(
    abs_rho = abs(Rho),
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
    class_median_abs_rho = median(abs_rho, na.rm = TRUE),
    class_pct_significant = mean(Adjusted_P_value < 0.05, na.rm = TRUE) * 100,
    # Rank within each class by absolute correlation
    class_rank = rank(-abs_rho)
  ) %>%
  ungroup() %>%
  # Filter classes with at least 10 members (adjust this threshold as needed)
  filter(class_size >= 10) %>%
  # Sort by median absolute Rho within class
  arrange(HMDB.Class, desc(abs_rho)) %>%
  # Add significance markers and labels
  mutate(
    # Flag for direction of correlation
    direction = ifelse(Rho > 0, "Positive", "Negative"),
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
    median_abs_rho = median(abs_rho, na.rm = TRUE),
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
p1 <-ggplot(plot_data, aes(x = metabolite_sort, y = abs_rho)) +
  # Add alternating backgrounds for classes
  geom_rect(data = class_summary,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = "gray95", alpha = 0.5, inherit.aes = FALSE) +
  # Add points with site colors and direction shapes
  geom_point(aes(color = Site, 
                 shape = direction,
                 size = p_size), 
             alpha = 0.85) +
  # Add FDR significance threshold
  geom_hline(yintercept = median(subset(plot_data, Adjusted_P_value == 0.05)$abs_rho, na.rm = TRUE), 
             linetype = "dashed", color = "black", size = 0.5) +
  # Add labels for top metabolites
  geom_text_repel(
    data = subset(plot_data, label != ""),
    aes(label = label),
    size = 5,
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
  scale_shape_manual(values = c("Positive" = 16, "Negative" = 25),
                     name = "Correlation") +
  scale_size_continuous(
    name = "P-value", 
    breaks = c(1, 2, 3, 4),
    labels = c("ns", "P < 0.05", "P < 0.01", "P < 0.001"),
    range = c(2, 5)
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.1),
    minor_breaks = seq(0, 1, by = 0.05),
    limits = c(0, max(plot_data$abs_rho) * 1.05),
    expand = expansion(mult = c(0.02, 0.1))
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
  labs(x="HMDB.Class", y = "Absolute Correlation (|Rho|)"
  )+ylim(c(0.2,0.6))






# 绘制是否为微生物来源代谢物的Rho大小

plot_microbial_source_data<-plot_data

plot_microbial_source_data$HMDB.Source.Microbial<-gsub("NA","FALSE",plot_microbial_source_data$HMDB.Source.Microbial)

HMDB.Source.Microbial_color<-c("","")

ggplot(data=plot_microbial_source_data,aes(x=HMDB.Source.Microbial,y=abs_rho,fill=HMDB.Source.Microbial))+ #”fill=“设置填充颜色
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
    quantiles = 3
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

# 统计每个类别的数量和百分比
class_summary <- data %>%
  dpcount(HMDB.Class) %>%
  mutate(
    total = sum(n),
    percentage = n / total * 100
  ) %>%
  arrange(desc(percentage))

# 为了使图表更易读，只显示前15个最常见的类别，其余归为"其他"类别
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

# 创建一个带有单一堆叠柱的数据框
stacked_data <- class_summary %>%
  mutate(
    x = 1,  # 所有类别共用一个x轴位置，形成堆叠
    ymin = lag(cumsum(percentage), default = 0),
    ymax = cumsum(percentage),
    pos = (ymin + ymax) / 2  # 标签位置
  )

# 生成足够的颜色
num_colors <- nrow(stacked_data)+1
if(num_colors <= 8) {
  colors <- brewer.pal(max(3, num_colors), "Set1")
} else {
  colors <- colorRampPalette(brewer.pal(8, "Set1"))(num_colors)
}

# 创建堆叠柱状图
p <- ggplot(stacked_data) +
  geom_rect(aes(
    xmin = x - 0.4, 
    xmax = x + 0.4, 
    ymin = ymin, 
    ymax = ymax,
    fill = HMDB.Class
  )) +
  scale_fill_manual(values = colors) +
  geom_text(
    aes(
      x = x, 
      y = pos, 
      label = sprintf("%s\n(%.1f%%)", HMDB.Class, percentage)
    ),
    size = 3
  ) +
  theme_minimal() +
  labs(
    title = "metabolite_annotation",
    y = "Persent (%)",
    fill = "HMDB.Class"
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  ) +
  scale_y_continuous(breaks = seq(0, 100, 10))




# 读取数据
# 假设你的数据在名为"metabolite_annotation.csv"的文件中
data <- plot_data



# 确保Site列中的值符合预期
expected_sites <- c("gut", "oral", "skin", "nasal")
data$Site <- tolower(data$Site)  # 转换为小写
data <- data %>% filter(Site %in% expected_sites)  # 只保留四个指定位点的数据

# 假设数据中已经包含了正确处理的类别信息（包括已经归纳的"others"）
# 对于每个位点，统计各个类别的数量和百分比
site_class_final <- data %>%
  group_by(Site, HMDB.Class) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  group_by(Site) %>%
  mutate(
    total = sum(count),
    percentage = count / total * 100
  ) %>%
  arrange(Site, desc(percentage)) %>%
  # 确保"others"类别（如果存在）排在最后
  mutate(HMDB.Class = fct_relevel(as.factor(HMDB.Class), "others", after = Inf))

# 生成足够的颜色
unique_classes <- unique(as.character(site_class_final$HMDB.Class))
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
colors<-colors[-6]
# 确保"others"类别（如果存在）使用灰色
if("others" %in% unique_classes) {
  names(colors) <- unique_classes
  colors["others"] <- "#999999" # 灰色
}

# 创建堆叠柱状图
p <- ggplot(site_class_final, aes(x = Site, y = percentage, fill = HMDB.Class)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(
    title = "各位点(Site)代谢物HMDB类别分布",
    x = "位点",
    y = "百分比 (%)",
    fill = "代谢物类别"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right"
  ) +
  scale_y_continuous(breaks = seq(0, 100, 10))


