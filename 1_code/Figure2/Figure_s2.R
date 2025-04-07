rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")
setwd("1_code/4_site_merge/")
library(tidymass)
library(tidyverse)
library(readxl)
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
metabolite_annotation<-subset(metabolite_annotation,!(HMDB.Name=="NA"))
metabolite_annotation_micro<-subset(metabolite_annotation,HMDB.Source.Microbial=="TRUE")
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



##
##adjust BMI, sex, and IRIS, ethnicity
library(tidyverse)
library(ggpubr)
library(rstatix)

gut_expression_data <-
  extract_expression_data(gut_object) %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

library(plyr)

gut_sample_info <-
  gut_object@sample_info

#######adjust BMI, sex, and IRIS, ethnicity
gut_expression_data <-
  lm_adjust(expression_data = gut_expression_data,
            sample_info = gut_sample_info,
            threads = 3)

gut_temp_object <- gut_object
gut_temp_object@expression_data <- gut_expression_data


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



##
##adjust BMI, sex, and IRIS, ethnicity
library(tidyverse)
library(ggpubr)
library(rstatix)

skin_expression_data <-
  extract_expression_data(skin_object) %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

library(plyr)

skin_sample_info <-
  skin_object@sample_info

#######adjust BMI, sex, and IRIS, ethnicity
skin_expression_data <-
  lm_adjust(expression_data = skin_expression_data,
            sample_info = skin_sample_info,
            threads = 3)

skin_temp_object <- skin_object
skin_temp_object@expression_data <- skin_expression_data



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



##
##adjust BMI, sex, and IRIS, ethnicity
library(tidyverse)
library(ggpubr)
library(rstatix)

nasal_expression_data <-
  extract_expression_data(nasal_object) %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

library(plyr)

nasal_sample_info <-
  nasal_object@sample_info

#######adjust BMI, sex, and IRIS, ethnicity
nasal_expression_data <-
  lm_adjust(expression_data = nasal_expression_data,
            sample_info = nasal_sample_info,
            threads = 3)

nasal_temp_object <- nasal_object
nasal_temp_object@expression_data <- nasal_expression_data

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

# 加载必要的包
library(Hmisc)      # 用于计算相关性
library(igraph)     # 用于构建网络
library(ggplot2)    # 用于绘图
library(reshape2)   # 用于数据重构
library(tidyverse)  # 数据处理

# 假设您的数据存储在以下数据框中:
# gut_data: gut菌群数据
# oral_data: oral菌群数据
# skin_data: skin菌群数据
# nasal_data: nasal菌群数据
# metabolites: 代谢物数据

# 函数：计算相关性并筛选显著结果
calculate_correlations <- function(microbiome_data,taxdata, metabolite_data,metabolite_anno, site, threshold = 0.3, p_value = 0.05) {
  # 确保样本行匹配
  common_samples <- intersect(rownames(microbiome_data), rownames(metabolite_data))
  
  # 子集化数据
  microbiome_subset <- microbiome_data[common_samples,]
  
  colnames(microbiome_subset)<-  paste(site,taxdata$Genus,sep = "_")
  
  metabolite_subset <- metabolite_data[common_samples,metabolite_anno$variable_id]
  
  
  colnames(metabolite_subset)<-  metabolite_anno$HMDB.Name
  # 计算相关性矩阵和p值
  cors <- matrix(NA, nrow = ncol(microbiome_subset), ncol = ncol(metabolite_subset))
  pvals <- matrix(NA, nrow = ncol(microbiome_subset), ncol = ncol(metabolite_subset))
  
  # 逐个计算相关性
  for(i in 1:ncol(microbiome_subset)) {
    for(j in 1:ncol(metabolite_subset)) {
      test_result <- cor.test(microbiome_subset[,i], metabolite_subset[,j], 
                              method = "spearman", exact = FALSE)
      cors[i,j] <- test_result$estimate
      pvals[i,j] <- test_result$p.value
    }
  }
  
  # 设置行列名
  rownames(cors) <- colnames(microbiome_subset)
  colnames(cors) <- colnames(metabolite_subset)
  rownames(pvals) <- colnames(microbiome_subset)
  colnames(pvals) <- colnames(metabolite_subset)
  
  # 转换为长数据格式
  cors_df <- melt(cors)
  pvals_df <- melt(pvals)
  
  # 合并相关系数和p值
  result_df <- data.frame(
    Microbiome = cors_df$Var1,
    Metabolite = cors_df$Var2,
    Correlation = cors_df$value,
    P_value = pvals_df$value,
    Site = site
  )
  
  # 筛选显著相关的结果
  significant_cors <- result_df %>%
    filter(abs(Correlation) >= threshold & P_value <= p_value)
  
  return(significant_cors)
}

# 对每个部位计算相关性
gut_cors <- calculate_correlations(t(gut_temp_object@expression_data),gut_temp_object@variable_info,t(metabolomics_temp_object@expression_data),metabolite_annotation, "Gut")
oral_cors <- calculate_correlations(t(oral_temp_object@expression_data),oral_temp_object@variable_info, t(metabolomics_temp_object@expression_data),metabolite_annotation ,"Oral")
skin_cors <- calculate_correlations(t(skin_temp_object@expression_data),skin_temp_object@variable_info ,t(metabolomics_temp_object@expression_data),metabolite_annotation, "Skin")
nasal_cors <- calculate_correlations(t(nasal_temp_object@expression_data),nasal_temp_object@variable_info, t(metabolomics_temp_object@expression_data),metabolite_annotation, "Nasal")

# 合并所有相关性结果
all_cors <- rbind(gut_cors, oral_cors, skin_cors, nasal_cors)



correlation_data<-all_cors


# 使用ggraph创建网络图
edges <- correlation_data %>%
  select(Microbiome, Metabolite, Correlation, Site)

# 创建节点列表
nodes <- data.frame(
  name = unique(c(correlation_data$Microbiome, correlation_data$Metabolite)),
  type = ifelse(unique(c(correlation_data$Microbiome, correlation_data$Metabolite)) %in% correlation_data$Microbiome, "Microbiome", "Metabolite"),
  site = ifelse(unique(c(correlation_data$Microbiome, correlation_data$Metabolite)) %in% correlation_data$Microbiome, 
                correlation_data$Site[match(unique(c(correlation_data$Microbiome, correlation_data$Metabolite)), correlation_data$Microbiome)], 
                "Metabolite")
)

# 创建igraph对象
graph <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
# 计算每个节点的degree
# 计算degree
node_degrees <- degree(graph)
nodes$degree <- node_degrees[match(nodes$name, names(node_degrees))]

# 筛选top 20节点
top_nodes <- nodes %>%
  group_by(site) %>%
  top_n(10, degree) %>%
  ungroup()

# 筛选包含top节点的边
filtered_edges <- edges %>%
  filter(Microbiome %in% top_nodes$name | Metabolite %in% top_nodes$name)

# 创建新的节点列表（包含所有在筛选后的边中出现的节点）
filtered_nodes <- nodes %>%
  filter(name %in% unique(c(filtered_edges$Microbiome, filtered_edges$Metabolite)))


# 获取每个site的前10个节点的名字
top_10_nodes <- filtered_nodes %>%
  group_by(site) %>%
  top_n(10, degree) %>%
  ungroup() %>%
  pull(name)

# 添加是否显示标签的列
filtered_nodes$show_label <- filtered_nodes$name %in% top_10_nodes

# 创建筛选后的igraph对象
filtered_graph <- graph_from_data_frame(d = filtered_edges, vertices = filtered_nodes, directed = FALSE)

# 重新计算degree
filtered_degrees <- degree(filtered_graph)
V(filtered_graph)$degree <- filtered_degrees
V(filtered_graph)$show_label <- filtered_nodes$show_label
# 创建ggraph可视化

library(ggraph)
p<-ggraph(filtered_graph, layout = "stress") +
  # 添加边
  geom_edge_link(aes(edge_alpha = abs(Correlation),
                     edge_width = abs(Correlation),
                     color = Correlation > 0),
                 show.legend = TRUE) +
  # 添加节点
  geom_node_point(aes(fill = site, 
                      size = degree,
                      shape = type),color="white") +
  # 添加节点标签
  geom_node_text(aes(label = ifelse(show_label, name, "")), 
                 repel = TRUE, 
                 size = 2,
                 max.overlaps = 20)+
  # 设置配色
  scale_edge_color_manual(values = c("TRUE" = "#FF9999", "FALSE" = "#9999FF"),
                          name = "Correlation",
                          labels = c("Negative", "Positive")) +
  scale_fill_manual(values = c("Gut" = "#edd064",
                               "Oral" = "#a1d5b9",
                               "Skin" = "#f2ccac",
                               "Nasal" = "#a17db4",
                               "Metabolite" = "#FF9999")) +
  # 设置节点大小
  scale_size_continuous(range = c(2, 10), name = "Degree")+
  # 设置节点形状
  scale_shape_manual(values = c("Microbiome" = 21, "Metabolite" = 22)) +
  # 设置边的透明度和宽度
  scale_edge_alpha(range = c(0.2, 0.8)) +
  scale_edge_width(range = c(0.3, 2)) +
  # 主题设置
  theme_graph() +
  theme(legend.position = "none") +
  guides(
    color = guide_legend(title = "Type"),
    size = guide_legend(title = "Node Type"),
    shape = guide_legend(title = "Node Type"),
    edge_alpha = guide_legend(title = "Correlation Strength"),
    edge_width = guide_legend(title = "Correlation Strength")
  ) +
  labs(title = "")

setwd(r4projects::get_project_wd())
ggsave(p,
       filename = "4_manuscript/Figures/Figure_2/figure_s2_network.pdf",
       width = 10,
       height = 10)





##########


# 导入必要的库
library(igraph)
library(dplyr)

# 基于原始网络分析重要节点
# 假设filtered_graph和filtered_nodes已经存在

# 1. 提取所有节点的重要性度量
node_importance <- data.frame(
  name = V(filtered_graph)$name,
  site = V(filtered_graph)$site,
  type = V(filtered_graph)$type,
  degree = V(filtered_graph)$degree,
  betweenness = betweenness(filtered_graph),
  closeness = closeness(filtered_graph),
  eigen_centrality = eigen_centrality(filtered_graph)$vector
)

# 2. 分别统计微生物和代谢物的重要节点
# 微生物节点(按部位和度数排序)
important_microbes <- node_importance %>%
  filter(type == "Microbiome") %>%
  group_by(site) %>%
  arrange(desc(degree)) %>%
  slice_head(n = 10) %>%
  ungroup()

# 代谢物节点(按度数排序)
important_metabolites <- node_importance %>%
  filter(type == "Metabolite") %>%
  arrange(desc(degree)) %>%
  slice_head(n = 20)



# 4. 统计每个部位对代谢物有最强影响的微生物
# 这里我们根据与代谢物的相关强度来评估影响力

# 首先，计算每个微生物与代谢物的平均相关性强度
microbe_influence <- all_cors %>%
  dplyr::group_by(Microbiome, Site) %>%
  dplyr::summarize(
    avg_correlation = mean(abs(Correlation)),
    total_connections = n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Site, desc(avg_correlation))

# 获取每个部位影响力最强的微生物(前5名)
top_influential_microbes <- microbe_influence %>%
  dplyr::group_by(Site) %>%
  dplyr::slice_head(n = 5) %>%
  dplyr::ungroup()

cat("\n每个部位对代谢物影响最强的微生物(前5名):\n")
print(top_influential_microbes, n = nrow(top_influential_microbes))

# 5. 找出被多个部位共同影响的代谢物
metabolites_by_site <- all_cors %>%
  dplyr::group_by(Metabolite) %>%
  dplyr::summarize(
    sites = list(unique(Site)),
    site_count = length(unique(Site)),
    total_microbes = length(unique(Microbiome)),
    avg_correlation = mean(abs(Correlation)),
    .groups = "drop"
  ) %>%
  dplyr:: arrange(desc(site_count), desc(avg_correlation))

# 获取被多个部位共同影响的代谢物(出现在2个以上部位)
shared_metabolites <- metabolites_by_site %>%
  dplyr::filter(site_count > 1) %>%
  dplyr::slice_head(n = 15)

cat("\n被多个部位微生物共同影响的代谢物(前15名):\n")
print(shared_metabolites, n = nrow(shared_metabolites))

# 6. 创建关于重要节点的可视化总结表格
library(ggplot2)
library(forcats)

# 可视化每个部位的顶级微生物
p1 <- important_microbes %>%
  group_by(site) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(name = fct_reorder(name, degree)) %>%
  ggplot(aes(x = degree, y = name, fill = site)) +
  geom_col() +
  scale_fill_manual(values = c("Gut" = "#edd064",
                               "Oral" = "#a1d5b9",
                               "Skin" = "#f2ccac",
                               "Nasal" = "#a17db4")) +
  facet_wrap(~site, scales = "free_y",nrow = 1) +
  labs(title = "",
       x = "Degree",
       y = "") +
  theme_bw()

# 保存结果
ggsave("4_manuscript/Figures/Figure_2/Figure_s2_important_microbes.pdf", 
       p1, width = 15, height = 4)

# 可视化重要的代谢物
p2 <- important_metabolites %>%
  slice_head(n = 10) %>%
  mutate(name = fct_reorder(name, degree)) %>%
  ggplot(aes(x = degree, y = name, fill = "Metabolite")) +
  geom_col(fill = "#FF9999") +
  labs(title = "",
       x = "Degree",
       y = "") +
  theme_bw()

ggsave("4_manuscript/Figures/Figure_2/Figure_s2_important_Metabolite.pdf", 
       p2, width = 4, height = 4)