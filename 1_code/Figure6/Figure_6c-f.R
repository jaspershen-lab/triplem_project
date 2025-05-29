rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")

library(tidyverse)
library(tidymass)

###load("data)
load("3_data_analysis/gut_microbiome/data_preparation/object_cross_section")

gut_object <- object_cross_section

load(
  "3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section"
)

metabolomics_object <- object_cross_section

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


load("3_data_analysis/oral_microbiome/data_preparation/object_cross_section")

oral_object <- object_cross_section

load(
  "3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section"
)

metabolite_annotation <- read_excel(
  "3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx"
)
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


##
gut_data <- gut_temp_object@expression_data
oral_data <- oral_temp_object@expression_data

metabolome_data <- metabolomics_temp_object@expression_data

group_comparison_results<-readRDS("1_code/Figure6/group_comparison_results")

sample_info$IRIS[21]<-"IR"
sample_info$IRIS[24]<-"IR"


# 可视化分组比较结果
# plot_group_comparison_results <- function(results) {
#   library(ggplot2)
#   library(dplyr)
#   library(patchwork)
#   library(reshape2)
#   
#   summary_df <- results$summary
#   group1_name <- results$group1_name
#   group2_name <- results$group2_name
#   
#   # 1. R²值比较散点图
#   p1 <- ggplot(summary_df, aes(x = group1_r2, y = group2_r2)) +
#     geom_point(aes(color = p_diff_adjusted < 0.05), alpha = 0.7, size = 2) +
#     geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
#     scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"),
#                        labels = c("Not significant", "Significant (p < 0.05)")) +
#     theme_minimal() +
#     labs(title = paste("R² Comparison:", group1_name, "vs", group2_name),
#          x = paste("R² in", group1_name, "group"),
#          y = paste("R² in", group2_name, "group"),
#          color = "Group Difference") +
#     coord_equal()
#   
#   # 2. R²差异的分布图
#   p2 <- ggplot(summary_df, aes(x = r2_difference)) +
#     geom_histogram(bins = 30, alpha = 0.7, fill = "steelblue") +
#     geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
#     theme_minimal() +
#     labs(title = "Distribution of R² Differences",
#          x = paste("R² Difference (", group1_name, "-", group2_name, ")"),
#          y = "Count")
#   
#   # 3. 交互特征数量比较
#   interaction_data <- summary_df %>%
#     select(metabolite, group1_n_interaction, group2_n_interaction) %>%
#     melt(id.vars = "metabolite", 
#          variable.name = "group", 
#          value.name = "n_interactions") %>%
#     mutate(group = ifelse(group == "group1_n_interaction", group1_name, group2_name))
#   
#   p3 <- ggplot(interaction_data, aes(x = group, y = n_interactions, fill = group)) +
#     geom_boxplot(alpha = 0.7) +
#     geom_jitter(width = 0.2, alpha = 0.5) +
#     scale_fill_manual(values=c("IR"="#E69F00", "IS"="#0072B2")) +
#     theme_minimal() +
#     labs(title = "Number of Selected Interaction Features",
#          x = "Group",
#          y = "Number of Interaction Features",
#          fill = "Group")+stat_compare_means(label = "p.format")
#   
#   # 4. 效应大小vs显著性
#   p4 <- ggplot(summary_df, aes(x = effect_size, y = -log10(p_diff_adjusted))) +
#     geom_point(aes(color = abs(effect_size) > 0.5 & p_diff_adjusted < 0.05), 
#                alpha = 0.7, size = 2) +
#     geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
#     geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray") +
#     scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"),
#                        labels = c("Not significant", "Significant & Large effect")) +
#     theme_minimal() +
#     labs(title = "Effect Size vs Significance",
#          x = "Effect Size (Cohen's d)",
#          y = "-log10(Adjusted p-value)",
#          color = "Classification")
#   
#   # 5. 火山图样式的R²差异图
#   p5 <- ggplot(summary_df, aes(x = r2_difference, y = -log10(p_diff_adjusted))) +
#     geom_point(aes(color = p_diff_adjusted < 0.05 & abs(r2_difference) > 0.1), 
#                alpha = 0.7, size = 2) +
#     geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
#     geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "gray") +
#     scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"),
#                        labels = c("Not significant", "Significant difference")) +
#     theme_minimal() +
#     labs(title = "Volcano Plot: R² Differences",
#          x = paste("R² Difference (", group1_name, "-", group2_name, ")"),
#          y = "-log10(Adjusted p-value)",
#          color = "Significance")
#   
#   # 6. 相关性热图：交互特征数量与R²表现
#   correlation_data <- summary_df %>%
#     select(group1_r2, group2_r2, group1_n_interaction, group2_n_interaction, 
#            r2_difference, interaction_difference) %>%
#     cor(use = "complete.obs")
#   
#   correlation_melted <- melt(correlation_data)
#   
#   p6 <- ggplot(correlation_melted, aes(x = Var1, y = Var2, fill = value)) +
#     geom_tile() +
#     geom_text(aes(label = round(value, 2)), color = "white", size = 3) +
#     scale_fill_gradient2(low = "blue", mid = "white", high = "red",
#                          midpoint = 0, limit = c(-1, 1)) +
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     labs(title = "Correlation Matrix",
#          x = "", y = "", fill = "Correlation")
#   
#   # 组合所有图表
#   combined_plots <- (p1 + p2) / (p3 + p4) / (p5 + p6) +
#     plot_layout(heights = c(1, 1, 1))
#   
#   return(list(
#     r2_comparison = p1,
#     r2_difference_dist = p2,
#     interaction_features = p3,
#     effect_size_significance = p4,
#     volcano_plot = p5,
#     correlation_heatmap = p6,
#     combined = combined_plots
#   ))
# }

# 生成详细的分组比较报告
generate_group_comparison_report <- function(results, top_n = 10) {
  summary_df <- results$summary
  group1_name <- results$group1_name
  group2_name <- results$group2_name
  
  # 总体统计
  cat("=== GROUP COMPARISON ANALYSIS REPORT ===\n\n")
  cat(sprintf("Group 1 (%s): %d samples\n", group1_name, length(results$group_info$group1_samples)))
  cat(sprintf("Group 2 (%s): %d samples\n", group2_name, length(results$group_info$group2_samples)))
  cat(sprintf("Total metabolites analyzed: %d\n\n", nrow(summary_df)))
  
  # 显著差异统计
  significant_diffs <- sum(summary_df$p_diff_adjusted < 0.05, na.rm = TRUE)
  cat(sprintf("Metabolites with significant group differences (p < 0.05): %d (%.1f%%)\n", 
              significant_diffs, 100 * significant_diffs / nrow(summary_df)))
  
  large_effect <- sum(abs(summary_df$effect_size) > 0.5 & summary_df$p_diff_adjusted < 0.05, na.rm = TRUE)
  cat(sprintf("Metabolites with large effect size (|d| > 0.5) and significant: %d (%.1f%%)\n\n", 
              large_effect, 100 * large_effect / nrow(summary_df)))
  
  # R²表现比较
  cat("=== R² PERFORMANCE COMPARISON ===\n")
  cat(sprintf("Mean R² in %s: %.3f ± %.3f\n", 
              group1_name, mean(summary_df$group1_r2, na.rm = TRUE), sd(summary_df$group1_r2, na.rm = TRUE)))
  cat(sprintf("Mean R² in %s: %.3f ± %.3f\n", 
              group2_name, mean(summary_df$group2_r2, na.rm = TRUE), sd(summary_df$group2_r2, na.rm = TRUE)))
  cat(sprintf("Mean R² difference: %.3f ± %.3f\n\n", 
              mean(summary_df$r2_difference, na.rm = TRUE), sd(summary_df$r2_difference, na.rm = TRUE)))
  
  # 交互特征使用情况
  cat("=== INTERACTION FEATURES USAGE ===\n")
  cat(sprintf("Mean interaction features in %s: %.1f ± %.1f\n", 
              group1_name, mean(summary_df$group1_n_interaction, na.rm = TRUE), 
              sd(summary_df$group1_n_interaction, na.rm = TRUE)))
  cat(sprintf("Mean interaction features in %s: %.1f ± %.1f\n", 
              group2_name, mean(summary_df$group2_n_interaction, na.rm = TRUE), 
              sd(summary_df$group2_n_interaction, na.rm = TRUE)))
  cat(sprintf("Mean interaction difference: %.1f ± %.1f\n\n", 
              mean(summary_df$interaction_difference, na.rm = TRUE), 
              sd(summary_df$interaction_difference, na.rm = TRUE)))
  
  # Top差异代谢物
  cat(sprintf("=== TOP %d METABOLITES WITH LARGEST GROUP DIFFERENCES ===\n", top_n))
  top_metabolites <- summary_df %>%
    arrange(p_diff_adjusted, desc(abs(r2_difference))) %>%
    head(top_n)
  
  for(i in 1:nrow(top_metabolites)) {
    row <- top_metabolites[i, ]
    cat(sprintf("%d. %s\n", i, row$metabolite))
    cat(sprintf("   %s R²: %.3f (95%% CI: %.3f-%.3f)\n", 
                group1_name, row$group1_r2, row$group1_ci_lower, row$group1_ci_upper))
    cat(sprintf("   %s R²: %.3f (95%% CI: %.3f-%.3f)\n", 
                group2_name, row$group2_r2, row$group2_ci_lower, row$group2_ci_upper))
    cat(sprintf("   Difference: %.3f (Effect size: %.2f, p = %.2e)\n", 
                row$r2_difference, row$effect_size, row$p_diff_adjusted))
    cat(sprintf("   Interaction features: %s=%d, %s=%d (diff=%d)\n\n", 
                group1_name, row$group1_n_interaction, 
                group2_name, row$group2_n_interaction, 
                row$interaction_difference))
  }
  
  # 返回筛选后的重要结果
  return(list(
    significant_metabolites = summary_df[summary_df$p_diff_adjusted < 0.05, ],
    large_effect_metabolites = summary_df[abs(summary_df$effect_size) > 0.5 & 
                                            summary_df$p_diff_adjusted < 0.05, ],
    summary_stats = list(
      n_significant = significant_diffs,
      n_large_effect = large_effect,
      mean_r2_group1 = mean(summary_df$group1_r2, na.rm = TRUE),
      mean_r2_group2 = mean(summary_df$group2_r2, na.rm = TRUE),
      mean_interaction_group1 = mean(summary_df$group1_n_interaction, na.rm = TRUE),
      mean_interaction_group2 = mean(summary_df$group2_n_interaction, na.rm = TRUE)
    )
  ))
}





# 生成可视化结果
visualization_results <- 
  plot_group_comparison_results(group_comparison_results)

# 显示组合图表
print(visualization_results$combined)

visualization_results$r2_comparison
visualization_results$combined

# 生成详细报告
detailed_report <- generate_group_comparison_report(group_comparison_results, top_n = 15)

# 查看显著差异的代谢物
significant_metabolites <- detailed_report$significant_metabolites
print(head(significant_metabolites, 10))


significant_metabolites<-
  merge(significant_metabolites,metabolite_annotation,by.x="metabolite",by.y="variable_id")


#挑选需要展示的代谢物

Figure_6c <-
  visualization_results$interaction_features 

Figure_6d<-
  visualization_results$volcano_plot +
  theme_bw() 

Figure_6c
Figure_6d

metabolites_index<-c("HMDB0001870","HMDB0000637","HMDB0001856","HMDB0059655","HMDB0002466","HMDB0000707","HMDB0000258","HMDB0002820","HMDB0000078","HMDB0002320","HMDB0001859","HMDB0001870","HMDB0000159","HMDB0000517","HMDB0000715","HMDB0000158","HMDB0002466","HMDB0000707","HMDB0001190","HMDB0000039","HMDB0002212")


significant_metabolites<-subset(significant_metabolites,HMDB%in%metabolites_index)



# plot_metabolites_mirror_style <- function(significant_metabolites, top_n = 30) {
#   library(ggplot2)
#   library(dplyr)
#   
#   # 准备数据
#   if(nrow(significant_metabolites) > top_n) {
#     plot_data <- significant_metabolites %>%
#       arrange(desc(abs(r2_difference))) %>%
#       head(top_n)
#   } else {
#     plot_data <- significant_metabolites %>%
#       arrange(desc(abs(r2_difference)))
#   }
#   
#   # 按原图风格排序：先正值（降序），后负值（从小到大，即从最负开始）
#   positive_data <- plot_data[plot_data$r2_difference > 0, ] %>%
#     arrange(desc(r2_difference))
#   negative_data <- plot_data[plot_data$r2_difference < 0, ] %>%
#     arrange(desc(r2_difference))  # 负值从最小（最负）到最大（接近0）
#   
#   # 重新组合数据
#   plot_data_ordered <- rbind(positive_data, negative_data)
#   
#   # 创建x轴位置
#   plot_data_ordered$x_pos <- 1:nrow(plot_data_ordered)
#   
#   # 为双向显示创建数据
#   plot_data_ordered$y_upper <- ifelse(plot_data_ordered$r2_difference > 0, 
#                                       plot_data_ordered$r2_difference, 0)
#   plot_data_ordered$y_lower <- ifelse(plot_data_ordered$r2_difference < 0, 
#                                       -plot_data_ordered$r2_difference, 0)  # 取绝对值显示在下方
#   
#   # 找到最大值用于设置y轴
#   max_val <- max(abs(plot_data_ordered$r2_difference))
#   
#   p <- ggplot(plot_data_ordered, aes(x = x_pos)) +
#     # 上方条形（IR优势，绿色）
#     geom_col(aes(y = y_upper), fill = "#E69F00", alpha = 0.9, width = 0.8) +
#     # 下方条形（IS优势，红色，向下显示）
#     geom_col(aes(y = -y_lower), fill = "#0072B2", alpha = 0.9, width = 0.8) +
#     # 零线
#     geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
#     # 设置y轴范围
#     scale_y_continuous(
#       limits = c(-max_val * 1.1, max_val * 1.1),
#       breaks = seq(-max_val, max_val, length.out = 7),
#       labels = function(x) sprintf("%.1f", abs(x))  # 显示绝对值
#     ) +
#     # x轴设置
#     scale_x_continuous(
#       breaks = plot_data_ordered$x_pos,
#       labels = plot_data_ordered$HMDB.Name,
#       expand = c(0.01, 0.01)
#     ) +
#     # 主题
#     theme_minimal() +
#     theme(
#       axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
#       axis.text.y = element_text(size = 10),
#       axis.title.y = element_text(size = 12, face = "bold"),
#       plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
#       panel.grid.major.x = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.grid.major.y = element_line(color = "gray90", size = 0.3),
#       axis.line.x = element_line(color = "black", size = 0.5),
#       plot.margin = margin(t = 20, r = 10, b = 10, l = 10)
#     ) +
#     # 标签
#     labs(
#       title = "Metabolite R² Differences: IR vs IS Groups",
#       x = "",
#       y = "R² Difference",
#       caption = ""
#     ) +
#     # 添加颜色图例标注
#     annotate("text", x = length(plot_data_ordered$metabolite) * 0.05, 
#              y = max_val * 1, label = "IR", 
#              color = "#E69F00", size = 4, fontface = "bold") +
#     annotate("text", x = length(plot_data_ordered$metabolite) * 0.05, 
#              y = -max_val * 0.9, label = "IS", 
#              color = "#0072B2", size = 4, fontface = "bold")
#   
#   return(p)
# }

Figure_6e <- plot_metabolites_mirror_style(significant_metabolites, top_n = 25)

Figure_6e

library(metpath)
library(tidyverse)
data("hmdb_pathway", package = "metpath")
data("query_id_hmdb", package = "metpath")

#get the class of pathways
pathway_class = 
  hmdb_pathway@pathway_class %>% 
  unlist()

remain_idx = which(unlist(pathway_class) == "Metabolic;primary_pathway")


hmdb_pathway = 
  hmdb_pathway[remain_idx]

compound_list<-hmdb_pathway@compound_list


# 需要筛选的 HMDB.ID 列表

# 生成详细报告
detailed_report <- generate_group_comparison_report(group_comparison_results, top_n = 15)

# 查看显著差异的代谢物
significant_metabolites <- detailed_report$significant_metabolites
print(head(significant_metabolites, 10))

significant_metabolites<-merge(significant_metabolites,metabolite_annotation,by.x="metabolite",by.y="variable_id")

selected_hmdb_ids <- metabolite_annotation$HMDB.ID

# 筛选 compound_list 中每个 data.frame 的元素
filtered_compound_list <- lapply(compound_list, function(df) {
  df[df$HMDB.ID %in% selected_hmdb_ids, ]
})

# 查看筛选结果
filtered_compound_list

hmdb_pathway@compound_list<-filtered_compound_list

result = 
  enrich_hmdb(query_id = significant_metabolites$HMDB, 
              query_type = "compound", 
              id_type = "HMDB",
              pathway_database = hmdb_pathway,
              only_primary_pathway = TRUE,
              p_cutoff = 0.05, 
              p_adjust_method = "none", 
              threads = 3)

enrich_results<-result@result
enrich_results$p_value<-enrich_results$p_value/enrich_results$mapped_number

enrich_results<-subset(enrich_results,p_value<0.05)

# 加载必需的包
library(ggplot2)
library(dplyr)

# 创建示例数据
data <- enrich_results

# 计算-log10(p_value)用于颜色映射
data$neg_log10_p <- -log10(data$p_value)

# 按mapped_percentage排序
data <- data %>%
  arrange(desc(mapped_percentage)) %>%
  slice_head(n = 12)

# 创建柱状图
Figure_6f<- ggplot(data, aes(x = reorder(pathway_name, mapped_percentage), 
                            y = mapped_percentage, 
                            fill = neg_log10_p)) +
  
  # 添加柱状图
  geom_col(alpha = 0.8, width = 0.7) +
  
  # 设置颜色渐变
  scale_fill_gradient(low = "lightblue", high = "darkred", 
                      name = "-log10(P-value)") +
  
  # 翻转坐标轴，让条目名称在Y轴上
  coord_flip() +
  
  # 添加标题和轴标签
  labs(title = "",
       x = "Pathway",
       y = "Mapped_percentage (%)") +
  
  # 自定义主题
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )


print(Figure_6c)
print(Figure_6d)
print(Figure_6e)
print(Figure_6f)

ggsave(
  Figure_6c,
  filename = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_6/figure_6c.pdf"
  ),
  width = 6,
  height = 5
)

ggsave(
  Figure_6d,
  filename = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_6/figure_6d.pdf"
  ),
  width = 8,
  height = 5
)

ggsave(
  Figure_6e,
  filename = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_6/figure_6e.pdf"
  ),
  width = 6,
  height = 5
)

ggsave(
  Figure_6f,
  filename = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_6/figure_6f.pdf"
  ),
  width = 6,
  height = 5
)

