lm_adjust <-
  function(expression_data,
           sample_info,
           threads = 5) {
    library(future)
    library(furrr)
    # plan(strategy = multisession(workers = threads))
    new_expression_data <-
      rownames(expression_data) %>%
      purrr::map(function(name) {
        # cat(name, " ")
        x = as.numeric(expression_data[name,])
        temp_data =
          data.frame(x = x,
                     sample_info)
        
        
        temp_data$Gender[temp_data$Gender == 'Female'] = 0
        temp_data$Gender[temp_data$Gender == 'Male'] = 1
        temp_data$Gender = as.numeric(temp_data$Gender)
        
        temp_data$IRIS[temp_data$IRIS == 'IR'] = 1
        temp_data$IRIS[temp_data$IRIS == 'IS'] = 2
        temp_data$IRIS[temp_data$IRIS == 'Unknown'] = 0
        temp_data$IRIS = as.numeric(temp_data$IRIS)
        
        adjusted_x <-
          lm(x ~ Gender + BMI + IRIS, data = temp_data) %>%
          residuals()
        adjusted_x
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    colnames(new_expression_data) <-
      colnames(expression_data)
    
    rownames(new_expression_data) <-
      rownames(expression_data)
    new_expression_data
  }





body_site_color = c(
  "gut" = "#edd064",
  "skin" = "#f2ccac",
  "oral" = "#a1d5b9",
  "nasal" = "#a17db4"
)


iris_color = c(
  "IR" = "#bf3f7c",
  "IS" = "#46c1be",
  "Unknown" = "#546672"
)


sex_color <-
  c(
    "Female" = wesanderson::wes_palettes$Rushmore1[3],
    "Male" = wesanderson::wes_palettes$Rushmore1[1]
  )

ethnicity_color <-
  c(
    "Caucasian" = wesanderson::wes_palettes$Darjeeling2[1],
    "Asian" = wesanderson::wes_palettes$Darjeeling2[2],
    "Hispanics" = wesanderson::wes_palettes$Darjeeling2[3],
    "Black" = wesanderson::wes_palettes$Darjeeling2[4]
  )



microbiome_genus_filt <- function(gut_microbiome_table, gut_microbiome_tax, prevalence) {
  # 计算原始表中的列数
  sample_num <- length(colnames(gut_microbiome_table))
  
  # 将tax数据的前6列绑定到microbiome表上
  gut_microbiome_table <- cbind(gut_microbiome_table, gut_microbiome_tax[, 1:6])
  
  # 使用aggregate函数按Genus对数据进行聚合求和
  gut_microbiome_table <- aggregate(gut_microbiome_table[, 1:sample_num], 
                                    by = list(gut_microbiome_table$Genus), 
                                    sum)
  
  # 设置行名为Genus名
  row.names(gut_microbiome_table) <- gut_microbiome_table$Group.1
  
  # 移除由aggregate自动创建的'Group.1'列
  gut_microbiome_table <- gut_microbiome_table[, -1]
  
  # 转置微生物组表，以便进行OTU存在性分析
  gut_microbiome_table <- t(gut_microbiome_table)
  
  # 计算每个OTU的存在性
  otu_presence <- colSums(gut_microbiome_table > 0)
  
  # 计算存在性阈值
  threshold <- prevalence * sample_num
  
  # 基于阈值过滤OTU
  gut_microbiome_table <- gut_microbiome_table[, otu_presence >= threshold]
  
  # 将处理后的表转换为数据框
  gut_microbiome_table <- data.frame(gut_microbiome_table)
  
  # 返回处理后的表
  return(gut_microbiome_table)
}


phylum_name =
  c(
    "Actinobacteria",
    "Bacteroidetes",
    "Cyanobacteria/Chloroplast"  ,
    "Firmicutes",
    "Lentisphaerae",
    "Proteobacteria",
    "Synergistetes",
    "Verrucomicrobia",
    "Campilobacterota",
    "Candidatus_Saccharibacteria",
    "Fusobacteria",
    "Plantae",
    "Tenericutes",
    "Spirochaetes"
  )

phylum_color =
  ggsci::pal_simpsons()(n = length(phylum_name))

site_color <-
  c(gut = "#edd064" , oral = "#a1d5b9" , skin = "#f2ccac", nasal = "#a17db4")


# 可视化分组比较结果
plot_group_comparison_results <- function(results) {
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(reshape2)
  
  summary_df <- results$summary
  group1_name <- results$group1_name
  group2_name <- results$group2_name
  
  # 1. R²值比较散点图
  p1 <- ggplot(summary_df, aes(x = group1_r2, y = group2_r2)) +
    geom_point(aes(color = p_diff_adjusted < 0.05), alpha = 0.7, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"),
                       labels = c("Not significant", "Significant (p < 0.05)")) +
    theme_bw() +
    labs(title = paste("R² Comparison:", group1_name, "vs", group2_name),
         x = paste("R² in", group1_name, "group"),
         y = paste("R² in", group2_name, "group"),
         color = "Group Difference") +
    coord_equal()
  
  # 2. R²差异的分布图
  p2 <- ggplot(summary_df, aes(x = r2_difference)) +
    geom_histogram(bins = 30, alpha = 0.7, fill = "steelblue", color = 'black') +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    theme_bw() +
    labs(title = "Distribution of R² Differences",
         x = paste("R² Difference (", group1_name, "-", group2_name, ")"),
         y = "Count")
  
  # 3. 交互特征数量比较
  interaction_data <- summary_df %>%
    select(metabolite, group1_n_interaction, group2_n_interaction) %>%
    melt(id.vars = "metabolite", 
         variable.name = "group", 
         value.name = "n_interactions") %>%
    mutate(group = ifelse(group == "group1_n_interaction", group1_name, group2_name))
  
  interaction_data$group <- factor(interaction_data$group, levels = c("IS", "IR"))
  
  p3 <- ggplot(interaction_data, aes(x = group, y = n_interactions, fill = group)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    scale_fill_manual(values=c("IR"="#E69F00", "IS"="#0072B2")) +
    theme_bw() +
    labs(title = "Number of Selected Interaction Features",
         x = "Group",
         y = "Number of Interaction Features",
         fill = "Group")+stat_compare_means(label = "p.format")
  
  # 4. 效应大小vs显著性
  p4 <- ggplot(summary_df, aes(x = effect_size, y = -log10(p_diff_adjusted))) +
    geom_point(aes(fill = abs(effect_size) > 0.5 & p_diff_adjusted < 0.05), 
               alpha = 0.7, size = 4, shape = 21) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray") +
    scale_fill_manual(values = c("FALSE" = "gray", "TRUE" = "red"),
                      labels = c("Not significant", "Significant & Large effect")) +
    theme_bw() +
    labs(title = "Effect Size vs Significance",
         x = "Effect Size (Cohen's d)",
         y = "-log10(Adjusted p-value)",
         fill = "Classification")
  
  # 5. 火山图样式的R²差异图
  p5 <- ggplot(summary_df, aes(x = r2_difference, y = -log10(p_diff_adjusted))) +
    geom_point(aes(fill = p_diff_adjusted < 0.05 & abs(r2_difference) > 0.1), 
               alpha = 0.7, size = 4, shape = 21) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "gray") +
    scale_fill_manual(values = c("FALSE" = "gray", "TRUE" = "red"),
                      labels = c("Not significant", "Significant difference")) +
    theme_bw() +
    labs(title = "Volcano Plot: R² Differences",
         x = paste("R² Difference (", group1_name, "-", group2_name, ")"),
         y = "-log10(Adjusted p-value)",
         color = "Significance")
  
  # 6. 相关性热图：交互特征数量与R²表现
  correlation_data <- summary_df %>%
    select(group1_r2, group2_r2, group1_n_interaction, group2_n_interaction, 
           r2_difference, interaction_difference) %>%
    cor(use = "complete.obs")
  
  correlation_melted <- melt(correlation_data)
  
  p6 <- ggplot(correlation_melted, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    geom_text(aes(label = round(value, 2)), color = "white", size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limit = c(-1, 1)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Correlation Matrix",
         x = "", y = "", fill = "Correlation")
  
  # 组合所有图表
  combined_plots <- (p1 + p2) / (p3 + p4) / (p5 + p6) +
    plot_layout(heights = c(1, 1, 1))
  
  return(list(
    r2_comparison = p1,
    r2_difference_dist = p2,
    interaction_features = p3,
    effect_size_significance = p4,
    volcano_plot = p5,
    correlation_heatmap = p6,
    combined = combined_plots
  ))
}


plot_metabolites_mirror_style <- function(significant_metabolites, top_n = 30) {
  library(ggplot2)
  library(dplyr)
  
  # 准备数据
  if(nrow(significant_metabolites) > top_n) {
    plot_data <- significant_metabolites %>%
      arrange(desc(abs(r2_difference))) %>%
      head(top_n)
  } else {
    plot_data <- significant_metabolites %>%
      arrange(desc(abs(r2_difference)))
  }
  
  # 按原图风格排序：先正值（降序），后负值（从小到大，即从最负开始）
  positive_data <- plot_data[plot_data$r2_difference > 0, ] %>%
    arrange(desc(r2_difference))
  negative_data <- plot_data[plot_data$r2_difference < 0, ] %>%
    arrange(desc(r2_difference))  # 负值从最小（最负）到最大（接近0）
  
  # 重新组合数据
  plot_data_ordered <- rbind(positive_data, negative_data)
  
  # 创建x轴位置
  plot_data_ordered$x_pos <- 1:nrow(plot_data_ordered)
  
  # 为双向显示创建数据
  plot_data_ordered$y_upper <- ifelse(plot_data_ordered$r2_difference > 0, 
                                      plot_data_ordered$r2_difference, 0)
  plot_data_ordered$y_lower <- ifelse(plot_data_ordered$r2_difference < 0, 
                                      -plot_data_ordered$r2_difference, 0)  # 取绝对值显示在下方
  
  # 找到最大值用于设置y轴
  max_val <- max(abs(plot_data_ordered$r2_difference))
  
  plot_data_ordered <-
    plot_data_ordered %>% 
    dplyr::arrange(r2_difference)
  
  p <- ggplot(plot_data_ordered, aes(x = x_pos)) +
    # 上方条形（IR优势，绿色）
    geom_col(aes(y = y_upper), fill = "#E69F00", alpha = 0.9, width = 0.8) +
    # 下方条形（IS优势，红色，向下显示）
    geom_col(aes(y = -y_lower), fill = "#0072B2", alpha = 0.9, width = 0.8) +
    # 零线
    geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
    # 设置y轴范围
    scale_y_continuous(
      limits = c(-max_val * 1.1, max_val * 1.1),
      breaks = seq(-max_val, max_val, length.out = 7),
      labels = function(x) sprintf("%.1f", abs(x))  # 显示绝对值
    ) +
    # x轴设置
    scale_x_continuous(
      breaks = plot_data_ordered$x_pos,
      labels = plot_data_ordered$HMDB.Name,
      expand = c(0.01, 0.01)
    ) +
    # 主题
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "gray90", size = 0.3),
      axis.line.x = element_line(color = "black", size = 0.5),
      plot.margin = margin(t = 20, r = 10, b = 10, l = 10)
    ) +
    # 标签
    labs(
      title = "Metabolite R² Differences: IR vs IS Groups",
      x = "",
      y = "R² Difference",
      caption = ""
    ) +
    # 添加颜色图例标注
    annotate("text", x = length(plot_data_ordered$metabolite) * 0.05, 
             y = max_val * 1, label = "IR", 
             color = "#E69F00", size = 4, fontface = "bold") +
    annotate("text", x = length(plot_data_ordered$metabolite) * 0.05, 
             y = -max_val * 0.9, label = "IS", 
             color = "#0072B2", size = 4, fontface = "bold")
  
  return(p)
}
