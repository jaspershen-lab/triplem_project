# 读取interaction的结果


gut_oral_interaction<-load("combined_results_with_interactions")



# 2. 分位数分割图 - 根据一个菌群的分位数来展示另一个菌群与代谢物的关系
plot_microbe_interaction_quantile <- function(gut_data, oral_data, metabolite_data,
                                              gut_feature, oral_feature, metabolite,
                                              n_quantiles = 3, use_log = FALSE) {
  # 获取共同样本
  gut_samples <- colnames(gut_data)
  oral_samples <- colnames(oral_data)
  meta_samples <- colnames(metabolite_data)
  common_samples <- Reduce(intersect, list(gut_samples, oral_samples, meta_samples))
  
  # 提取数据
  gut_values <- as.numeric(gut_data[gut_feature, common_samples])
  oral_values <- as.numeric(oral_data[oral_feature, common_samples])
  metabolite_values <- as.numeric(metabolite_data[metabolite, common_samples])
  
  # 可选对数转换
  if(use_log) {
    # 处理零值
    gut_values[gut_values == 0] <- min(gut_values[gut_values > 0]) / 2
    oral_values[oral_values == 0] <- min(oral_values[oral_values > 0]) / 2
    
    gut_values <- log10(gut_values)
    oral_values <- log10(oral_values)
  }
  
  # 创建数据框
  plot_data <- data.frame(
    gut = gut_values,
    oral = oral_values,
    metabolite = metabolite_values
  )
  
  # 根据肠道菌群分位数分组
  quantile_breaks <- quantile(oral_values, probs = seq(0, 1, length.out = n_quantiles + 1))
  plot_data$oral_quantile <- cut(oral_values, 
                                breaks = quantile_breaks, 
                                labels = paste0("Q", 1:n_quantiles),
                                include.lowest = TRUE)
  
  # 对于每个分位数计算口腔菌群与代谢物的相关性
  quantile_cors <- lapply(levels(plot_data$oral_quantile), function(q) {
    subset_data <- plot_data[plot_data$oral_quantile == q, ]
    cor_value <- cor(subset_data$gut, subset_data$metabolite, method = "spearman")
    data.frame(quantile = q, correlation = cor_value)
  })
  quantile_cors <- do.call(rbind, quantile_cors)
  
  # 创建散点图
  p <- ggplot(plot_data, aes(x = gut, y = metabolite, color = oral_quantile)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, aes(group = oral_quantile)) +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = paste("Interaction effect of", gut_feature, "and", oral_feature),
      subtitle = paste0("On metabolite: ", metabolite,
                        "\nCorrelation by quantile: ", 
                        paste(quantile_cors$quantile, ":", round(quantile_cors$correlation, 2), collapse = ", ")),
      x = paste0(oral_feature, if(use_log) " (log10)" else ""),
      y = paste0(metabolite),
      color = paste0(gut_feature, "\nQuantile")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      legend.position = "right"
    )+facet_wrap(~oral_quantile)
  
  return(p)
}

# 3. 交互热图 - 将两个菌群按分位数分割，展示每个组合对应的代谢物水平
plot_microbe_interaction_heatmap <- function(gut_data, oral_data, metabolite_data,
                                             gut_feature, oral_feature, metabolite,
                                             n_quantiles = 4, use_log = FALSE) {
  # 获取共同样本
  gut_samples <- colnames(gut_data)
  oral_samples <- colnames(oral_data)
  meta_samples <- colnames(metabolite_data)
  common_samples <- Reduce(intersect, list(gut_samples, oral_samples, meta_samples))
  
  # 提取数据
  gut_values <- gut_data[gut_feature, common_samples]
  oral_values <- oral_data[oral_feature, common_samples]
  metabolite_values <- metabolite_data[metabolite, common_samples]
  
  # 可选对数转换
  if(use_log) {
    # 处理零值
    gut_values[gut_values == 0] <- min(gut_values[gut_values > 0]) / 2
    oral_values[oral_values == 0] <- min(oral_values[oral_values > 0]) / 2
    
    gut_values <- log10(gut_values)
    oral_values <- log10(oral_values)
  }
  
  # 创建数据框
  plot_data <- data.frame(
    gut = gut_values,
    oral = oral_values,
    metabolite = metabolite_values
  )
  
  # 根据分位数分组
  gut_breaks <- quantile(gut_values, probs = seq(0, 1, length.out = n_quantiles + 1))
  oral_breaks <- quantile(oral_values, probs = seq(0, 1, length.out = n_quantiles + 1))
  
  plot_data$gut_quantile <- cut(gut_values, 
                                breaks = gut_breaks, 
                                labels = paste0("G", 1:n_quantiles),
                                include.lowest = TRUE)
  
  plot_data$oral_quantile <- cut(oral_values, 
                                 breaks = oral_breaks, 
                                 labels = paste0("O", 1:n_quantiles),
                                 include.lowest = TRUE)
  
  # 计算每个组合的平均代谢物水平
  heatmap_data <- aggregate(metabolite ~ gut_quantile + oral_quantile, 
                            data = plot_data, 
                            FUN = mean)
  
  # 计算组合样本数量
  count_data <- as.data.frame(table(plot_data$gut_quantile, plot_data$oral_quantile))
  names(count_data) <- c("gut_quantile", "oral_quantile", "count")
  
  # 合并数据
  heatmap_data <- merge(heatmap_data, count_data, by = c("gut_quantile", "oral_quantile"))
  
  # 创建热图
  p <- ggplot(heatmap_data, aes(x = oral_quantile, y = gut_quantile, fill = metabolite)) +
    geom_tile(color = "white") +
    geom_text(aes(label = paste0(round(metabolite, 2), "\n(n=", count, ")")), 
              size = 3, color = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = median(heatmap_data$metabolite)) +
    labs(
      title = paste("Interaction effect between", gut_feature, "and", oral_feature),
      subtitle = paste0("Average metabolite level (", metabolite, ") by quantile groups"),
      x = paste0(oral_feature, " Quantiles"),
      y = paste0(gut_feature, " Quantiles"),
      fill = "Mean\nMetabolite\nLevel"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "right"
    )
  
  return(p)
}



