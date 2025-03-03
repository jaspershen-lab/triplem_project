plot_microbe_interaction_basic <- function(gut_data, oral_data, metabolite_data, 
                                           gut_feature, oral_feature, metabolite,
                                           use_log = FALSE) {
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
  
  # 计算相关性
  gut_cor <- cor(gut_values, metabolite_values, method = "spearman")
  oral_cor <- cor(oral_values, metabolite_values, method = "spearman")
  
  # 创建散点图
  p <- ggplot(plot_data, aes(x = gut, y = oral, color = metabolite)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(metabolite_values)) +
    labs(
      title = paste("Interaction between", gut_feature, "and", oral_feature),
      subtitle = paste0("Impact on metabolite: ", metabolite, 
                        "\nGut correlation: ", round(gut_cor, 2),
                        ", Oral correlation: ", round(oral_cor, 2)),
      x = paste0(gut_feature, if(use_log) " (log10)" else ""),
      y = paste0(oral_feature, if(use_log) " (log10)" else ""),
      color = "Metabolite\nLevel"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      legend.position = "right"
    )
  
  return(p)
}

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
    )+facet_wrap(~oral_quantile)+xlim(c(0,2.5))
  
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

# 4. 3D散点图 - 使用plotly创建交互式3D图
plot_microbe_interaction_3d <- function(gut_data, oral_data, metabolite_data,
                                        gut_feature, oral_feature, metabolite,
                                        use_log = FALSE) {
  # 需要安装并加载plotly
  if(!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is needed for this function. Please install it.")
  }
  require(plotly)
  
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
    sample = common_samples,
    gut = gut_values,
    oral = oral_values,
    metabolite = metabolite_values
  )
  
  # 创建3D散点图
  p <- plot_ly(plot_data, x = ~gut, y = ~oral, z = ~metabolite, 
               type = "scatter3d", mode = "markers",
               marker = list(size = 5, color = ~metabolite, colorscale = "Viridis",
                             colorbar = list(title = "Metabolite Level")),
               text = ~sample,
               hoverinfo = "text") %>%
    layout(
      title = paste("3D Interaction Plot of", gut_feature, "and", oral_feature, "on", metabolite),
      scene = list(
        xaxis = list(title = paste0(gut_feature, if(use_log) " (log10)" else "")),
        yaxis = list(title = paste0(oral_feature, if(use_log) " (log10)" else "")),
        zaxis = list(title = metabolite)
      )
    )
  
  return(p)
}

# 5. 创建交互响应面图
plot_microbe_interaction_surface <- function(gut_data, oral_data, metabolite_data,
                                             gut_feature, oral_feature, metabolite,
                                             use_log = FALSE) {
  # 需要安装并加载plotly
  if(!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is needed for this function. Please install it.")
  }
  require(plotly)
  
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
  
  # 拟合一个平滑模型来创建响应面
  # 使用广义加性模型或局部多项式回归
  if(requireNamespace("mgcv", quietly = TRUE)) {
    require(mgcv)
    model <- gam(metabolite ~ s(gut, oral), data = plot_data)
    
    # 创建预测网格
    gut_seq <- seq(min(plot_data$gut), max(plot_data$gut), length.out = 50)
    oral_seq <- seq(min(plot_data$oral), max(plot_data$oral), length.out = 50)
    
    grid_data <- expand.grid(gut = gut_seq, oral = oral_seq)
    grid_data$predicted <- predict(model, newdata = grid_data)
    
    # 重塑为矩阵形式
    z_matrix <- matrix(grid_data$predicted, nrow = length(gut_seq), ncol = length(oral_seq))
    
    # 创建响应面
    p <- plot_ly() %>%
      add_trace(
        x = gut_seq,
        y = oral_seq,
        z = z_matrix,
        type = "surface",
        colorscale = "Viridis",
        contours = list(
          z = list(
            show = TRUE,
            usecolormap = TRUE,
            highlightcolor = "#ff0000",
            project = list(z = TRUE)
          )
        )
      ) %>%
      add_trace(
        data = plot_data,
        x = ~gut, y = ~oral, z = ~metabolite,
        type = "scatter3d",
        mode = "markers",
        marker = list(size = 3, color = "red", symbol = 104),
        showlegend = FALSE
      ) %>%
      layout(
        title = paste("Response Surface for Interaction of", gut_feature, "and", oral_feature),
        scene = list(
          xaxis = list(title = paste0(gut_feature, if(use_log) " (log10)" else "")),
          yaxis = list(title = paste0(oral_feature, if(use_log) " (log10)" else "")),
          zaxis = list(title = metabolite)
        )
      )
  } else {
    # 如果没有mgcv包，则使用loess
    # 创建一个网格
    gut_seq <- seq(min(plot_data$gut), max(plot_data$gut), length.out = 50)
    oral_seq <- seq(min(plot_data$oral), max(plot_data$oral), length.out = 50)
    
    grid_data <- expand.grid(gut = gut_seq, oral = oral_seq)
    
    # 拟合loess模型
    loess_model <- loess(metabolite ~ gut * oral, data = plot_data, span = 0.5)
    grid_data$predicted <- predict(loess_model, newdata = grid_data)
    
    # 重塑为矩阵形式
    z_matrix <- matrix(grid_data$predicted, nrow = length(gut_seq), ncol = length(oral_seq))
    
    # 创建响应面
    p <- plot_ly() %>%
      add_trace(
        x = gut_seq,
        y = oral_seq,
        z = z_matrix,
        type = "surface",
        colorscale = "Viridis",
        contours = list(
          z = list(
            show = TRUE,
            usecolormap = TRUE,
            highlightcolor = "#ff0000",
            project = list(z = TRUE)
          )
        )
      ) %>%
      add_trace(
        data = plot_data,
        x = ~gut, y = ~oral, z = ~metabolite,
        type = "scatter3d",
        mode = "markers",
        marker = list(size = 3, color = "red", symbol = 104),
        showlegend = FALSE
      ) %>%
      layout(
        title = paste("Response Surface for Interaction of", gut_feature, "and", oral_feature),
        scene = list(
          xaxis = list(title = paste0(gut_feature, if(use_log) " (log10)" else "")),
          yaxis = list(title = paste0(oral_feature, if(use_log) " (log10)" else "")),
          zaxis = list(title = metabolite)
        )
      )
  }
  
  return(p)
}

# 6. 创建差异效应散点图
plot_microbe_differential_effect <- function(gut_data, oral_data, metabolite_data,
                                             gut_feature, oral_feature, metabolite,
                                             n_quantiles = 3, use_log = FALSE) {
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
    sample = common_samples,
    gut = gut_values,
    oral = oral_values,
    metabolite = metabolite_values
  )
  
  # 根据肠道菌群分位数分组
  gut_breaks <- quantile(gut_values, probs = seq(0, 1, length.out = n_quantiles + 1))
  plot_data$gut_quantile <- cut(gut_values, 
                                breaks = gut_breaks, 
                                labels = paste0("G", 1:n_quantiles),
                                include.lowest = TRUE)
  
  # 根据口腔菌群分位数分组
  oral_breaks <- quantile(oral_values, probs = seq(0, 1, length.out = n_quantiles + 1))
  plot_data$oral_quantile <- cut(oral_values, 
                                 breaks = oral_breaks, 
                                 labels = paste0("O", 1:n_quantiles),
                                 include.lowest = TRUE)
  
  # 为肠道和口腔菌群的每个分位数组合计算平均代谢物水平
  mean_data <- aggregate(metabolite ~ gut_quantile + oral_quantile, 
                         data = plot_data, 
                         FUN = mean)
  
  # 计算最高和最低组合的差异
  spread_data <- matrix(mean_data$metabolite, 
                        nrow = n_quantiles, 
                        ncol = n_quantiles,
                        dimnames = list(levels(plot_data$gut_quantile), 
                                        levels(plot_data$oral_quantile)))
  
  # 计算效应差异
  effect_range <- max(mean_data$metabolite) - min(mean_data$metabolite)
  
  # 计算肠道菌群在不同口腔菌群分位数下的效应差异
  gut_effects <- lapply(levels(plot_data$oral_quantile), function(oq) {
    subset_data <- mean_data[mean_data$oral_quantile == oq, ]
    max_effect <- max(subset_data$metabolite) - min(subset_data$metabolite)
    data.frame(
      oral_quantile = oq,
      gut_effect = max_effect,
      gut_effect_percent = max_effect / effect_range * 100
    )
  })
  gut_effects <- do.call(rbind, gut_effects)
  
  # 计算口腔菌群在不同肠道菌群分位数下的效应差异
  oral_effects <- lapply(levels(plot_data$gut_quantile), function(gq) {
    subset_data <- mean_data[mean_data$gut_quantile == gq, ]
    max_effect <- max(subset_data$metabolite) - min(subset_data$metabolite)
    data.frame(
      gut_quantile = gq,
      oral_effect = max_effect,
      oral_effect_percent = max_effect / effect_range * 100
    )
  })
  oral_effects <- do.call(rbind, oral_effects)
  
  # 创建差异效应散点图
  p1 <- ggplot(gut_effects, aes(x = oral_quantile, y = gut_effect_percent)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", gut_effect_percent)), 
              vjust = -0.5, size = 3) +
    labs(
      title = paste("Effect of", gut_feature, "across", oral_feature, "quantiles"),
      subtitle = paste0("Impact on metabolite: ", metabolite),
      x = paste0(oral_feature, " Quantile"),
      y = "Gut Effect (% of Total Range)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
  
  p2 <- ggplot(oral_effects, aes(x = gut_quantile, y = oral_effect_percent)) +
    geom_bar(stat = "identity", fill = "coral", alpha = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", oral_effect_percent)), 
              vjust = -0.5, size = 3) +
    labs(
      title = paste("Effect of", oral_feature, "across", gut_feature, "quantiles"),
      subtitle = paste0("Impact on metabolite: ", metabolite),
      x = paste0(gut_feature, " Quantile"),
      y = "Oral Effect (% of Total Range)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
  
  # 使用patchwork组合图
  if(requireNamespace("patchwork", quietly = TRUE)) {
    require(patchwork)
    p <- p1 / p2
  } else {
    # 如果没有patchwork，返回列表
    p <- list(gut_effect = p1, oral_effect = p2)
  }
  
  return(p)
}

 7. 使用示例
 plot_examples <- function() {
   gut_feature <- "Bacteroides"
   oral_feature <- "Streptococcus"
   metabolite <- "butyrate"
   
  # 基本散点图
  p1 <- plot_microbe_interaction_basic(gut_data, oral_data, metabolite_data,
                                       gut_feature, oral_feature, metabolite)
   
 # 分位数散点图
   p2 <- plot_microbe_interaction_quantile(gut_data, oral_data, metabolite_data,
                                          gut_feature, oral_feature, metabolite)
   
   # 交互热图
   p3 <- plot_microbe_interaction_heatmap(gut_data, oral_data, metabolite_data,
                                         gut_feature, oral_feature, metabolite)
  
   # 差异效应图
   p4 <- plot_microbe_differential_effect(gut_data, oral_data, metabolite_data,
                                         gut_feature, oral_feature, metabolite)
   
   # 返回图列表
   return(list(basic = p1, quantile = p2, heatmap = p3, differential = p4))
 }