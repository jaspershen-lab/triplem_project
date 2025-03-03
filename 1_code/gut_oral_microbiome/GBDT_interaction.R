# 新增函数：创建肠道菌群和口腔菌群的交互特征
create_microbiome_interactions <- function(gut_data, oral_data, 
                                           method = "multiplication",
                                           feature_selection = TRUE,
                                           max_features = 50,
                                           correlation_threshold = 0.2) {
  # 获取样本ID并确保数据格式一致
  gut_samples <- colnames(gut_data)
  oral_samples <- colnames(oral_data)
  
  # 确保两个数据集有相同的样本
  common_samples <- intersect(gut_samples, oral_samples)
  message(sprintf("Common samples for interaction analysis: %d", length(common_samples)))
  
  # 提取共同样本的数据
  gut_matched <- gut_data[, common_samples]
  oral_matched <- oral_data[, common_samples]
  
  # 转置数据，使行为样本，列为特征
  gut_df <- as.data.frame(t(gut_matched))
  oral_df <- as.data.frame(t(oral_matched))
  
  # 添加前缀
  colnames(gut_df) <- paste0("gut_", colnames(gut_df))
  colnames(oral_df) <- paste0("oral_", colnames(oral_df))
  
  # 如果需要特征选择，则选择最重要的特征
  if (feature_selection) {
    # 对于肠道特征，可以基于方差或平均丰度选择
    gut_variances <- apply(gut_df, 2, var)
    gut_means <- apply(gut_df, 2, mean)
    gut_importance <- gut_variances * gut_means  # 结合方差和平均丰度
    top_gut_indices <- order(gut_importance, decreasing = TRUE)[1:min(max_features, ncol(gut_df))]
    selected_gut_features <- colnames(gut_df)[top_gut_indices]
    
    # 对于口腔特征，同样基于方差或平均丰度选择
    oral_variances <- apply(oral_df, 2, var)
    oral_means <- apply(oral_df, 2, mean)
    oral_importance <- oral_variances * oral_means
    top_oral_indices <- order(oral_importance, decreasing = TRUE)[1:min(max_features, ncol(oral_df))]
    selected_oral_features <- colnames(oral_df)[top_oral_indices]
    
    message(sprintf("Selected %d gut features and %d oral features for interaction analysis", 
                    length(selected_gut_features), length(selected_oral_features)))
    
    # 筛选数据框
    gut_df_selected <- gut_df[, selected_gut_features, drop = FALSE]
    oral_df_selected <- oral_df[, selected_oral_features, drop = FALSE]
  } else {
    gut_df_selected <- gut_df
    oral_df_selected <- oral_df
    selected_gut_features <- colnames(gut_df)
    selected_oral_features <- colnames(oral_df)
  }
  
  # 创建交互特征
  message("Creating interaction features...")
  interaction_list <- list()
  interaction_names <- character()
  
  # 计算特征之间的相关性，只选择相关的特征对创建交互
  if (correlation_threshold > 0) {
    # 合并数据计算相关性
    combined_df <- cbind(gut_df_selected, oral_df_selected)
    correlation_matrix <- cor(combined_df, method = "spearman")
    
    # 只保留肠道-口腔特征对之间的相关性
    gut_indices <- which(grepl("^gut_", colnames(combined_df)))
    oral_indices <- which(grepl("^oral_", colnames(combined_df)))
    
    # 创建一个计数器来跟踪添加的交互特征
    interaction_count <- 0
    
    # 遍历所有肠道和口腔特征对
    for (i in gut_indices) {
      gut_feature <- colnames(combined_df)[i]
      for (j in oral_indices) {
        oral_feature <- colnames(combined_df)[j]
        # 检查相关性是否超过阈值
        if (abs(correlation_matrix[i, j]) >= correlation_threshold) {
          # 根据指定方法创建交互项
          if (method == "multiplication") {
            interaction_feature <- gut_df_selected[, gut_feature] * oral_df_selected[, oral_feature]
          } else if (method == "ratio") {
            # 避免除以零
            denominator <- oral_df_selected[, oral_feature]
            denominator[denominator == 0] <- min(denominator[denominator > 0]) / 10
            interaction_feature <- gut_df_selected[, gut_feature] / denominator
          } else if (method == "difference") {
            interaction_feature <- gut_df_selected[, gut_feature] - oral_df_selected[, oral_feature]
          } else {
            stop("Unsupported interaction method. Use 'multiplication', 'ratio', or 'difference'.")
          }
          
          # 将交互特征添加到列表
          interaction_name <- paste0("int_", gsub("gut_", "", gut_feature), "_", 
                                     gsub("oral_", "", oral_feature))
          interaction_list[[interaction_name]] <- interaction_feature
          interaction_names <- c(interaction_names, interaction_name)
          interaction_count <- interaction_count + 1
        }
      }
    }
    
    message(sprintf("Created %d interaction features based on correlation threshold %.2f", 
                    interaction_count, correlation_threshold))
  } else {
    # 不使用相关性筛选，创建所有可能的交互项
    for (gut_feature in selected_gut_features) {
      for (oral_feature in selected_oral_features) {
        # 根据指定方法创建交互项
        if (method == "multiplication") {
          interaction_feature <- gut_df_selected[, gut_feature] * oral_df_selected[, oral_feature]
        } else if (method == "ratio") {
          # 避免除以零
          denominator <- oral_df_selected[, oral_feature]
          denominator[denominator == 0] <- min(denominator[denominator > 0]) / 10
          interaction_feature <- gut_df_selected[, gut_feature] / denominator
        } else if (method == "difference") {
          interaction_feature <- gut_df_selected[, gut_feature] - oral_df_selected[, oral_feature]
        } else {
          stop("Unsupported interaction method. Use 'multiplication', 'ratio', or 'difference'.")
        }
        
        # 将交互特征添加到列表
        interaction_name <- paste0("int_", gsub("gut_", "", gut_feature), "_", 
                                   gsub("oral_", "", oral_feature))
        interaction_list[[interaction_name]] <- interaction_feature
        interaction_names <- c(interaction_names, interaction_name)
      }
    }
    
    message(sprintf("Created %d interaction features between %d gut features and %d oral features", 
                    length(interaction_list), length(selected_gut_features), 
                    length(selected_oral_features)))
  }
  
  # 将交互特征列表转换为数据框
  interaction_df <- as.data.frame(interaction_list)
  
  # 确保行名一致
  rownames(interaction_df) <- rownames(gut_df)
  
  # 返回交互特征数据框以及原始筛选后的数据
  return(list(
    interaction_features = interaction_df,
    gut_features = gut_df_selected,
    oral_features = oral_df_selected,
    interaction_method = method,
    feature_names = list(
      gut = selected_gut_features,
      oral = selected_oral_features,
      interaction = interaction_names
    )
  ))
}

# 修改现有的preprocess_combined_data函数，添加交互特征支持
preprocess_combined_data_with_interactions <- function(gut_data, oral_data, metabolite_data,
                                                       include_interactions = TRUE,
                                                       interaction_method = "multiplication",
                                                       max_features = 30,
                                                       correlation_threshold = 0.3) {
  # 获取样本ID
  gut_samples <- colnames(gut_data)
  oral_samples <- colnames(oral_data)
  meta_samples <- colnames(metabolite_data)
  
  # 检查样本ID
  message("Initial sample counts:")
  message(sprintf("Gut microbiome samples: %d", length(gut_samples)))
  message(sprintf("Oral microbiome samples: %d", length(oral_samples)))
  message(sprintf("Metabolite samples: %d", length(meta_samples)))
  
  # 找出共同样本
  common_samples <- Reduce(intersect, list(gut_samples, oral_samples, meta_samples))
  message(sprintf("Common samples across all datasets: %d", length(common_samples)))
  
  # 提取共同样本的数据
  gut_matched <- gut_data[, common_samples]
  oral_matched <- oral_data[, common_samples]
  metabolite_matched <- metabolite_data[, common_samples]
  
  # 转换和检查数据维度
  gut_df <- as.data.frame(t(gut_matched))
  oral_df <- as.data.frame(t(oral_matched))
  
  # 添加来源前缀
  colnames(gut_df) <- paste0("gut_", colnames(gut_df))
  colnames(oral_df) <- paste0("oral_", colnames(oral_df))
  
  # 检查维度
  message(sprintf("Gut features: %d", ncol(gut_df)))
  message(sprintf("Oral features: %d", ncol(oral_df)))
  
  # 基本的微生物组合并
  combined_microbiome <- cbind(gut_df, oral_df)
  
  # 添加交互特征（如果需要）
  if (include_interactions) {
    interactions_result <- create_microbiome_interactions(
      gut_data = gut_matched,
      oral_data = oral_matched,
      method = interaction_method,
      feature_selection = TRUE,
      max_features = max_features,
      correlation_threshold = correlation_threshold
    )
    
    # 合并交互特征
    interaction_df <- interactions_result$interaction_features
    message(sprintf("Adding %d interaction features", ncol(interaction_df)))
    
    # 确保行名匹配
    if (!identical(rownames(combined_microbiome), rownames(interaction_df))) {
      warning("Row names (samples) do not match between microbiome data and interaction features!")
    }
    
    # 合并所有特征
    combined_microbiome <- cbind(combined_microbiome, interaction_df)
  }
  
  # 确保样本顺序一致
  rownames(combined_microbiome) <- common_samples
  
  return(list(
    combined_microbiome = combined_microbiome,
    metabolite = t(metabolite_matched),
    common_samples = common_samples,
    interactions_info = if(include_interactions) interactions_result else NULL
  ))
}

# 扩展analyze_combined_metabolite_ev函数，支持交互特征
analyze_combined_metabolite_ev_with_interactions <- function(gut_data, oral_data, metabolite_data, 
                                                             n_cores = NULL, seed = 42,
                                                             do_feature_selection = TRUE,
                                                             correlation_method = "spearman",
                                                             p_threshold = 0.05,
                                                             p_adjust_method = "BH",
                                                             rho_threshold = 0.1,
                                                             include_interactions = TRUE,
                                                             interaction_method = "multiplication",
                                                             max_interaction_features = 30,
                                                             interaction_correlation_threshold = 0.3) {
  # 设置全局种子
  set.seed(seed)
  
  # 设置参数
  n_boots <- 100
  n_folds <- 10
  
  # 设置并行
  if(is.null(n_cores)) {
    n_cores <- detectCores() - 1
  }
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # 数据预处理（包括交互特征）
  processed_data <- preprocess_combined_data_with_interactions(
    gut_data, oral_data, metabolite_data,
    include_interactions = include_interactions,
    interaction_method = interaction_method,
    max_features = max_interaction_features,
    correlation_threshold = interaction_correlation_threshold
  )
  
  # 准备数据
  X <- as.matrix(processed_data$combined_microbiome)
  Y <- processed_data$metabolite
  
  # 打印信息
  message("\nAnalysis settings:")
  message(sprintf("Feature selection: %s", if(do_feature_selection) "Yes" else "No"))
  if(do_feature_selection) {
    message(sprintf("Correlation method: %s", correlation_method))
    message(sprintf("P-value threshold: %.3f", p_threshold))
    message(sprintf("P-value adjustment method: %s", p_adjust_method))
    message(sprintf("Correlation coefficient threshold: %.2f", rho_threshold))
  }
  message(sprintf("Total features (combined microbes): %d", ncol(X)))
  message(sprintf("Gut features: %d", sum(grepl("^gut_", colnames(X)))))
  message(sprintf("Oral features: %d", sum(grepl("^oral_", colnames(X)))))
  message(sprintf("Interaction features: %d", sum(grepl("^int_", colnames(X)))))
  message(sprintf("Samples: %d", nrow(X)))
  message(sprintf("Metabolites: %d", ncol(Y)))
  message(sprintf("Using %d cores", n_cores))
  
  # GBDT参数
  gbdt_params <- list(
    n.trees = 100,
    interaction.depth = 15,
    shrinkage = 0.01,
    n.minobsinnode = 8,
    bag.fraction = 0.8,
    train.fraction = 0.8
  )
  
  # 导出函数到并行环境
  clusterExport(cl, c("single_cv", "calculate_r2"), envir = environment())
  
  # 存储结果
  results <- list()
  
  # 设置进度条
  pb <- progress_bar$new(
    format = "[:bar] :percent | Metabolite :current/:total | Elapsed: :elapsed | ETA: :eta",
    total = ncol(Y)
  )
  
  # 分析每个代谢物
  for(i in 1:ncol(Y)) {
    start_time <- Sys.time()
    current_y <- Y[, i]
    
    # 特征选择
    if(do_feature_selection) {
      feature_selection <- select_relevant_features(
        X, current_y,
        correlation_method = correlation_method,
        p_threshold = p_threshold,
        p_adjust_method = p_adjust_method,
        rho_threshold = rho_threshold
      )
      selected_features <- feature_selection$selected_features
      if (length(selected_features) > 0) {
        X_selected <- X[, selected_features, drop = FALSE]
      } else {
        X_selected <- X
      }
      
      # 统计选中的特征中来自不同来源的比例
      n_gut <- sum(grepl("^gut_", selected_features))
      n_oral <- sum(grepl("^oral_", selected_features))
      n_interaction <- sum(grepl("^int_", selected_features))
    } else {
      X_selected <- X
      n_gut <- sum(grepl("^gut_", colnames(X)))
      n_oral <- sum(grepl("^oral_", colnames(X)))
      n_interaction <- sum(grepl("^int_", colnames(X)))
    }
    
    # 只有在有选定特征时才继续分析
    if(ncol(X_selected) > 0) {
      # 并行Bootstrap
      boot_results <- foreach(b = 1:n_boots,
                              .combine = 'c',
                              .packages = c("gbm", "caret")) %dopar% {
                                local_seed <- seed * 1000 + b + i * n_boots
                                set.seed(local_seed)
                                boot_idx <- sample(1:nrow(X_selected), nrow(X_selected), replace = TRUE)
                                boot_X <- X_selected[boot_idx, , drop = FALSE]
                                boot_y <- current_y[boot_idx]
                                
                                single_cv(boot_X, boot_y, n_folds, gbdt_params, seed = local_seed)
                              }
      
      # 计算统计量
      mean_r2 <- mean(boot_results)
      ci <- quantile(boot_results, probs = c(0.025, 0.975))
      t_stat <- mean_r2 / (sd(boot_results) / sqrt(n_boots))
      p_value <- 2 * pt(-abs(t_stat), df = n_boots - 1)
      
      # 使用所有选定特征训练一个最终的GBDT模型，以计算特征重要性
      final_data <- data.frame(current_y, X_selected)
      final_model <- gbm(
        current_y ~ ., data = final_data,
        distribution = "gaussian",
        n.trees = gbdt_params$n.trees,
        interaction.depth = gbdt_params$interaction.depth,
        shrinkage = gbdt_params$shrinkage,
        n.minobsinnode = gbdt_params$n.minobsinnode,
        bag.fraction = gbdt_params$bag.fraction,
        train.fraction = gbdt_params$train.fraction,
        verbose = FALSE
      )
      feature_importance <- summary(final_model, plot = FALSE)
    } else {
      boot_results <- rep(0, n_boots)
      mean_r2 <- 0
      ci <- c(0, 0)
      p_value <- 1
      feature_importance <- NULL
    }
    
    # 存储结果
    results[[i]] <- list(
      metabolite = colnames(Y)[i],
      mean_r2 = mean_r2,
      ci_lower = ci[1],
      ci_upper = ci[2],
      p_value = p_value,
      all_r2s = boot_results,
      feature_selection = if(do_feature_selection) feature_selection else NULL,
      n_selected_features = if(do_feature_selection) length(selected_features) else ncol(X),
      n_gut_features = n_gut,
      n_oral_features = n_oral,
      n_interaction_features = n_interaction,
      feature_importance = feature_importance
    )
    
    # 更新进度条和打印结果
    pb$tick()
    end_time <- Sys.time()
    time_taken <- difftime(end_time, start_time, units = "mins")
    message(sprintf("\nMetabolite %d/%d (%s)", i, ncol(Y), colnames(Y)[i]))
    message(sprintf("Selected features: %d (Gut: %d, Oral: %d, Interaction: %d)", 
                    if(do_feature_selection) length(selected_features) else ncol(X),
                    n_gut, n_oral, n_interaction))
    message(sprintf("R² = %.3f (95%% CI: %.3f-%.3f, p = %.3e)", 
                    mean_r2, ci[1], ci[2], p_value))
    message(sprintf("Time: %.2f mins", time_taken))
  }
  
  # 关闭并行集群
  stopCluster(cl)
  
  # 整理结果为数据框
  summary_df <- do.call(rbind, lapply(results, function(x) {
    data.frame(
      metabolite = x$metabolite,
      r2_mean = x$mean_r2,
      r2_ci_lower = x$ci_lower,
      r2_ci_upper = x$ci_upper,
      p_value = x$p_value,
      n_selected_features = x$n_selected_features,
      n_gut_features = x$n_gut_features,
      n_oral_features = x$n_oral_features,
      n_interaction_features = x$n_interaction_features
    )
  }))
  
  # 添加交互特征分析函数
  analyze_interaction_importance <- function(results) {
    # 提取所有代谢物的特征重要性
    all_importance <- do.call(rbind, lapply(1:length(results), function(i) {
      imp <- results[[i]]$feature_importance
      if(!is.null(imp)) {
        data.frame(
          metabolite = results[[i]]$metabolite,
          feature = imp$var,
          importance = imp$rel.inf
        )
      }
    }))
    
    # 添加特征类型
    all_importance$feature_type <- "other"
    all_importance$feature_type[grepl("^gut_", all_importance$feature)] <- "gut"
    all_importance$feature_type[grepl("^oral_", all_importance$feature)] <- "oral"
    all_importance$feature_type[grepl("^int_", all_importance$feature)] <- "interaction"
    
    # 按代谢物和特征类型汇总重要性
    importance_by_type <- aggregate(importance ~ metabolite + feature_type, 
                                    data = all_importance, sum)
    
    # 提取交互特征详情
    interaction_details <- all_importance[all_importance$feature_type == "interaction", ]
    
    # 提取交互特征的组成部分
    if(nrow(interaction_details) > 0) {
      # 从交互特征名称中提取肠道和口腔特征
      pattern <- "int_(.+)_(.+)"
      interaction_parts <- do.call(rbind, lapply(interaction_details$feature, function(feat) {
        if(grepl(pattern, feat)) {
          parts <- regmatches(feat, regexec(pattern, feat))[[1]]
          data.frame(
            feature = feat,
            gut_part = parts[2],
            oral_part = parts[3]
          )
        } else {
          data.frame(
            feature = feat,
            gut_part = NA,
            oral_part = NA
          )
        }
      }))
      
      # 合并交互部分和重要性
      interaction_analysis <- merge(interaction_details, interaction_parts, by = "feature")
    } else {
      interaction_analysis <- data.frame()
    }
    
    return(list(
      importance_by_type = importance_by_type,
      interaction_details = interaction_analysis
    ))
  }
  
  # 扩展可视化函数
  plot_combined_feature_selection_results_with_interactions <- function(results) {
    # 准备数据
    feature_summary <- do.call(rbind, lapply(results$detailed_results, function(x) {
      data.frame(
        metabolite = x$metabolite,
        n_features = x$n_selected_features,
        n_gut = x$n_gut_features,
        n_oral = x$n_oral_features,
        n_interaction = x$n_interaction_features,
        r2 = x$mean_r2
      )
    }))
    
    # 分析交互特征重要性
    interaction_importance <- analyze_interaction_importance(results$detailed_results)
    
    # 1. 特征数量与R²的关系图
    p1 <- ggplot(feature_summary, aes(x = n_features, y = r2)) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "loess", se = TRUE) +
      theme_minimal() +
      labs(title = "Number of Selected Features vs R²",
           x = "Number of Selected Features",
           y = "R² Score")
    
    # 2. 特征类型比例图
    feature_types <- reshape2::melt(feature_summary[, c("metabolite", "n_gut", "n_oral", "n_interaction", "r2")],
                                    id.vars = c("metabolite", "r2"),
                                    variable.name = "feature_type",
                                    value.name = "count")
    feature_types$feature_type <- gsub("n_", "", feature_types$feature_type)
    
    p2 <- ggplot(feature_types, aes(x = feature_type, y = count, fill = feature_type)) +
      geom_boxplot(alpha = 0.7) +
      theme_minimal() +
      scale_fill_brewer(palette = "Set1") +
      labs(title = "Distribution of Selected Features by Type",
           x = "Feature Type",
           y = "Number of Features",
           fill = "Feature Type")
    
    # 3. 交互特征对模型性能的影响
    p3 <- ggplot(feature_summary, aes(x = n_interaction / n_features, y = r2)) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "loess", se = TRUE) +
      theme_minimal() +
      labs(title = "Proportion of Interaction Features vs R²",
           x = "Proportion of Interaction Features",
           y = "R² Score")
    
    # 4. 特征重要性按类型
    if(nrow(interaction_importance$importance_by_type) > 0) {
      p4 <- ggplot(interaction_importance$importance_by_type, 
                   aes(x = feature_type, y = importance, fill = feature_type)) +
        geom_boxplot(alpha = 0.7) +
        theme_minimal() +
        scale_fill_brewer(palette = "Set1") +
        labs(title = "Feature Importance by Type",
             x = "Feature Type",
             y = "Relative Importance (%)",
             fill = "Feature Type")
    } else {
      p4 <- ggplot() + 
        geom_text(aes(x = 0, y = 0, label = "No interaction features selected")) +
        theme_void() +
        labs(title = "Feature Importance by Type")
    }
    
    # 组合图表
    combined_plots <- (p1 + p2) / (p3 + p4) +
      plot_layout(heights = c(1, 1))
    
    return(list(
      feature_vs_r2 = p1,
      feature_distribution = p2,
      interaction_proportion = p3,
      importance_by_type = p4,
      combined = combined_plots,
      interaction_analysis = interaction_importance
    ))
  }
  
  # 添加可视化分析
  viz_results <- plot_combined_feature_selection_results_with_interactions(list(
    detailed_results = results,
    summary = summary_df
  ))
  
  return(list(
    detailed_results = results,
    summary = summary_df,
    common_samples = processed_data$common_samples,
    visualization = viz_results,
    interaction_info = processed_data$interactions_info
  ))
}

#使用示例
 combined_results_with_interactions <- analyze_combined_metabolite_ev_with_interactions(
   gut_data = gut_temp_object@expression_data,
   oral_data = oral_temp_object@expression_data,
   metabolite_data = metabolomics_temp_object@expression_data,
   do_feature_selection = TRUE,
   correlation_method = "spearman",
   p_threshold = 0.05,
   p_adjust_method = "none",
   rho_threshold = 0.1,
   include_interactions = TRUE,
   interaction_method = "multiplication",
   max_interaction_features = 30,
   interaction_correlation_threshold = 0.2
)