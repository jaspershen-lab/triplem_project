# 交叉验证函数（适配小样本）
single_cv <- function(X, y, n_folds = 5, gbdt_params, seed = NULL) {
  # 设置随机种子
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  # 检查样本数是否足够进行交叉验证
  n_samples <- length(y)
  if(n_samples < n_folds) {
    # 如果样本数少于折数，使用留一法交叉验证
    n_folds <- n_samples
  }
  
  # 检查是否有足够的样本进行训练
  min_train_samples <- max(gbdt_params$n.minobsinnode + 2, 5)
  if(n_samples < min_train_samples) {
    warning(paste("Too few samples for reliable cross-validation:", n_samples))
    return(0)  # 返回0表示无法训练
  }
  
  # 创建fold索引
  fold_ids <- sample(rep(1:n_folds, length.out = length(y)))
  all_predictions <- numeric(length(y))
  
  for(fold in 1:n_folds) {
    # 划分训练集和测试集
    test_idx <- which(fold_ids == fold)
    train_idx <- which(fold_ids != fold)
    
    # 确保训练集有足够的样本
    if(length(train_idx) < gbdt_params$n.minobsinnode + 1) {
      warning(paste("Insufficient training samples in fold", fold))
      all_predictions[test_idx] <- mean(y[train_idx])  # 使用均值预测
      next
    }
    
    # 准备数据
    train_data <- data.frame(X[train_idx, , drop = FALSE])
    test_data <- data.frame(X[test_idx, , drop = FALSE])
    train_data$target <- y[train_idx]
    
    # 检查特征数量，如果特征过多相对于样本数，进行降维
    if(ncol(train_data) - 1 > length(train_idx) / 2) {
      # 简单的特征选择：选择与目标变量相关性最高的特征
      correlations <- abs(cor(train_data[, -ncol(train_data)], train_data$target, use = "complete.obs"))
      top_features <- order(correlations, decreasing = TRUE)[1:min(floor(length(train_idx) / 2), ncol(train_data) - 1)]
      selected_cols <- c(colnames(train_data)[top_features], "target")
      train_data <- train_data[, selected_cols]
      test_data <- test_data[, colnames(train_data)[-ncol(train_data)]]
    }
    
    tryCatch({
      # 训练模型
      model <- gbm(
        target ~ .,
        data = train_data,
        distribution = "gaussian",
        n.trees = gbdt_params$n.trees,
        interaction.depth = gbdt_params$interaction.depth,
        shrinkage = gbdt_params$shrinkage,
        n.minobsinnode = gbdt_params$n.minobsinnode,
        bag.fraction = min(gbdt_params$bag.fraction, (length(train_idx) - 1) / length(train_idx)),
        train.fraction = 1.0,
        verbose = FALSE
      )
      
      # 预测
      all_predictions[test_idx] <- predict(model, test_data, n.trees = gbdt_params$n.trees)
    }, error = function(e) {
      warning(paste("Error in fold", fold, ":", e$message))
      all_predictions[test_idx] <- mean(y[train_idx])  # 使用均值预测
    })
  }
  
  return(calculate_r2(y, all_predictions))
}

# R²计算函数
calculate_r2 <- function(actual, predicted) {
  1 - sum((actual - predicted)^2) / sum((actual - mean(actual))^2)
}

# 特征选择函数
select_relevant_features <- function(X, y, correlation_method = "spearman", 
                                     p_threshold = 0.05,
                                     p_adjust_method = "none",
                                     rho_threshold = 0.1) {
  # 计算每个特征与目标变量的相关性
  correlations <- sapply(1:ncol(X), function(i) {
    result <- cor.test(X[,i], y, method = correlation_method)
    c(correlation = result$estimate,
      p_value = result$p.value)
  })
  
  # 转换为矩阵便于处理
  correlations<-as.data.frame(t(correlations))
  colnames(correlations) <- c("correlation", "p_value")
  rownames(correlations) <- colnames(X)
  
  # 进行多重检验校正
  adjusted_p_values <- p.adjust(correlations[,"p_value"], method = p_adjust_method)
  
  # 添加校正后的p值到结果中
  correlations <- cbind(correlations, adjusted_p_value = adjusted_p_values)
  
  # 同时根据校正后的p值和相关系数绝对值进行筛选
  significant_features <- which(adjusted_p_values < p_threshold & 
                                  abs(correlations[,"correlation"]) >= rho_threshold)
  
  # 按相关性绝对值排序
  if(length(significant_features) > 0) {
    abs_cors <- abs(correlations[significant_features, "correlation"])
    significant_features <- significant_features[order(abs_cors, decreasing = TRUE)]
  }
  
  return(list(
    selected_indices = significant_features,
    selected_features = colnames(X)[significant_features],
    correlations = correlations
  ))
}

# 新增函数：按分组预处理样本信息
preprocess_sample_groups <- function(sample_info, gut_data, oral_data, metabolite_data, 
                                     group_column = "IRIS", 
                                     group1 = "IR", group2 = "IS") {
  # 检查分组信息
  if(!group_column %in% colnames(sample_info)) {
    stop(paste("Group column", group_column, "not found in sample_info"))
  }
  
  # 获取各数据集的样本ID
  gut_samples <- colnames(gut_data)
  oral_samples <- colnames(oral_data)
  meta_samples <- colnames(metabolite_data)
  
  # 找出共同样本
  common_samples <- Reduce(intersect, list(gut_samples, oral_samples, meta_samples))
  
  # 匹配样本信息
  # 假设sample_info的行名或某一列包含样本ID
  if("sample_id" %in% colnames(sample_info)) {
    sample_ids <- sample_info$sample_id
  } else {
    sample_ids <- rownames(sample_info)
  }
  
  # 找出有分组信息且在所有数据集中都存在的样本
  valid_samples <- intersect(common_samples, sample_ids)
  
  # 提取这些样本的分组信息
  sample_groups <- sample_info[sample_ids %in% valid_samples, group_column]
  names(sample_groups) <- sample_ids[sample_ids %in% valid_samples]
  
  # 按分组筛选样本
  group1_samples <- names(sample_groups)[sample_groups == group1]
  group2_samples <- names(sample_groups)[sample_groups == group2]
  
  message(sprintf("Total valid samples: %d", length(valid_samples)))
  message(sprintf("%s group samples: %d", group1, length(group1_samples)))
  message(sprintf("%s group samples: %d", group2, length(group2_samples)))
  
  return(list(
    group1_samples = group1_samples,
    group2_samples = group2_samples,
    all_valid_samples = valid_samples,
    group_info = sample_groups
  ))
}

# 创建肠道菌群和口腔菌群的交互特征（按分组）
create_microbiome_interactions_by_group <- function(gut_data, oral_data, 
                                                    target_samples,
                                                    method = "multiplication",
                                                    feature_selection = TRUE,
                                                    max_features = 50,
                                                    correlation_threshold = 0.2) {
  # 提取目标样本的数据
  gut_matched <- gut_data[, target_samples]
  oral_matched <- oral_data[, target_samples]
  
  # 转置数据，使行为样本，列为特征
  gut_df <- as.data.frame(t(gut_matched))
  oral_df <- as.data.frame(t(oral_matched))
  
  # 添加前缀
  colnames(gut_df) <- paste0("gut_", colnames(gut_df))
  colnames(oral_df) <- paste0("oral_", colnames(oral_df))
  
  # 如果需要特征选择，则选择最重要的特征
  if (feature_selection) {
    # 对于肠道特征，基于方差和平均丰度选择
    gut_variances <- apply(gut_df, 2, var)
    gut_means <- apply(gut_df, 2, mean)
    gut_importance <- gut_variances * gut_means
    top_gut_indices <- order(gut_importance, decreasing = TRUE)[1:min(max_features, ncol(gut_df))]
    selected_gut_features <- colnames(gut_df)[top_gut_indices]
    
    # 对于口腔特征，同样基于方差和平均丰度选择
    oral_variances <- apply(oral_df, 2, var)
    oral_means <- apply(oral_df, 2, mean)
    oral_importance <- oral_variances * oral_means
    top_oral_indices <- order(oral_importance, decreasing = TRUE)[1:min(max_features, ncol(oral_df))]
    selected_oral_features <- colnames(oral_df)[top_oral_indices]
    
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
  interaction_list <- list()
  interaction_names <- character()
  
  # 计算特征之间的相关性，只选择相关的特征对创建交互
  if (correlation_threshold > 0) {
    combined_df <- cbind(gut_df_selected, oral_df_selected)
    correlation_matrix <- cor(combined_df, method = "spearman")
    
    gut_indices <- which(grepl("^gut_", colnames(combined_df)))
    oral_indices <- which(grepl("^oral_", colnames(combined_df)))
    
    interaction_count <- 0
    
    for (i in gut_indices) {
      gut_feature <- colnames(combined_df)[i]
      for (j in oral_indices) {
        oral_feature <- colnames(combined_df)[j]
        if (abs(correlation_matrix[i, j]) >= correlation_threshold) {
          if (method == "multiplication") {
            interaction_feature <- gut_df_selected[, gut_feature] * oral_df_selected[, oral_feature]
          } else if (method == "ratio") {
            denominator <- oral_df_selected[, oral_feature]
            denominator[denominator == 0] <- min(denominator[denominator > 0]) / 10
            interaction_feature <- gut_df_selected[, gut_feature] / denominator
          } else if (method == "difference") {
            interaction_feature <- gut_df_selected[, gut_feature] - oral_df_selected[, oral_feature]
          }
          
          interaction_name <- paste0("int_", gsub("gut_", "", gut_feature), "_", 
                                     gsub("oral_", "", oral_feature))
          interaction_list[[interaction_name]] <- interaction_feature
          interaction_names <- c(interaction_names, interaction_name)
          interaction_count <- interaction_count + 1
        }
      }
    }
  } else {
    for (gut_feature in selected_gut_features) {
      for (oral_feature in selected_oral_features) {
        if (method == "multiplication") {
          interaction_feature <- gut_df_selected[, gut_feature] * oral_df_selected[, oral_feature]
        } else if (method == "ratio") {
          denominator <- oral_df_selected[, oral_feature]
          denominator[denominator == 0] <- min(denominator[denominator > 0]) / 10
          interaction_feature <- gut_df_selected[, gut_feature] / denominator
        } else if (method == "difference") {
          interaction_feature <- gut_df_selected[, gut_feature] - oral_df_selected[, oral_feature]
        }
        
        interaction_name <- paste0("int_", gsub("gut_", "", gut_feature), "_", 
                                   gsub("oral_", "", oral_feature))
        interaction_list[[interaction_name]] <- interaction_feature
        interaction_names <- c(interaction_names, interaction_name)
      }
    }
  }
  
  # 将交互特征列表转换为数据框
  interaction_df <- as.data.frame(interaction_list)
  rownames(interaction_df) <- rownames(gut_df)
  
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

# 按分组预处理合并数据
preprocess_combined_data_by_group <- function(gut_data, oral_data, metabolite_data,
                                              target_samples,
                                              include_interactions = TRUE,
                                              interaction_method = "multiplication",
                                              max_features = 30,
                                              correlation_threshold = 0.3) {
  # 提取目标样本的数据
  gut_matched <- gut_data[, target_samples]
  oral_matched <- oral_data[, target_samples]
  metabolite_matched <- metabolite_data[, target_samples]
  
  # 转换和检查数据维度
  gut_df <- as.data.frame(t(gut_matched))
  oral_df <- as.data.frame(t(oral_matched))
  
  # 添加来源前缀
  colnames(gut_df) <- paste0("gut_", colnames(gut_df))
  colnames(oral_df) <- paste0("oral_", colnames(oral_df))
  
  # 基本的微生物组合并
  combined_microbiome <- cbind(gut_df, oral_df)
  
  # 添加交互特征（如果需要）
  if (include_interactions) {
    interactions_result <- create_microbiome_interactions_by_group(
      gut_data = gut_matched,
      oral_data = oral_matched,
      target_samples = target_samples,
      method = interaction_method,
      feature_selection = TRUE,
      max_features = max_features,
      correlation_threshold = correlation_threshold
    )
    
    interaction_df <- interactions_result$interaction_features
    
    # 确保行名匹配
    if (!identical(rownames(combined_microbiome), rownames(interaction_df))) {
      warning("Row names (samples) do not match between microbiome data and interaction features!")
    }
    
    # 合并所有特征
    combined_microbiome <- cbind(combined_microbiome, interaction_df)
  }
  
  # 确保样本顺序一致
  rownames(combined_microbiome) <- target_samples
  
  return(list(
    combined_microbiome = combined_microbiome,
    metabolite = t(metabolite_matched),
    target_samples = target_samples,
    interactions_info = if(include_interactions) interactions_result else NULL
  ))
}

# 分组比较分析主函数
analyze_group_differences_metabolite_interactions <- function(gut_data, oral_data, metabolite_data, 
                                                              sample_info,
                                                              group_column = "IRIS",
                                                              group1 = "IR", group2 = "IS",
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
  # 加载必要的库
  library(gbm)
  library(parallel)
  library(doParallel)
  library(foreach)
  library(progress)
  
  # 设置全局种子
  set.seed(seed)
  
  # 设置并行
  if(is.null(n_cores)) {
    n_cores <- detectCores() - 1
  }
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # 预处理分组信息
  group_info <- preprocess_sample_groups(sample_info, gut_data, oral_data, metabolite_data,
                                         group_column, group1, group2)
  
  # 分别处理两个分组的数据
  message("\nProcessing Group 1 data...")
  processed_data_group1 <- preprocess_combined_data_by_group(
    gut_data, oral_data, metabolite_data,
    target_samples = group_info$group1_samples,
    include_interactions = include_interactions,
    interaction_method = interaction_method,
    max_features = max_interaction_features,
    correlation_threshold = interaction_correlation_threshold
  )
  
  message("\nProcessing Group 2 data...")
  processed_data_group2 <- preprocess_combined_data_by_group(
    gut_data, oral_data, metabolite_data,
    target_samples = group_info$group2_samples,
    include_interactions = include_interactions,
    interaction_method = interaction_method,
    max_features = max_interaction_features,
    correlation_threshold = interaction_correlation_threshold
  )
  
  # 准备数据
  X1 <- as.matrix(processed_data_group1$combined_microbiome)
  Y1 <- processed_data_group1$metabolite
  X2 <- as.matrix(processed_data_group2$combined_microbiome)
  Y2 <- processed_data_group2$metabolite
  
  # 确保两组有相同的代谢物
  common_metabolites <- intersect(colnames(Y1), colnames(Y2))
  Y1 <- Y1[, common_metabolites, drop = FALSE]
  Y2 <- Y2[, common_metabolites, drop = FALSE]
  
  # 根据样本数量调整参数
  n_samples_group1 <- nrow(X1)
  n_samples_group2 <- nrow(X2)
  min_samples <- min(n_samples_group1, n_samples_group2)
  
  # 设置参数
  n_boots <- min(50, max(20, min_samples))  # 根据样本数调整bootstrap次数
  n_folds <- min(5, max(3, floor(min_samples / 5)))  # 调整交叉验证折数
  
  # 动态调整参数以适应小样本
  gbdt_params <- list(
    n.trees = min(30, max(50, min_samples * 2)),
    interaction.depth = min(10, max(3, floor(min_samples / 5))),
    shrinkage = 0.05,
    n.minobsinnode = max(2, min(5, floor(min_samples / 8))),
    bag.fraction = min(0.8, max(0.5, (min_samples - 5) / min_samples)),
    train.fraction = 1.0  # 使用全部数据进行训练，因为样本少
  )
  
  message(sprintf("Adjusted GBDT parameters for small sample sizes:"))
  message(sprintf("  n.trees: %d", gbdt_params$n.trees))
  message(sprintf("  interaction.depth: %d", gbdt_params$interaction.depth))
  message(sprintf("  n.minobsinnode: %d", gbdt_params$n.minobsinnode))
  message(sprintf("  bag.fraction: %.2f", gbdt_params$bag.fraction))
  
  # 导出函数到并行环境
  clusterExport(cl, c("single_cv", "calculate_r2", "select_relevant_features"), envir = environment())
  
  # 存储结果
  results_group1 <- list()
  results_group2 <- list()
  comparison_results <- list()
  
  # 设置进度条
  pb <- progress_bar$new(
    format = "[:bar] :percent | Metabolite :current/:total | Elapsed: :elapsed | ETA: :eta",
    total = length(common_metabolites)
  )
  
  message(sprintf("\nAnalyzing %d common metabolites...", length(common_metabolites)))
  message(sprintf("Group 1 (%s): %d samples, %d features", group1, nrow(X1), ncol(X1)))
  message(sprintf("Group 2 (%s): %d samples, %d features", group2, nrow(X2), ncol(X2)))
  message(sprintf("Bootstrap iterations: %d", n_boots))
  message(sprintf("Cross-validation folds: %d", n_folds))
  
  # 分析每个代谢物
  for(i in 1:length(common_metabolites)) {
    metabolite_name <- common_metabolites[i]
    current_y1 <- Y1[, metabolite_name]
    current_y2 <- Y2[, metabolite_name]
    
    # 分析组1
    if(do_feature_selection) {
      feature_selection1 <- select_relevant_features(
        X1, current_y1, correlation_method, p_threshold, p_adjust_method, rho_threshold
      )
      if (length(feature_selection1$selected_features) > 0) {
        X1_selected <- X1[, feature_selection1$selected_features, drop = FALSE]
      } else {
        X1_selected <- X1
      }
    } else {
      X1_selected <- X1
      feature_selection1 <- NULL
    }
    
    # 分析组2
    if(do_feature_selection) {
      feature_selection2 <- select_relevant_features(
        X2, current_y2, correlation_method, p_threshold, p_adjust_method, rho_threshold
      )
      if (length(feature_selection2$selected_features) > 0) {
        X2_selected <- X2[, feature_selection2$selected_features, drop = FALSE]
      } else {
        X2_selected <- X2
      }
    } else {
      X2_selected <- X2
      feature_selection2 <- NULL
    }
    
    # Bootstrap分析 - 组1
    if(ncol(X1_selected) > 0) {
      boot_results1 <- foreach(b = 1:n_boots,
                               .combine = 'c',
                               .packages = c("gbm", "caret")) %dopar% {
                                 local_seed <- seed * 1000 + b + i * n_boots
                                 set.seed(local_seed)
                                 boot_idx <- sample(1:nrow(X1_selected), nrow(X1_selected), replace = TRUE)
                                 boot_X <- X1_selected[boot_idx, , drop = FALSE]
                                 boot_y <- current_y1[boot_idx]
                                 
                                 single_cv(boot_X, boot_y, n_folds, gbdt_params, seed = local_seed)
                               }
      
      mean_r2_1 <- mean(boot_results1)
      ci_1 <- quantile(boot_results1, probs = c(0.025, 0.975))
      t_stat_1 <- mean_r2_1 / (sd(boot_results1) / sqrt(n_boots))
      p_value_1 <- 2 * pt(-abs(t_stat_1), df = n_boots - 1)
    } else {
      boot_results1 <- rep(0, n_boots)
      mean_r2_1 <- 0
      ci_1 <- c(0, 0)
      p_value_1 <- 1
    }
    
    # Bootstrap分析 - 组2
    if(ncol(X2_selected) > 0) {
      boot_results2 <- foreach(b = 1:n_boots,
                               .combine = 'c',
                               .packages = c("gbm", "caret")) %dopar% {
                                 local_seed <- seed * 2000 + b + i * n_boots
                                 set.seed(local_seed)
                                 boot_idx <- sample(1:nrow(X2_selected), nrow(X2_selected), replace = TRUE)
                                 boot_X <- X2_selected[boot_idx, , drop = FALSE]
                                 boot_y <- current_y2[boot_idx]
                                 
                                 single_cv(boot_X, boot_y, n_folds, gbdt_params, seed = local_seed)
                               }
      
      mean_r2_2 <- mean(boot_results2)
      ci_2 <- quantile(boot_results2, probs = c(0.025, 0.975))
      t_stat_2 <- mean_r2_2 / (sd(boot_results2) / sqrt(n_boots))
      p_value_2 <- 2 * pt(-abs(t_stat_2), df = n_boots - 1)
    } else {
      boot_results2 <- rep(0, n_boots)
      mean_r2_2 <- 0
      ci_2 <- c(0, 0)
      p_value_2 <- 1
    }
    
    # 组间差异检验
    diff_r2 <- mean_r2_1 - mean_r2_2
    pooled_sd <- sqrt((var(boot_results1) + var(boot_results2)) / 2)
    t_diff <- diff_r2 / (pooled_sd * sqrt(2/n_boots))
    p_diff <- 2 * pt(-abs(t_diff), df = 2*n_boots - 2)
    
    # 计算特征类型统计
    if(do_feature_selection) {
      # 组1特征统计
      selected_features1 <- if(!is.null(feature_selection1)) feature_selection1$selected_features else colnames(X1)
      n_gut1 <- sum(grepl("^gut_", selected_features1))
      n_oral1 <- sum(grepl("^oral_", selected_features1))
      n_interaction1 <- sum(grepl("^int_", selected_features1))
      
      # 组2特征统计
      selected_features2 <- if(!is.null(feature_selection2)) feature_selection2$selected_features else colnames(X2)
      n_gut2 <- sum(grepl("^gut_", selected_features2))
      n_oral2 <- sum(grepl("^oral_", selected_features2))
      n_interaction2 <- sum(grepl("^int_", selected_features2))
    } else {
      n_gut1 <- sum(grepl("^gut_", colnames(X1)))
      n_oral1 <- sum(grepl("^oral_", colnames(X1)))
      n_interaction1 <- sum(grepl("^int_", colnames(X1)))
      n_gut2 <- sum(grepl("^gut_", colnames(X2)))
      n_oral2 <- sum(grepl("^oral_", colnames(X2)))
      n_interaction2 <- sum(grepl("^int_", colnames(X2)))
    }
    
    # 存储结果
    results_group1[[i]] <- list(
      metabolite = metabolite_name,
      group = group1,
      mean_r2 = mean_r2_1,
      ci_lower = ci_1[1],
      ci_upper = ci_1[2],
      p_value = p_value_1,
      all_r2s = boot_results1,
      n_selected_features = length(if(do_feature_selection && !is.null(feature_selection1)) feature_selection1$selected_features else colnames(X1)),
      n_gut_features = n_gut1,
      n_oral_features = n_oral1,
      n_interaction_features = n_interaction1
    )
    
    results_group2[[i]] <- list(
      metabolite = metabolite_name,
      group = group2,
      mean_r2 = mean_r2_2,
      ci_lower = ci_2[1],
      ci_upper = ci_2[2],
      p_value = p_value_2,
      all_r2s = boot_results2,
      n_selected_features = length(if(do_feature_selection && !is.null(feature_selection2)) feature_selection2$selected_features else colnames(X2)),
      n_gut_features = n_gut2,
      n_oral_features = n_oral2,
      n_interaction_features = n_interaction2
    )
    
    comparison_results[[i]] <- list(
      metabolite = metabolite_name,
      r2_diff = diff_r2,
      r2_group1 = mean_r2_1,
      r2_group2 = mean_r2_2,
      p_diff = p_diff,
      effect_size = diff_r2 / pooled_sd,
      interaction_diff = n_interaction1 - n_interaction2,
      interaction_prop_diff = (n_interaction1/length(if(do_feature_selection && !is.null(feature_selection1)) feature_selection1$selected_features else colnames(X1))) - 
        (n_interaction2/length(if(do_feature_selection && !is.null(feature_selection2)) feature_selection2$selected_features else colnames(X2)))
    )
    
    # 更新进度条
    pb$tick()
    message(sprintf("\nMetabolite %d/%d (%s)", i, length(common_metabolites), metabolite_name))
    message(sprintf("%s: R² = %.3f (95%% CI: %.3f-%.3f)", group1, mean_r2_1, ci_1[1], ci_1[2]))
    message(sprintf("%s: R² = %.3f (95%% CI: %.3f-%.3f)", group2, mean_r2_2, ci_2[1], ci_2[2]))
    message(sprintf("Difference: %.3f (p = %.3e)", diff_r2, p_diff))
  }
  
  # 关闭并行集群
  stopCluster(cl)
  
  # 整理结果为数据框
  summary_df <- do.call(rbind, lapply(1:length(common_metabolites), function(i) {
    data.frame(
      metabolite = common_metabolites[i],
      group1_r2 = results_group1[[i]]$mean_r2,
      group1_ci_lower = results_group1[[i]]$ci_lower,
      group1_ci_upper = results_group1[[i]]$ci_upper,
      group1_p = results_group1[[i]]$p_value,
      group1_n_interaction = results_group1[[i]]$n_interaction_features,
      group2_r2 = results_group2[[i]]$mean_r2,
      group2_ci_lower = results_group2[[i]]$ci_lower,
      group2_ci_upper = results_group2[[i]]$ci_upper,
      group2_p = results_group2[[i]]$p_value,
      group2_n_interaction = results_group2[[i]]$n_interaction_features,
      r2_difference = comparison_results[[i]]$r2_diff,
      p_difference = comparison_results[[i]]$p_diff,
      effect_size = comparison_results[[i]]$effect_size,
      interaction_difference = comparison_results[[i]]$interaction_diff
    )
  }))
  
  # 多重检验校正
  summary_df$p_diff_adjusted <- p.adjust(summary_df$p_difference, method = "BH")
  summary_df$group1_p_adjusted <- p.adjust(summary_df$group1_p, method = "BH")
  summary_df$group2_p_adjusted <- p.adjust(summary_df$group2_p, method = "BH")
  
  return(list(
    group1_results = results_group1,
    group2_results = results_group2,
    comparison_results = comparison_results,
    summary = summary_df,
    group_info = group_info,
    group1_name = group1,
    group2_name = group2
  ))
}

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
    theme_minimal() +
    labs(title = paste("R² Comparison:", group1_name, "vs", group2_name),
         x = paste("R² in", group1_name, "group"),
         y = paste("R² in", group2_name, "group"),
         color = "Group Difference") +
    coord_equal()
  
  # 2. R²差异的分布图
  p2 <- ggplot(summary_df, aes(x = r2_difference)) +
    geom_histogram(bins = 30, alpha = 0.7, fill = "steelblue") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    theme_minimal() +
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
  
  p3 <- ggplot(interaction_data, aes(x = group, y = n_interactions, fill = group)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    scale_fill_manual(values=c("IR"="#E69F00", "IS"="#0072B2")) +
    theme_minimal() +
    labs(title = "Number of Selected Interaction Features",
         x = "Group",
         y = "Number of Interaction Features",
         fill = "Group")+stat_compare_means(label = "p.format")
  
  # 4. 效应大小vs显著性
  p4 <- ggplot(summary_df, aes(x = effect_size, y = -log10(p_diff_adjusted))) +
    geom_point(aes(color = abs(effect_size) > 0.5 & p_diff_adjusted < 0.05), 
               alpha = 0.7, size = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray") +
    scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"),
                       labels = c("Not significant", "Significant & Large effect")) +
    theme_minimal() +
    labs(title = "Effect Size vs Significance",
         x = "Effect Size (Cohen's d)",
         y = "-log10(Adjusted p-value)",
         color = "Classification")
  
  # 5. 火山图样式的R²差异图
  p5 <- ggplot(summary_df, aes(x = r2_difference, y = -log10(p_diff_adjusted))) +
    geom_point(aes(color = p_diff_adjusted < 0.05 & abs(r2_difference) > 0.1), 
               alpha = 0.7, size = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "gray") +
    scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"),
                       labels = c("Not significant", "Significant difference")) +
    theme_minimal() +
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
    theme_minimal() +
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


sample_info$IRIS[21]<-"IR"
sample_info$IRIS[24]<-"IR"

# 使用示例
# 假设你有sample_info数据框，其中包含IRIS列
group_comparison_results <- analyze_group_differences_metabolite_interactions(
  gut_data = gut_temp_object@expression_data,
  oral_data = oral_temp_object@expression_data,
  metabolite_data = metabolomics_temp_object@expression_data,
  sample_info = sample_info,  # 你的样本信息表格
  group_column = "IRIS",      # 分组列名
  group1 = "IR",              # 第一组名称
  group2 = "IS",              # 第二组名称
  do_feature_selection = TRUE,
  correlation_method = "spearman",
  p_threshold = 0.05,
  p_adjust_method = "none",
  rho_threshold = 0.3,
  include_interactions = TRUE,
  interaction_method = "multiplication",
  max_interaction_features = 30,
  interaction_correlation_threshold = 0.3
)

saveRDS(group_comparison_results,"group_comparison_results")

# 生成可视化结果
visualization_results <- plot_group_comparison_results(group_comparison_results)

# 显示组合图表
print(visualization_results$combined)

# 生成详细报告
detailed_report <- generate_group_comparison_report(group_comparison_results, top_n = 15)

# 查看显著差异的代谢物
significant_metabolites <- detailed_report$significant_metabolites
print(head(significant_metabolites, 10))


significant_metabolites<-merge(significant_metabolites,metabolite_annotation,by.x="metabolite",by.y="variable_id")


#挑选需要展示的代谢物




metabolites_index<-c("HMDB0001870","HMDB0000637","HMDB0001856","HMDB0059655","HMDB0002466","HMDB0000707","HMDB0000258","HMDB0002820","HMDB0000078","HMDB0002320","HMDB0001859","HMDB0001870","HMDB0000159","HMDB0000517","HMDB0000715","HMDB0000158","HMDB0002466","HMDB0000707","HMDB0001190","HMDB0000039","HMDB0002212")


significant_metabolites<-subset(significant_metabolites,HMDB%in%metabolites_index)



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
    theme_minimal() +
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

p2 <- plot_metabolites_mirror_style(significant_metabolites, top_n = 25)
print(p2)



