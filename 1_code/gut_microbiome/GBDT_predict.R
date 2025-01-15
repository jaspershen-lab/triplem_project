# 添加新的模型训练函数
train_final_model <- function(X, y, selected_features = NULL, gbdt_params) {
  # 如果有选定的特征，只使用这些特征
  if (!is.null(selected_features)) {
    X <- X[, selected_features, drop = FALSE]
  }
  
  # 准备训练数据
  train_data <- data.frame(X)
  train_data$target <- y
  
  # 训练最终模型
  model <- gbm(
    target ~ .,
    data = train_data,
    distribution = "gaussian",
    n.trees = gbdt_params$n.trees,
    interaction.depth = gbdt_params$interaction.depth,
    shrinkage = gbdt_params$shrinkage,
    n.minobsinnode = gbdt_params$n.minobsinnode,
    verbose = FALSE
  )
  
  return(model)
}

# 修改主分析函数
analyze_metabolite_ev <- function(microbiome_data, metabolite_data, 
                                  model_save_dir = "models",
                                  n_cores = NULL, 
                                  seed = 42,
                                  do_feature_selection = TRUE,
                                  correlation_method = "spearman",
                                  p_threshold = 0.05,
                                  p_adjust_method = "BH",
                                  rho_threshold = 0.3) {
  # 创建保存模型的目录
  if (!dir.exists(model_save_dir)) {
    dir.create(model_save_dir, recursive = TRUE)
  }
  
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
  
  # 数据预处理
  processed_data <- preprocess_data(microbiome_data, metabolite_data)
  
  # 转置数据
  X <- t(processed_data$microbiome)
  Y <- t(processed_data$metabolite)
  
  # GBDT参数
  gbdt_params <- list(
    n.trees = 50,
    interaction.depth = 10,
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
    metabolite_name <- colnames(Y)[i]
    
    # 在所有样本上进行特征选择
    if(do_feature_selection) {
      feature_selection <- select_relevant_features(
        X, current_y,
        correlation_method = correlation_method,
        p_threshold = p_threshold,
        p_adjust_method = p_adjust_method,
        rho_threshold = rho_threshold
      )
      selected_features <- feature_selection$selected_features
      X_selected <- X[, selected_features, drop = FALSE]
    } else {
      X_selected <- X
      selected_features <- colnames(X)
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
      
      # 训练最终模型
      final_model <- train_final_model(X, current_y, selected_features, gbdt_params)
      
      # 保存模型和相关信息
      model_info <- list(
        model = final_model,
        selected_features = selected_features,
        feature_importance = summary(final_model, plot = FALSE),
        performance = list(
          mean_r2 = mean_r2,
          ci_lower = ci[1],
          ci_upper = ci[2],
          p_value = p_value
        ),
        parameters = gbdt_params,
        feature_selection_params = list(
          method = correlation_method,
          p_threshold = p_threshold,
          p_adjust_method = p_adjust_method,
          rho_threshold = rho_threshold
        )
      )
      
      # 保存模型
      model_file <- file.path(model_save_dir, paste0(make.names(metabolite_name), "_model.rds"))
      saveRDS(model_info, model_file)
      
    } else {
      boot_results <- rep(0, n_boots)
      mean_r2 <- 0
      ci <- c(0, 0)
      p_value <- 1
    }
    
    # 存储结果
    results[[i]] <- list(
      metabolite = metabolite_name,
      mean_r2 = mean_r2,
      ci_lower = ci[1],
      ci_upper = ci[2],
      p_value = p_value,
      all_r2s = boot_results,
      feature_selection = if(do_feature_selection) feature_selection else NULL,
      n_selected_features = if(do_feature_selection) length(selected_features) else ncol(X),
      model_file = if(ncol(X_selected) > 0) model_file else NULL
    )
    
    # 更新进度条和打印结果
    pb$tick()
    end_time <- Sys.time()
    time_taken <- difftime(end_time, start_time, units = "mins")
    message(sprintf("\nMetabolite %d/%d (%s)", i, ncol(Y), metabolite_name))
    message(sprintf("Selected features: %d", length(selected_features)))
    message(sprintf("R² = %.3f (95%% CI: %.3f-%.3f, p = %.3e)", 
                    mean_r2, ci[1], ci[2], p_value))
    message(sprintf("Time: %.2f mins", time_taken))
    if(ncol(X_selected) > 0) {
      message(sprintf("Model saved to: %s", model_file))
    }
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
      model_file = x$model_file
    )
  }))
  
  # 添加可视化结果
  if(do_feature_selection) {
    viz_results <- plot_feature_selection_results(list(
      detailed_results = results,
      summary = summary_df
    ))
  } else {
    viz_results <- NULL
  }
  
  return(list(
    detailed_results = results,
    summary = summary_df,
    common_samples = processed_data$common_samples,
    visualization = viz_results,
    model_directory = model_save_dir
  ))
}


# 使用示例:
# results <- analyze_metabolite_ev(
#   microbiome_data, 
#   metabolite_data,
#   model_save_dir = "metabolite_models",
#   do_feature_selection = TRUE,
#   correlation_method = "spearman",
#   p_threshold = 0.05,
#   p_adjust_method = "none",
#   rho_threshold = 0.1
# )
#
