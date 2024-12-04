# 加载包
library(gbm)         # GBDT模型
library(caret)       # 交叉验证
library(dplyr)       # 数据处理
library(foreach)     # 并行计算
library(doParallel)  # 并行后端
library(progress)    # 进度条

# preprocess_data和calculate_r2函数保持不变

# 修改single_cv函数，添加种子参数
single_cv <- function(X, y, n_folds = 5, gbdt_params, seed = NULL) {
  # 如果提供了种子，设置随机种子
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  # 创建fold索引
  fold_ids <- sample(rep(1:n_folds, length.out = length(y)))
  all_predictions <- numeric(length(y))
  
  for(fold in 1:n_folds) {
    # 划分训练集和测试集
    test_idx <- which(fold_ids == fold)
    train_idx <- which(fold_ids != fold)
    
    # 准备训练数据
    train_data <- data.frame(X[train_idx, , drop = FALSE])
    train_data$target <- y[train_idx]
    
    # 训练模型
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
    
    # 预测
    test_data <- data.frame(X[test_idx, , drop = FALSE])
    all_predictions[test_idx] <- predict(model, test_data, n.trees = gbdt_params$n.trees)
  }
  
  return(calculate_r2(y, all_predictions))
}

# 添加种子参数的analyze_metabolite_ev函数
analyze_metabolite_ev <- function(microbiome_data, metabolite_data, n_cores = NULL, seed = 42) {
  # 设置全局种子
  set.seed(seed)
  
  # 设置参数
  n_boots <- 100  # 设置为100次bootstrap
  n_folds <- 10    # 保持5折交叉验证
  
  # 设置并行
  if(is.null(n_cores)) {
    n_cores <- detectCores() - 1
  }
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # 数据预处理
  processed_data <- preprocess_data(microbiome_data, metabolite_data)
  
  # 转置数据
  X <- t(processed_data$microbiome)  # 样本 × 特征
  Y <- t(processed_data$metabolite)  # 样本 × 代谢物
  
  # 打印维度信息
  message("\nData dimensions:")
  message(sprintf("Features (microbes): %d", ncol(X)))
  message(sprintf("Samples: %d", nrow(X)))
  message(sprintf("Metabolites: %d", ncol(Y)))
  message(sprintf("Using %d cores", n_cores))
  message(sprintf("Random seed: %d", seed))
  
  # GBDT参数
  gbdt_params <- list(
    n.trees = 30,          
    interaction.depth = 10,
    shrinkage = 0.1,
    n.minobsinnode = 10,
    bag.fraction = 0.5,
    train.fraction = 0.8
  )
  
  # 导出所需函数到并行环境
  clusterExport(cl, c("single_cv", "calculate_r2"), envir = environment())
  
  # 存储结果
  results <- list()
  
  # 设置进度条
  pb <- progress_bar$new(
    format = "[:bar] :percent | Metabolite :current/:total | Elapsed: :elapsed | ETA: :eta",
    total = ncol(Y)
  )
  
  # 对每个代谢物进行分析
  for(i in 1:ncol(Y)) {
    start_time <- Sys.time()
    
    # 获取当前代谢物数据
    current_y <- Y[, i]
    
    # 并行Bootstrap迭代，确保每次迭代使用不同但可重复的种子
    boot_r2s <- foreach(b = 1:n_boots,
                        .combine = 'c',
                        .packages = c("gbm", "caret")) %dopar% {
                          # 为每次bootstrap设置唯一但可重复的种子
                          local_seed <- seed * 1000 + b + i * n_boots
                          
                          # 抽样训练集
                          set.seed(local_seed)
                          boot_idx <- sample(1:nrow(X), nrow(X), replace = TRUE)
                          boot_X <- X[boot_idx, , drop = FALSE]
                          boot_y <- current_y[boot_idx]
                          
                          # 执行交叉验证
                          single_cv(boot_X, boot_y, n_folds, gbdt_params, seed = local_seed)
                        }
    
    # 计算统计量
    mean_r2 <- mean(boot_r2s)
    ci <- quantile(boot_r2s, probs = c(0.025, 0.975))
    t_stat <- mean_r2 / (sd(boot_r2s) / sqrt(n_boots))
    p_value <- 2 * pt(-abs(t_stat), df = n_boots - 1)
    
    # 存储结果
    results[[i]] <- list(
      metabolite = colnames(Y)[i],
      mean_r2 = mean_r2,
      ci_lower = ci[1],
      ci_upper = ci[2],
      p_value = p_value,
      all_r2s = boot_r2s
    )
    
    # 更新进度条
    pb$tick()
    
    # 打印当前结果
    end_time <- Sys.time()
    time_taken <- difftime(end_time, start_time, units = "mins")
    message(sprintf("\nMetabolite %d/%d (%s): R² = %.3f (95%% CI: %.3f-%.3f, p = %.3e) Time: %.2f mins", 
                    i, ncol(Y), colnames(Y)[i], mean_r2, ci[1], ci[2], p_value, time_taken))
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
      p_value = x$p_value
    )
  }))
  
  return(list(
    detailed_results = results,
    summary = summary_df,
    common_samples = processed_data$common_samples
  ))
}

# plot_ev_results函数保持不变

# 结果可视化函数
plot_ev_results <- function(results) {
  library(ggplot2)
  
  # 准备数据
  plot_data <- results$summary %>%
    arrange(desc(r2_mean))
  
  # 1. 显著性代谢物分布图
  p1 <- ggplot(plot_data, aes(x = r2_mean)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "white") +
    theme_minimal() +
    labs(title = "Distribution of Explained Variance (R²)",
         x = "R²",
         y = "Count")
  
  # 2. Top 20代谢物的R²和置信区间
  p2 <- plot_data %>%
    head(20) %>%
    ggplot(aes(x = reorder(metabolite, r2_mean), y = r2_mean)) +
    geom_point() +
    geom_errorbar(aes(ymin = r2_ci_lower, ymax = r2_ci_upper), width = 0.2) +
    coord_flip() +
    theme_minimal() +
    labs(title = "Top 20 Most Predictable Metabolites",
         x = "Metabolite",
         y = "R² (with 95% CI)")
  
  # 3. 火山图
  p3 <- ggplot(plot_data, 
               aes(x = r2_mean, y = -log10(p_value))) +
    geom_point(alpha = 0.6) +
    theme_minimal() +
    labs(title = "R² vs Statistical Significance",
         x = "R²",
         y = "-log10(P-value)")
  
  return(list(
    r2_dist = p1,
    top_metabolites = p2,
    volcano = p3
  ))
}