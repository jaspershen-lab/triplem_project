# load data
rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")

library(tidyverse)
library(tidymass)
library(plyr)
library(microbiomedataset)
###load("data)
load("3_data_analysis/skin_microbiome/data_preparation/object_cross_section")

gut_object<-object_cross_section

load("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section")

metabolomics_object<-object_cross_section

dir.create("3_data_analysis/skin_microbiome/Lasso/cross_section/",recursive = TRUE)

setwd("3_data_analysis/skin_microbiome/Lasso/cross_section/")




####only remain the genus level


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
  which(non_zero_per > 0.2)

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


microbiome_data<-gut_temp_object@expression_data
metabolite_data<-metabolomics_temp_object@expression_data

# 加载必要的包
library(glmnet)
library(caret)
library(dplyr)
library(parallel)

# 计算R²的辅助函数
calculate_r2 <- function(pred, actual) {
  # 确保输入为数值向量
  pred <- as.numeric(pred)
  actual <- as.numeric(actual)
  
  # 检查是否有NA值
  if(any(is.na(pred)) || any(is.na(actual))) {
    return(NA)
  }
  
  # 检查是否有足够的变异
  if(sd(actual) == 0 || sd(pred) == 0) {
    return(0)
  }
  
  # 计算R²
  cor_val <- cor(pred, actual)
  return(cor_val^2)
}

# 计算p值的辅助函数
calculate_pvalue <- function(observed_r2, null_r2s) {
  # 检查输入
  if(is.na(observed_r2)) return(NA)
  if(all(is.na(null_r2s))) return(NA)
  
  # 移除NA值
  null_r2s <- null_r2s[!is.na(null_r2s)]
  
  # 如果没有有效的null分布，返回NA
  if(length(null_r2s) == 0) return(NA)
  
  # 计算p值
  p_value <- mean(null_r2s >= observed_r2, na.rm = TRUE)
  
  # 如果p值为0，使用1除以（排列次数+1）作为上限
  if(p_value == 0) {
    p_value <- 1 / (length(null_r2s) + 1)
  }
  
  return(p_value)
}

# 嵌套交叉验证函数
nested_cv_lasso <- function(X, y, outer_folds = 5, inner_folds = 5, 
                            lambda_factor = 5, n_permutations = 1000) {
  
  # 创建外层折
  set.seed(123)
  outer_fold_indices <- createFolds(1:length(y), k = outer_folds)
  
  # 存储结果
  outer_results <- list()
  
  for(i in 1:outer_folds) {
    # 分割训练集和测试集
    test_indices <- outer_fold_indices[[i]]
    train_indices <- setdiff(1:length(y), test_indices)
    
    X_train <- X[train_indices, ]
    y_train <- y[train_indices]
    X_test <- X[test_indices, ]
    y_test <- y[test_indices]
    
    # 在训练集上进行内层交叉验证以选择lambda
    cv_fit <- cv.glmnet(X_train, y_train, alpha = 1, nfolds = inner_folds)
    
    # 使用更大的lambda
    selected_lambda <- cv_fit$lambda.min * lambda_factor
    
    # 在全部训练集上训练最终模型
    final_model <- glmnet(X_train, y_train, alpha = 1, lambda = selected_lambda)
    
    # 预测测试集
    test_pred <- predict(final_model, newx = X_test)
    test_r2 <- calculate_r2(test_pred, y_test)
    test_rmse <- sqrt(mean((as.numeric(test_pred) - y_test)^2))
    
    # 排列测试：随机打乱y值多次，计算随机R²分布
    perm_r2 <- numeric(n_permutations)
    for(p in 1:n_permutations) {
      # 只打乱训练集的y值
      y_perm <- y_train[sample(length(y_train))]
      
      tryCatch({
        perm_model <- glmnet(X_train, y_perm, alpha = 1, lambda = selected_lambda)
        perm_pred <- predict(perm_model, newx = X_test)
        perm_r2[p] <- calculate_r2(perm_pred, y_test)
      }, error = function(e) {
        perm_r2[p] <- NA
      })
    }
    
    # 计算p值
    p_value <- calculate_pvalue(test_r2, perm_r2)
    
    # 获取选择的特征
    coef_matrix <- coef(final_model)
    selected_features <- rownames(coef_matrix)[which(coef_matrix != 0)]
    
    # 存储结果
    outer_results[[i]] <- list(
      test_r2 = test_r2,
      test_rmse = test_rmse,
      n_features = length(selected_features) - 1,  # 减去截距项
      selected_features = selected_features[-1],   # 减去截距项
      lambda = selected_lambda,
      p_value = p_value,
      perm_r2_dist = perm_r2,
      n_valid_perms = sum(!is.na(perm_r2))  # 记录有效的排列次数
    )
  }
  
  return(outer_results)
}

# 主分析函数
analyze_all_metabolites <- function(metabolite_data, microbiome_data, 
                                    lambda_factor = 5,
                                    outer_folds = 5, 
                                    inner_folds = 5,
                                    n_permutations = 1000,
                                    ncores = 1) {
  
  # 数据预处理
  common_samples <- intersect(colnames(metabolite_data), colnames(microbiome_data))
  metabolite_subset <- metabolite_data[, common_samples]
  microbiome_subset <- microbiome_data[, common_samples]
  
  X <- t(microbiome_subset)
  y <- t(metabolite_subset)
  
  # 存储结果
  results_list <- list()
  
  # 并行处理设置
  if(ncores > 1) {
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
  }
  
  # 对每个代谢物进行分析
  for(i in 1:nrow(metabolite_data)) {
    tryCatch({
      cat(sprintf("\n分析代谢物 %d/%d: %s\n", 
                  i, nrow(metabolite_data), 
                  rownames(metabolite_data)[i]))
      
      # 运行嵌套交叉验证
      cv_results <- nested_cv_lasso(X, y[,i], 
                                    outer_folds = outer_folds,
                                    inner_folds = inner_folds,
                                    lambda_factor = lambda_factor,
                                    n_permutations = n_permutations)
      
      # 计算平均指标
      mean_r2 <- mean(sapply(cv_results, function(x) x$test_r2), na.rm = TRUE)
      mean_rmse <- mean(sapply(cv_results, function(x) x$test_rmse), na.rm = TRUE)
      mean_features <- mean(sapply(cv_results, function(x) x$n_features), na.rm = TRUE)
      valid_p_values <- sapply(cv_results, function(x) x$p_value)
      valid_p_values <- valid_p_values[!is.na(valid_p_values)]
      
      # 如果有有效的p值，取最小值；否则为NA
      min_p_value <- if(length(valid_p_values) > 0) min(valid_p_values) else NA
      
      # 计算有效排列测试的平均数量
      mean_valid_perms <- mean(sapply(cv_results, function(x) x$n_valid_perms))
      
      # 存储结果
      results_list[[i]] <- list(
        metabolite = rownames(metabolite_data)[i],
        mean_test_r2 = mean_r2,
        mean_test_rmse = mean_rmse,
        mean_n_features = mean_features,
        min_p_value = min_p_value,
        mean_valid_perms = mean_valid_perms,
        significant = !is.na(min_p_value) && min_p_value < 0.05,
        cv_results = cv_results
      )
      
      # 打印当前结果
      cat(sprintf("平均测试 R² = %.3f\n", mean_r2))
      cat(sprintf("平均特征数 = %.1f\n", mean_features))
      cat(sprintf("最小 p-value = %s\n", 
                  if(is.na(min_p_value)) "NA" else sprintf("%.3f", min_p_value)))
      cat(sprintf("平均有效排列次数 = %.1f\n", mean_valid_perms))
      
    }, error = function(e) {
      cat(sprintf("处理代谢物 %s 时出错: %s\n", 
                  rownames(metabolite_data)[i], 
                  conditionMessage(e)))
      
      # 记录错误情况
      results_list[[i]] <- list(
        metabolite = rownames(metabolite_data)[i],
        error = conditionMessage(e)
      )
    })
  }
  
  if(ncores > 1) {
    stopCluster(cl)
  }
  
  # 整理结果为数据框
  valid_results <- !sapply(results_list, function(x) exists("error", where = x))
  results_df <- data.frame(
    metabolite = sapply(results_list[valid_results], function(x) x$metabolite),
    test_r2 = sapply(results_list[valid_results], function(x) x$mean_test_r2),
    test_rmse = sapply(results_list[valid_results], function(x) x$mean_test_rmse),
    n_features = sapply(results_list[valid_results], function(x) x$mean_n_features),
    p_value = sapply(results_list[valid_results], function(x) x$min_p_value),
    valid_perms = sapply(results_list[valid_results], function(x) x$mean_valid_perms),
    significant = sapply(results_list[valid_results], function(x) x$significant)
  )
  
  # 应用FDR校正（只对有效的p值进行校正）
  valid_p <- !is.na(results_df$p_value)
  results_df$p_adj <- NA
  results_df$p_adj[valid_p] <- p.adjust(results_df$p_value[valid_p], method = "BH")
  results_df$significant_adj <- !is.na(results_df$p_adj) & results_df$p_adj < 0.05
  
  return(list(
    results_df = results_df,
    detailed_results = results_list
  ))
}


results <- analyze_all_metabolites(
  metabolite_data, 
  microbiome_data,
  lambda_factor = 1.5,
  n_permutations = 100,
  ncores = 10
)


saveRDS(results,"skin_lasso_results")