# load data
rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")

library(tidyverse)
library(tidymass)
library(plyr)
library(microbiomedataset)
library(readxl)
###load("data)
load("3_data_analysis/gut_microbiome/data_preparation/object_cross_section")

gut_object<-object_cross_section

load("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section")

metabolomics_object<-object_cross_section
metabolite_annotation<-read_excel("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx")

setwd("1_code/gut_oral_microbiome/")





# 加载包
library(gbm)         # GBDT模型
library(caret)       # 交叉验证
library(dplyr)       # 数据处理
library(foreach)     # 并行计算
library(doParallel)  # 并行后端
library(progress)    # 进度条
library(ggplot2)     # 可视化
library(tidyr)       # 数据整理
library(patchwork)   # 组合图表

# 数据预处理函数
preprocess_data <- function(microbiome_data, metabolite_data) {
  # 获取样本ID
  micro_samples <- colnames(microbiome_data)
  meta_samples <- colnames(metabolite_data)
  
  # 检查样本ID
  message("Initial sample counts:")
  message(sprintf("Microbiome samples: %d", length(micro_samples)))
  message(sprintf("Metabolite samples: %d", length(meta_samples)))
  
  # 找出共同样本
  common_samples <- intersect(micro_samples, meta_samples)
  message(sprintf("Common samples: %d", length(common_samples)))
  
  # 提取共同样本的数据
  microbiome_matched <- microbiome_data[, common_samples]
  metabolite_matched <- metabolite_data[, common_samples]
  
  return(list(
    microbiome = microbiome_matched,
    metabolite = metabolite_matched,
    common_samples = common_samples
  ))
}

# R²计算函数
calculate_r2 <- function(actual, predicted) {
  1 - sum((actual - predicted)^2) / sum((actual - mean(actual))^2)
}

# 特征选择函数
select_relevant_features <- function(X, y, correlation_method = "spearman", 
                                     p_threshold = 0.05,
                                     p_adjust_method = "none",
                                     rho_threshold = 0.1) {  # 添加相关系数阈值参数
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

# 交叉验证函数
single_cv <- function(X, y, n_folds = 5, gbdt_params, seed = NULL) {
  # 设置随机种子
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
    
    # 准备数据
    train_data <- data.frame(X[train_idx, , drop = FALSE])
    test_data <- data.frame(X[test_idx, , drop = FALSE])
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

# 可视化函数
plot_feature_selection_results <- function(results) {
  # 1. 准备特征选择摘要数据
  feature_summary <- do.call(rbind, lapply(results$detailed_results, function(x) {
    if(!is.null(x$feature_selection)) {
      data.frame(
        metabolite = x$metabolite,
        n_features = length(x$feature_selection$selected_features),
        r2 = x$mean_r2
      )
    }
  }))
  
  # 2. 提取所有相关性数据
  all_correlations <- do.call(rbind, lapply(results$detailed_results, function(x) {
    if(!is.null(x$feature_selection)) {
      correlations <- x$feature_selection$correlations
      data.frame(
        metabolite = x$metabolite,
        feature = rownames(correlations),
        correlation = correlations[,"correlation"],
        p_value = correlations[,"p_value"],
        adjusted_p_value = correlations[,"adjusted_p_value"],
        significant = correlations[,"adjusted_p_value"] < 0.05
      )
    }
  }))
  
  # 3. 统计每个特征被选中的次数
  feature_frequency <- all_correlations %>%
    filter(significant) %>%
    count(feature) %>%
    arrange(desc(n))
  
  # 创建可视化
  
  # 1. 特征数量与R²的关系图
  p1 <- ggplot(feature_summary, aes(x = n_features, y = r2)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "loess", se = TRUE) +
    theme_minimal() +
    labs(title = "Number of Selected Features vs R²",
         x = "Number of Selected Features",
         y = "R² Score")
  
  # 2. 选中特征数量分布图
  p2 <- ggplot(feature_summary, aes(x = n_features)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "white") +
    theme_minimal() +
    labs(title = "Distribution of Selected Feature Counts",
         x = "Number of Selected Features",
         y = "Count")
  
  # 3. Top 20最常被选中的特征
  p3 <- feature_frequency %>%
    head(20) %>%
    ggplot(aes(x = reorder(feature, n), y = n)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    theme_minimal() +
    labs(title = "Top 20 Most Frequently Selected Features",
         x = "Feature",
         y = "Number of Metabolites")
  
  # 4. 相关性强度分布图
  p4 <- ggplot(all_correlations %>% filter(significant), 
               aes(x = correlation)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "white") +
    theme_minimal() +
    labs(title = "Distribution of Significant Correlations",
         x = "Correlation Coefficient",
         y = "Count")
  
  # 5. 热图展示top特征与top代谢物的关系
  top_metabolites <- results$summary %>%
    arrange(desc(r2_mean)) %>%
    head(15) %>%
    pull(metabolite)
  
  top_features <- feature_frequency %>%
    head(15) %>%
    pull(feature)
  
  heatmap_data <- all_correlations %>%
    filter(metabolite %in% top_metabolites,
           feature %in% top_features)
  
  p5 <- ggplot(heatmap_data, 
               aes(x = feature, y = metabolite, fill = correlation)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Correlation Heatmap: Top Features vs Top Metabolites",
         x = "Features", y = "Metabolites")
  
  # 使用patchwork组合图表
  combined_plots <- (p1 + p2) / (p3 + p4) / p5 +
    plot_layout(heights = c(1, 1, 1.2))
  
  return(list(
    feature_vs_r2 = p1,
    feature_dist = p2,
    top_features = p3,
    correlation_dist = p4,
    correlation_heatmap = p5,
    combined = combined_plots
  ))
}

# 主分析函数
analyze_metabolite_ev <- function(microbiome_data, metabolite_data, n_cores = NULL, seed = 42,
                                  do_feature_selection = TRUE,
                                  correlation_method = "spearman",
                                  p_threshold = 0.05,
                                  p_adjust_method = "BH",
                                  rho_threshold = 0.3) {
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
  
  # 打印信息
  message("\nAnalysis settings:")
  message(sprintf("Feature selection: %s", if(do_feature_selection) "Yes" else "No"))
  if(do_feature_selection) {
    message(sprintf("Correlation method: %s", correlation_method))
    message(sprintf("P-value threshold: %.3f", p_threshold))
    message(sprintf("P-value adjustment method: %s", p_adjust_method))
    message(sprintf("Correlation coefficient threshold: %.2f", rho_threshold))
  }
  message(sprintf("Initial features (microbes): %d", ncol(X)))
  message(sprintf("Samples: %d", nrow(X)))
  message(sprintf("Metabolites: %d", ncol(Y)))
  message(sprintf("Using %d cores", n_cores))
  
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
    } else {
      # 如果没有选定特征，设置默认值
      boot_results <- rep(0, n_boots)
      mean_r2 <- 0
      ci <- c(0, 0)
      p_value <- 1
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
      n_selected_features = if(do_feature_selection) length(selected_features) else ncol(X)
    )
    
    # 更新进度条和打印结果
    pb$tick()
    end_time <- Sys.time()
    time_taken <- difftime(end_time, start_time, units = "mins")
    message(sprintf("\nMetabolite %d/%d (%s)", i, ncol(Y), colnames(Y)[i]))
    message(sprintf("Selected features: %d", length(selected_features)))
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
      n_selected_features = x$n_selected_features
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
    visualization = viz_results
  ))
}



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
# 运行预测（自动使用并行计算）


results <- analyze_metabolite_ev(
  microbiome_data, 
  metabolite_data,
  do_feature_selection = TRUE,
  correlation_method = "spearman",
  p_threshold = 0.05,
  p_adjust_method = "none",
  rho_threshold = 0.1
)


gut_results<-results


###load("data)
load("3_data_analysis/oral_microbiome/data_preparation/object_cross_section")

oral_object<-object_cross_section


####only remain the genus level


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



oral_sample_info <-
  oral_object@sample_info

#######adjust BMI, sex, and IRIS, ethnicity
oral_expression_data <-
  lm_adjust(expression_data = oral_expression_data,
            sample_info = oral_sample_info,
            threads = 3)

oral_temp_object <- oral_object
oral_temp_object@expression_data <- oral_expression_data
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


microbiome_data<-oral_temp_object@expression_data
metabolite_data<-metabolomics_temp_object@expression_data

results <- analyze_metabolite_ev(
  microbiome_data, 
  metabolite_data,
  do_feature_selection = TRUE,
  correlation_method = "spearman",
  p_threshold = 0.05,
  p_adjust_method = "none",
  rho_threshold = 0.1
)


oral_results<-results



##### 合并 gut 和 oral 的GBDT数据

oral_results_summary<-oral_results$summary
gut_results_summary<-gut_results$summary

gut_oral_results_summary<-cbind(gut_results_summary[,c(1,2,6)],oral_results_summary[,c(2,6)])
gut_oral_results_summary$HMDB

colnames(gut_oral_results_summary)<-c("metabolite","gut_R2","gut_features","oral_R2","oral_features")

gut_oral_results_summary<-merge(gut_oral_results_summary,metabolite_annotation[,c("variable_id","HMDB.Name","HMDB.Source.Microbial")],by.x="metabolite",by.y="variable_id")


gut_oral_results_summary <- gut_oral_results_summary %>%
  mutate(
    group = case_when(
      gut_R2 > 0.1 & oral_R2 <0.1 ~ "gut",
      oral_R2 > 0.1 & gut_R2 <0.1 ~ "oral",
      gut_R2 > 0.1 & oral_R2 > 0.1 ~ "co-influence",
      TRUE ~ "none"  # 这个处理其他情况，比如两个都是0的情况
    )
  )
# Create a subset for labeled points
labeled_data <- subset(gut_oral_results_summary, gut_R2 > 0.2 & oral_R2 > 0.1)



ggplot(gut_oral_results_summary, aes(x = gut_R2, y = oral_R2, fill = group, shape = HMDB.Source.Microbial)) +
  geom_point(size = 2, color = "white") +
  scale_shape_manual(values = c(21, 22, 24)) +
  scale_fill_manual(values = c("#f5ccca", "#edd064", "grey70", "#a1d5b9")) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21, size = 4)), 
    shape = guide_legend(override.aes = list(fill = "grey50", size = 4))
  ) +
  theme_pubclean()+ylim(c(0,0.3))+geom_text_repel(
    data = labeled_data,
    aes(label = HMDB.Name),
    box.padding = 0.5,   # Adjust label distance from points
    max.overlaps = 10,   # Allowable label overlaps
    force = 10,          # Increase label dispersion force
    segment.color = "grey50"
  ) +theme(
    legend.position = "bottom",
    axis.text.x = element_text(colour="black", size=14),
    axis.text.y = element_text(size=14, face="plain"),
    axis.title.y = element_text(size=14, face="plain"),
    axis.title.x = element_text(size=14, face="plain"),
    plot.title = element_text(size=15, face="bold", hjust=0.5)
  ) 
