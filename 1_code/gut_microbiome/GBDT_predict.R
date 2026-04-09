library(parallel)
# Add.

# Load packages.
library(gbm)         # GBDT.
library(caret)
library(dplyr)       # Process the data.
library(foreach)     # ParallelCalculate.
library(doParallel)  # Parallelafter.
library(progress)
library(ggplot2)
library(tidyr)       # Data.
library(patchwork)

# Data preprocessing function.
preprocess_data <- function(microbiome_data, metabolite_data) {
  # Getsample IDs.
  micro_samples <- colnames(microbiome_data)
  meta_samples <- colnames(metabolite_data)
  
  # Checksample IDs.
  message("Initial sample counts:")
  message(sprintf("Microbiome samples: %d", length(micro_samples)))
  message(sprintf("Metabolite samples: %d", length(meta_samples)))
  
  # samples.
  common_samples <- intersect(micro_samples, meta_samples)
  message(sprintf("Common samples: %d", length(common_samples)))
  
  # Extractsamplesdata.
  microbiome_matched <- microbiome_data[, common_samples]
  metabolite_matched <- metabolite_data[, common_samples]
  
  return(list(
    microbiome = microbiome_matched,
    metabolite = metabolite_matched,
    common_samples = common_samples
  ))
}

# R-squared helper function.
calculate_r2 <- function(actual, predicted) {
  1 - sum((actual - predicted)^2) / sum((actual - mean(actual))^2)
}

# Feature selection function.
select_relevant_features <- function(X, y, correlation_method = "spearman", 
                                     p_threshold = 0.05,
                                     p_adjust_method = "none",
                                     rho_threshold = 0.1) {  # Addcorrelation coefficientthreshold.
  # Calculatefeaturesand variablescorrelation.
  correlations <- sapply(1:ncol(X), function(i) {
    result <- cor.test(X[,i], y, method = correlation_method)
    c(correlation = result$estimate,
      p_value = result$p.value)
  })
  
  # ConvertmatrixProcess.
  correlations<-as.data.frame(t(correlations))
  colnames(correlations) <- c("correlation", "p_value")
  rownames(correlations) <- colnames(X)
  

  adjusted_p_values <- p.adjust(correlations[,"p_value"], method = p_adjust_method)
  
  # Addadjusted p-valuesresultsin.
  correlations <- cbind(correlations, adjusted_p_value = adjusted_p_values)
  
  # adjusted p-valuescorrelation coefficientfor Filter.
  significant_features <- which(adjusted_p_values < p_threshold & 
                                  abs(correlations[,"correlation"]) >= rho_threshold)
  
  # correlationfor sort.
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

# Cross-validation function.
single_cv <- function(X, y, n_folds = 5, gbdt_params, seed = NULL) {
  # Setrandom seed.
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  # Createfold.
  fold_ids <- sample(rep(1:n_folds, length.out = length(y)))
  all_predictions <- numeric(length(y))
  
  for(fold in 1:n_folds) {

    test_idx <- which(fold_ids == fold)
    train_idx <- which(fold_ids != fold)
    
    # data.
    train_data <- data.frame(X[train_idx, , drop = FALSE])
    test_data <- data.frame(X[test_idx, , drop = FALSE])
    train_data$target <- y[train_idx]
    

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
    

    test_data <- data.frame(X[test_idx, , drop = FALSE])
    all_predictions[test_idx] <- predict(model, test_data, n.trees = gbdt_params$n.trees)
  }
  
  return(calculate_r2(y, all_predictions))
}

plot_feature_selection_results <- function(results) {
  # 1. feature selectiondata.
  feature_summary <- do.call(rbind, lapply(results$detailed_results, function(x) {
    if(!is.null(x$feature_selection)) {
      data.frame(
        metabolite = x$metabolite,
        n_features = length(x$feature_selection$selected_features),
        r2 = x$mean_r2
      )
    }
  }))
  
  # 2. Extractcorrelationdata.
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
  
  # 3. Summarizefeaturesin .
  feature_frequency <- all_correlations %>%
    filter(significant) %>%
    count(feature) %>%
    arrange(desc(n))
  
  # Create.
  
  # 1. featuresand R-squared.
  p1 <- ggplot(feature_summary, aes(x = n_features, y = r2)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "loess", se = TRUE) +
    theme_minimal() +
    labs(title = "Number of Selected Features vs R²",
         x = "Number of Selected Features",
         y = "R² Score")
  
  # 2. in features.
  p2 <- ggplot(feature_summary, aes(x = n_features)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "white") +
    theme_minimal() +
    labs(title = "Distribution of Selected Feature Counts",
         x = "Number of Selected Features",
         y = "Count")
  
  # 3. Top 20in features.
  p3 <- feature_frequency %>%
    head(20) %>%
    ggplot(aes(x = reorder(feature, n), y = n)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    theme_minimal() +
    labs(title = "Top 20 Most Frequently Selected Features",
         x = "Feature",
         y = "Number of Metabolites")
  
  # 4. correlation.
  p4 <- ggplot(all_correlations %>% filter(significant), 
               aes(x = correlation)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "white") +
    theme_minimal() +
    labs(title = "Distribution of Significant Correlations",
         x = "Correlation Coefficient",
         y = "Count")
  
  # 5. heatmaptopfeaturesand topmetabolites.
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
  
  # Usepatchwork.
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

train_final_model <- function(X, y, selected_features = NULL, gbdt_params) {
  # If features,Usefeatures.
  if (!is.null(selected_features)) {
    X <- X[, selected_features, drop = FALSE]
  }
  
  # data.
  train_data <- data.frame(X)
  train_data$target <- y
  

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

# main analysis function.
analyze_metabolite_ev <- function(microbiome_data, metabolite_data, 
                                  model_save_dir = "models",
                                  n_cores = NULL, 
                                  seed = 42,
                                  do_feature_selection = TRUE,
                                  correlation_method = "spearman",
                                  p_threshold = 0.05,
                                  p_adjust_method = "BH",
                                  rho_threshold = 0.3) {
  # CreateSave.
  if (!dir.exists(model_save_dir)) {
    dir.create(model_save_dir, recursive = TRUE)
  }
  
  # Set.
  set.seed(seed)
  
  # Set.
  n_boots <- 100
  n_folds <- 10
  
  # Setparallel.
  if(is.null(n_cores)) {
    n_cores <- detectCores() - 1
  }
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Data preprocessing.
    processed_data <- preprocess_data(microbiome_data, metabolite_data)
  
  # Transposedata.
  X <- t(processed_data$microbiome)
  Y <- t(processed_data$metabolite)
  
  # GBDT.
  gbdt_params <- list(
    n.trees = 50,
    interaction.depth = 10,
    shrinkage = 0.03,
    n.minobsinnode = 8,
    bag.fraction = 0.8,
    train.fraction = 0.8
  )
  
  # parallel.
  clusterExport(cl, c("single_cv", "calculate_r2"), envir = environment())
  
  # results.
  results <- list()
  
  # Set.
  pb <- progress_bar$new(
    format = "[:bar] :percent | Metabolite :current/:total | Elapsed: :elapsed | ETA: :eta",
    total = ncol(Y)
  )
  
  # Analyzemetabolites.
  for(i in 1:ncol(Y)) {
    start_time <- Sys.time()
    current_y <- Y[, i]
    metabolite_name <- colnames(Y)[i]
    
    # samplesfeature selection.
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
    
    # featuresAnalyze.
    if(ncol(X_selected) > 0) {
      # ParallelBootstrap.
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
      
      # CalculateSummarize.
      mean_r2 <- mean(boot_results)
      ci <- quantile(boot_results, probs = c(0.025, 0.975))
      t_stat <- mean_r2 / (sd(boot_results) / sqrt(n_boots))
      p_value <- 2 * pt(-abs(t_stat), df = n_boots - 1)
      

      final_model <- train_final_model(X, current_y, selected_features, gbdt_params)
      
      # Save.
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
      
      # Save.
      model_file <- file.path(model_save_dir, paste0(make.names(metabolite_name), "_model.rds"))
      saveRDS(model_info, model_file)
      
    } else {
      boot_results <- rep(0, n_boots)
      mean_r2 <- 0
      ci <- c(0, 0)
      p_value <- 1
    }
    
    # results.
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
    
    # results.
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
  
  # Stop the parallel cluster.
  stopCluster(cl)
  
  # Organize the results into a data frame.
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
  
  # Add visualization results.
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


