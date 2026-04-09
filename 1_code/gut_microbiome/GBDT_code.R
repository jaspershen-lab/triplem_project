# Translated comment.
library(gbm)         # translated comment
library(caret)       # translated comment
library(dplyr)       # translated comment
library(foreach)     # translated comment
library(doParallel)  # translated comment
library(progress)    # translated comment

# Translated comment.

preprocess_data <- function(microbiome_data, metabolite_data) {
  # Translated comment.
  micro_samples <- colnames(microbiome_data)
  meta_samples <- colnames(metabolite_data)
  
  # Translated comment.
  message("Initial sample counts:")
  message(sprintf("Microbiome samples: %d", length(micro_samples)))
  message(sprintf("Metabolite samples: %d", length(meta_samples)))
  
  # Translated comment.
  common_samples <- intersect(micro_samples, meta_samples)
  message(sprintf("Common samples: %d", length(common_samples)))
  
  # Translated comment.
  microbiome_matched <- microbiome_data[, common_samples]
  metabolite_matched <- metabolite_data[, common_samples]
  
  return(list(
    microbiome = microbiome_matched,
    metabolite = metabolite_matched,
    common_samples = common_samples
  ))
}

calculate_r2 <- function(actual, predicted) {
  1 - sum((actual - predicted)^2) / sum((actual - mean(actual))^2)
}
# Translated comment.
single_cv <- function(X, y, n_folds = 5, gbdt_params, seed = NULL) {
  # Translated comment.
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  # Translated comment.
  fold_ids <- sample(rep(1:n_folds, length.out = length(y)))
  all_predictions <- numeric(length(y))
  
  for(fold in 1:n_folds) {
    # Translated comment.
    test_idx <- which(fold_ids == fold)
    train_idx <- which(fold_ids != fold)
    
    # Translated comment.
    train_data <- data.frame(X[train_idx, , drop = FALSE])
    train_data$target <- y[train_idx]
    
    # Translated comment.
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
    
    # Translated comment.
    test_data <- data.frame(X[test_idx, , drop = FALSE])
    all_predictions[test_idx] <- predict(model, test_data, n.trees = gbdt_params$n.trees)
  }
  
  return(calculate_r2(y, all_predictions))
}

# Translated comment.
analyze_metabolite_ev <- function(microbiome_data, metabolite_data, n_cores = NULL, seed = 42) {
  # Translated comment.
  set.seed(seed)
  
  # Translated comment.
  n_boots <- 100  # translated comment
  n_folds <- 10    # translated comment
  
  # Translated comment.
  if(is.null(n_cores)) {
    n_cores <- detectCores() - 1
  }
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Translated comment.
  processed_data <- preprocess_data(microbiome_data, metabolite_data)
  
  # Translated comment.
  X <- t(processed_data$microbiome)  # translated comment
  Y <- t(processed_data$metabolite)  # translated comment
  
  # Translated comment.
  message("\nData dimensions:")
  message(sprintf("Features (microbes): %d", ncol(X)))
  message(sprintf("Samples: %d", nrow(X)))
  message(sprintf("Metabolites: %d", ncol(Y)))
  message(sprintf("Using %d cores", n_cores))
  message(sprintf("Random seed: %d", seed))
  
  # Translated comment.
  gbdt_params <- list(
    n.trees = 50,          
    interaction.depth = 10,
    shrinkage = 0.01,
    n.minobsinnode = 10,
    bag.fraction = 0.5,
    train.fraction = 0.8
  )
  
  # Translated comment.
  clusterExport(cl, c("single_cv", "calculate_r2"), envir = environment())
  
  # Translated comment.
  results <- list()
  
  # Translated comment.
  pb <- progress_bar$new(
    format = "[:bar] :percent | Metabolite :current/:total | Elapsed: :elapsed | ETA: :eta",
    total = ncol(Y)
  )
  
  # Translated comment.
  for(i in 1:ncol(Y)) {
    start_time <- Sys.time()
    
    # Translated comment.
    current_y <- Y[, i]
    
    # Translated comment.
    boot_r2s <- foreach(b = 1:n_boots,
                        .combine = 'c',
                        .packages = c("gbm", "caret")) %dopar% {
                          # Translated comment.
                          local_seed <- seed * 1000 + b + i * n_boots
                          
                          # Translated comment.
                          set.seed(local_seed)
                          boot_idx <- sample(1:nrow(X), nrow(X), replace = TRUE)
                          boot_X <- X[boot_idx, , drop = FALSE]
                          boot_y <- current_y[boot_idx]
                          
                          # Translated comment.
                          single_cv(boot_X, boot_y, n_folds, gbdt_params, seed = local_seed)
                        }
    
    # Translated comment.
    mean_r2 <- mean(boot_r2s)
    ci <- quantile(boot_r2s, probs = c(0.025, 0.975))
    t_stat <- mean_r2 / (sd(boot_r2s) / sqrt(n_boots))
    p_value <- 2 * pt(-abs(t_stat), df = n_boots - 1)
    
    # Translated comment.
    results[[i]] <- list(
      metabolite = colnames(Y)[i],
      mean_r2 = mean_r2,
      ci_lower = ci[1],
      ci_upper = ci[2],
      p_value = p_value,
      all_r2s = boot_r2s
    )
    
    # Translated comment.
    pb$tick()
    
    # Translated comment.
    end_time <- Sys.time()
    time_taken <- difftime(end_time, start_time, units = "mins")
    message(sprintf("\nMetabolite %d/%d (%s): R² = %.3f (95%% CI: %.3f-%.3f, p = %.3e) Time: %.2f mins", 
                    i, ncol(Y), colnames(Y)[i], mean_r2, ci[1], ci[2], p_value, time_taken))
  }
  
  # Translated comment.
  stopCluster(cl)
  
  # Translated comment.
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

# Translated comment.

# Translated comment.
plot_ev_results <- function(results) {
  library(ggplot2)
  
  # Translated comment.
  plot_data <- results$summary %>%
    arrange(desc(r2_mean))
  
  # Translated comment.
  p1 <- ggplot(plot_data, aes(x = r2_mean)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "white") +
    theme_minimal() +
    labs(title = "Distribution of Explained Variance (R²)",
         x = "R²",
         y = "Count")
  
  # Translated comment.
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
  
  # Translated comment.
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