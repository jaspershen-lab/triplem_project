rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")
library(ggbeeswarm)
setwd("1_code/4_site_merge/")
##### Merge gut  oral GBDTdata.
gut_GBDT_results<-readRDS("../../3_data_analysis/gut_microbiome/GBDT/cross_section/gut_GBDT_results")
oral_GBDT_results<-readRDS("../../3_data_analysis/oral_microbiome/GBDT/cross_section/oral_GBDT_results")
metabolite_annotation<-read_excel("../../3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx")

oral_results_summary<-oral_GBDT_results$summary
gut_results_summary<-gut_GBDT_results$summary

gut_oral_results_summary<-cbind(gut_results_summary[,c(1,2,6)],oral_results_summary[,c(2,6)])
gut_oral_results_summary$HMDB

colnames(gut_oral_results_summary)<-c("metabolite","gut_R2","gut_features","oral_R2","oral_features")

gut_oral_results_summary<-merge(gut_oral_results_summary,metabolite_annotation[,c("variable_id","HMDB.Name","HMDB.Source.Microbial")],by.x="metabolite",by.y="variable_id")


preprocess_combined_data <- function(gut_data, oral_data, metabolite_data) {
  # Getsample IDs.
  gut_samples <- colnames(gut_data)
  oral_samples <- colnames(oral_data)
  meta_samples <- colnames(metabolite_data)
  
  # Checksample IDs.
  message("Initial sample counts:")
  message(sprintf("Gut microbiome samples: %d", length(gut_samples)))
  message(sprintf("Oral microbiome samples: %d", length(oral_samples)))
  message(sprintf("Metabolite samples: %d", length(meta_samples)))
  
  # samples.
  common_samples <- Reduce(intersect, list(gut_samples, oral_samples, meta_samples))
  message(sprintf("Common samples across all datasets: %d", length(common_samples)))
  
  # Extractsamplesdata.
  gut_matched <- gut_data[, common_samples]
  oral_matched <- oral_data[, common_samples]
  metabolite_matched <- metabolite_data[, common_samples]
  
  # ConvertCheckdata.
  gut_df <- as.data.frame(t(gut_matched))
  oral_df <- as.data.frame(t(oral_matched))
  
  # Add.
  colnames(gut_df) <- paste0("gut_", colnames(gut_df))
  colnames(oral_df) <- paste0("oral_", colnames(oral_df))
  
  # Check.
  message(sprintf("Gut features: %d", ncol(gut_df)))
  message(sprintf("Oral features: %d", ncol(oral_df)))
  
  # Mergedata.
  combined_microbiome <- cbind(gut_df, oral_df)
  
  # Addfeatures.
  colnames(combined_microbiome)[1:(ncol(combined_microbiome)-1)] <- 
    paste0(combined_microbiome$source, "", 
           colnames(combined_microbiome)[1:(ncol(combined_microbiome)-1)])
  
  # Ensure sample order matches.
  rownames(combined_microbiome) <- common_samples
  
  return(list(
    combined_microbiome = combined_microbiome,
    metabolite = t(metabolite_matched),
    common_samples = common_samples
  ))
}

# main analysis functionUseMergedata.
analyze_combined_metabolite_ev <- function(gut_data, oral_data, metabolite_data, 
                                           n_cores = NULL, seed = 42,
                                           do_feature_selection = TRUE,
                                           correlation_method = "spearman",
                                           p_threshold = 0.05,
                                           p_adjust_method = "BH",
                                           rho_threshold = 0.1) {
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
  processed_data <- preprocess_combined_data(gut_data, oral_data, metabolite_data)
  
  # data.
  X <- as.matrix(processed_data$combined_microbiome)
  Y <- processed_data$metabolite
  

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
  message(sprintf("Samples: %d", nrow(X)))
  message(sprintf("Metabolites: %d", ncol(Y)))
  message(sprintf("Using %d cores", n_cores))
  
  # GBDT.
  gbdt_params <- list(
    n.trees = 100,
    interaction.depth = 15,
    shrinkage = 0.01,
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
    
    # Feature selection.
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
      
      # Summarizein featuresin oralgutproportion.
      n_gut <- sum(grepl("^gut_", selected_features))
      n_oral <- sum(grepl("^oral_", selected_features))
    } else {
      X_selected <- X
      n_gut <- sum(grepl("^gut_", colnames(X)))
      n_oral <- sum(grepl("^oral_", colnames(X)))
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
    } else {
      boot_results <- rep(0, n_boots)
      mean_r2 <- 0
      ci <- c(0, 0)
      p_value <- 1
    }
    
    # results.
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
      n_oral_features = n_oral
    )
    
    # results.
    pb$tick()
    end_time <- Sys.time()
    time_taken <- difftime(end_time, start_time, units = "mins")
    message(sprintf("\nMetabolite %d/%d (%s)", i, ncol(Y), colnames(Y)[i]))
    message(sprintf("Selected features: %d (Gut: %d, Oral: %d)", 
                    length(selected_features), n_gut, n_oral))
    message(sprintf("R² = %.3f (95%% CI: %.3f-%.3f, p = %.3e)", 
                    mean_r2, ci[1], ci[2], p_value))
    message(sprintf("Time: %.2f mins", time_taken))
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
      n_gut_features = x$n_gut_features,
      n_oral_features = x$n_oral_features
    )
  }))
  
  # AddAnalyze.
  if(do_feature_selection) {
    viz_results <- plot_combined_feature_selection_results(list(
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

plot_combined_feature_selection_results <- function(results) {
  # data.
  feature_summary <- do.call(rbind, lapply(results$detailed_results, function(x) {
    data.frame(
      metabolite = x$metabolite,
      n_features = x$n_selected_features,
      n_gut = x$n_gut_features,
      n_oral = x$n_oral_features,
      r2 = x$mean_r2
    )
  }))
  
  # 1. featuresand R-squared.
  p1 <- ggplot(feature_summary, aes(x = n_features, y = r2)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "loess", se = TRUE) +
    theme_minimal() +
    labs(title = "Number of Selected Features vs R²",
         x = "Number of Selected Features",
         y = "R² Score")
  
  # 2. oralvsgutfeaturesproportion.
  p2 <- ggplot(feature_summary, aes(x = n_gut / (n_gut + n_oral), y = r2)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "loess", se = TRUE) +
    theme_minimal() +
    labs(title = "Proportion of Gut Features vs R²",
         x = "Proportion of Gut Features",
         y = "R² Score")
  
  # 3. features.
  feature_source_data <- data.frame(
    Source = rep(c("Gut", "Oral"), nrow(feature_summary)),
    Count = c(feature_summary$n_gut, feature_summary$n_oral),
    Metabolite = rep(feature_summary$metabolite, 2)
  )
  
  p3 <- ggplot(feature_source_data, aes(x = Source, y = Count)) +
    geom_boxplot(fill = "steelblue", alpha = 0.7) +
    theme_minimal() +
    labs(title = "Distribution of Selected Features by Source",
         y = "Number of Selected Features")
  

  combined_plots <- (p1 + p2) / p3 +
    plot_layout(heights = c(1, 0.8))
  
  return(list(
    feature_vs_r2 = p1,
    source_proportion = p2,
    source_distribution = p3,
    combined = combined_plots
  ))
}







 
gut_oral_interaction<-readRDS("../../1_code/gut_oral_microbiome/combined_results_with_interactions")
 
 
 gut_oral_results_summary_co_influence<- gut_oral_results_summary
 merge_model<-gut_oral_interaction$summary[,c("metabolite","r2_mean")]
 gut_oral_results_summary_co_influence<-merge(gut_oral_results_summary_co_influence,merge_model,by="metabolite")
   
 gut_oral_results_summary_co_influence$R2_diff<-gut_oral_results_summary_co_influence$r2_mean-(gut_oral_results_summary_co_influence$gut_R2+gut_oral_results_summary_co_influence$oral_R2)
 
 
 
 gut_oral_results_summary_co_influence<-subset(gut_oral_results_summary_co_influence,group="co-influence")
 
 
 p1 <- ggplot(gut_oral_results_summary_co_influence, aes(x = R2_diff)) +
   geom_histogram(binwidth = 0.02, fill = "steelblue", color = "black", alpha = 0.7) +
   theme_minimal() +
   labs(title = "Distribution of R2_diff",
        x = "R2_diff",
        y = "Count") +
   theme(plot.title = element_text(hjust = 0.5))
 
 
 
 # Read the data.
 data <- gut_oral_results_summary_co_influence
 
 # Process the data.
 library(tidyr)
 library(dplyr)
 library(ggplot2)
 # Handle duplicated HMDB.Name entries.
 data$HMDB.Name <- make.unique(data$HMDB.Name, sep = "_")
 # Sort the data and select the top 30 metabolites.
 data_sorted <- data %>%
   arrange(desc(r2_mean)) %>%
   slice(1:30)
 
 # Create long-format data for the stacked plot.
 data_long <- data_sorted %>%
   dplyr::select(HMDB.Name, gut_R2, oral_R2, r2_mean) %>%
   gather(key = "source", value = "value", c(gut_R2, oral_R2))
 
 # Create a separate long-format table for r2_mean.
 r2_mean_long <- data_sorted %>%
   dplyr::select(HMDB.Name, r2_mean) %>%
   dplyr::mutate(source = "r2_mean") %>%
   dplyr::rename(value = r2_mean)
 
 # Create the factor level order.
 level_order <- data_sorted$HMDB.Name
 
 # Convert HMDB.Name to a factor and set the level order.
 data_long$HMDB.Name <- factor(data_long$HMDB.Name, levels = level_order)
 r2_mean_long$HMDB.Name <- factor(r2_mean_long$HMDB.Name, levels = level_order)
 
 # Create the stacked bar chart and r2_mean comparison plot.
 ggplot() +
   # Stacked gut_R2 and oral_R2.
   geom_col(data = data_long, 
            aes(x = HMDB.Name, y = value, fill = source),
            position = "stack",
            width = 0.4) +
   # Bars for r2_mean.
   geom_col(data = r2_mean_long,
            aes(x = HMDB.Name, y = value, fill = source),
            width = 0.4,
            position = position_nudge(x = 0.4)) +  # Shift the r2_mean bars to the right.
   scale_fill_manual(values = c("gut_R2" = "#edd064", 
                                "oral_R2" = "#a1d5b9",
                                "r2_mean" = "grey50"),
                     labels = c("Gut", "Oral", "Gut+Oral")) +
   theme_minimal() +
   theme(
     axis.text.x = element_text(angle = 45, hjust = 0.95, size = 8),
     plot.title = element_text(hjust = 0.5),
     legend.title = element_blank(),
     legend.position = "top",
     panel.grid.major.x = element_blank()
   ) +
   labs(
     x = "",
     y = "Explained variance (%)"
   ) +
   scale_y_continuous(limits = c(0, 0.65))
 
 # Summarize.
 data_sorted$sum_R2 <- data_sorted$gut_R2 + data_sorted$oral_R2
 print("Summary of gut_R2 + oral_R2 vs r2_mean:")
 data_sorted$difference <- data_sorted$sum_R2 - data_sorted$r2_mean
 summary(data_sorted$difference)
 
 #