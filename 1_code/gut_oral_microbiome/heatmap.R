  preprocess_combined_data <- function(gut_data, oral_data, metabolite_data) {
    # Translated comment.
    gut_samples <- colnames(gut_data)
    oral_samples <- colnames(oral_data)
    meta_samples <- colnames(metabolite_data)
    
    # Translated comment.
    message("Initial sample counts:")
    message(sprintf("Gut microbiome samples: %d", length(gut_samples)))
    message(sprintf("Oral microbiome samples: %d", length(oral_samples)))
    message(sprintf("Metabolite samples: %d", length(meta_samples)))
    
    # Translated comment.
    common_samples <- Reduce(intersect, list(gut_samples, oral_samples, meta_samples))
    message(sprintf("Common samples across all datasets: %d", length(common_samples)))
    
    # Translated comment.
    gut_matched <- gut_data[, common_samples]
    oral_matched <- oral_data[, common_samples]
    metabolite_matched <- metabolite_data[, common_samples]
    
    # Translated comment.
    gut_df <- as.data.frame(t(gut_matched))
    oral_df <- as.data.frame(t(oral_matched))
    
    # Translated comment.
    colnames(gut_df) <- paste0("gut_", colnames(gut_df))
    colnames(oral_df) <- paste0("oral_", colnames(oral_df))
    
    # Translated comment.
    message(sprintf("Gut features: %d", ncol(gut_df)))
    message(sprintf("Oral features: %d", ncol(oral_df)))
    
    # Translated comment.
    combined_microbiome <- cbind(gut_df, oral_df)
    
    # Translated comment.
    if ("source" %in% colnames(combined_microbiome)) {
      colnames(combined_microbiome)[1:(ncol(combined_microbiome)-1)] <- 
        paste0(combined_microbiome$source, "_", 
               colnames(combined_microbiome)[1:(ncol(combined_microbiome)-1)])
    }
    
    # Translated comment.
    rownames(combined_microbiome) <- common_samples
    
    return(list(
      combined_microbiome = combined_microbiome,
      metabolite = t(metabolite_matched),
      common_samples = common_samples
    ))
  }
  
  analyze_combined_metabolite_ev <- function(gut_data, oral_data, metabolite_data, 
                                             n_cores = NULL, seed = 42,
                                             do_feature_selection = TRUE,
                                             correlation_method = "spearman",
                                             p_threshold = 0.05,
                                             p_adjust_method = "BH",
                                             rho_threshold = 0.1) {
    # Translated comment.
    set.seed(seed)
    
    # Translated comment.
    n_boots <- 100
    n_folds <- 10
    
    # Translated comment.
    if(is.null(n_cores)) {
      n_cores <- detectCores() - 1
    }
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    
    # Translated comment.
    processed_data <- preprocess_combined_data(gut_data, oral_data, metabolite_data)
    
    # Translated comment.
    X <- as.matrix(processed_data$combined_microbiome)
    Y <- processed_data$metabolite
    
    # Translated comment.
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
    
    # Translated comment.
    gbdt_params <- list(
      n.trees = 100,
      interaction.depth = 15,
      shrinkage = 0.01,
      n.minobsinnode = 8,
      bag.fraction = 0.8,
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
      current_y <- Y[, i]
      
      # Translated comment.
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
        
        # Translated comment.
        n_gut <- sum(grepl("^gut_", selected_features))
        n_oral <- sum(grepl("^oral_", selected_features))
      } else {
        X_selected <- X
        n_gut <- sum(grepl("^gut_", colnames(X)))
        n_oral <- sum(grepl("^oral_", colnames(X)))
      }
      
      # Translated comment.
      if(ncol(X_selected) > 0) {
        # Translated comment.
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
        
        # Translated comment.
        mean_r2 <- mean(boot_results)
        ci <- quantile(boot_results, probs = c(0.025, 0.975))
        t_stat <- mean_r2 / (sd(boot_results) / sqrt(n_boots))
        p_value <- 2 * pt(-abs(t_stat), df = n_boots - 1)
        
        # Translated comment.
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
      
      # Translated comment.
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
        feature_importance = feature_importance  # translated comment
      )
      
      # Translated comment.
      pb$tick()
      end_time <- Sys.time()
      time_taken <- difftime(end_time, start_time, units = "mins")
      message(sprintf("\nMetabolite %d/%d (%s)", i, ncol(Y), colnames(Y)[i]))
      message(sprintf("Selected features: %d (Gut: %d, Oral: %d)", 
                      if(do_feature_selection) length(selected_features) else ncol(X),
                      n_gut, n_oral))
      message(sprintf("R² = %.3f (95%% CI: %.3f-%.3f, p = %.3e)", 
                      mean_r2, ci[1], ci[2], p_value))
      message(sprintf("Time: %.2f mins", time_taken))
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
        p_value = x$p_value,
        n_selected_features = x$n_selected_features,
        n_gut_features = x$n_gut_features,
        n_oral_features = x$n_oral_features
      )
    }))
    
    # Translated comment.
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
    # Translated comment.
    feature_summary <- do.call(rbind, lapply(results$detailed_results, function(x) {
      data.frame(
        metabolite = x$metabolite,
        n_features = x$n_selected_features,
        n_gut = x$n_gut_features,
        n_oral = x$n_oral_features,
        r2 = x$mean_r2
      )
    }))
    
    # Translated comment.
    p1 <- ggplot(feature_summary, aes(x = n_features, y = r2)) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "loess", se = TRUE) +
      theme_minimal() +
      labs(title = "Number of Selected Features vs R²",
           x = "Number of Selected Features",
           y = "R² Score")
    
    # Translated comment.
    p2 <- ggplot(feature_summary, aes(x = n_gut / (n_gut + n_oral), y = r2)) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "loess", se = TRUE) +
      theme_minimal() +
      labs(title = "Proportion of Gut Features vs R²",
           x = "Proportion of Gut Features",
           y = "R² Score")
    
    # Translated comment.
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
    
    # Translated comment.
    combined_plots <- (p1 + p2) / p3 +
      plot_layout(heights = c(1, 0.8))
    
    return(list(
      feature_vs_r2 = p1,
      source_proportion = p2,
      source_distribution = p3,
      combined = combined_plots
    ))
  }
  
  combined_results <- analyze_combined_metabolite_ev(
    gut_data = gut_temp_object@expression_data,
    oral_data = oral_temp_object@expression_data,
    metabolite_data = metabolomics_temp_object@expression_data,
    do_feature_selection = TRUE,
    correlation_method = "spearman",
    p_threshold = 0.05,
    p_adjust_method = "none",
    rho_threshold = 0.1
  )
  
  
   library(ggplot2)
   feature_imp_df <- combined_results$detailed_results[[1]]$feature_importance
   ggplot(feature_imp_df, aes(x = reorder(var, rel.inf), y = rel.inf)) +
     geom_bar(stat = "identity", fill = "steelblue") +
     coord_flip() +
     labs(title = paste("Feature Importance for", combined_results$detailed_results[[1]]$metabolite),
          x = "Bacterial Feature", y = "Importance Score") +
     theme_minimal()
   
   
   
   # Translated comment.
   get_feature_importance_df <- function(combined_results) {
     # Translated comment.
     importance_list <- list()
     
     # Translated comment.
     for (i in seq_along(combined_results$detailed_results)) {
       metabolite_name <- combined_results$detailed_results[[i]]$metabolite
       feature_importance <- combined_results$detailed_results[[i]]$feature_importance
       
       # Translated comment.
       if (is.null(feature_importance)) next
       
       # Translated comment.
       feature_importance$metabolite <- metabolite_name
       
       # Translated comment.
       importance_list[[i]] <- feature_importance
     }
     
     # Translated comment.
     importance_df <- do.call(rbind, importance_list)
     
     # Translated comment.
     colnames(importance_df) <- c("species", "importance", "metabolite")
     
     return(importance_df)
   }
   
   # Translated comment.
   feature_importance_df <- get_feature_importance_df(combined_results)
   
   # Translated comment.
   head(feature_importance_df)
   
   
   
   
   gut_oral_results_summary_co_influence<- gut_oral_results_summary
   merge_model<-combined_results$summary[,c("metabolite","r2_mean")]
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
   
   
   
   # Translated comment.
   data <- gut_oral_results_summary_co_influence
   
   # Translated comment.
   library(tidyr)
   library(dplyr)
   library(ggplot2)
   # Translated comment.
   data$HMDB.Name <- make.unique(data$HMDB.Name, sep = "_")
   # Translated comment.
   data_sorted <- data %>%
     arrange(desc(r2_mean)) %>%
     slice(1:30)
   
   # Translated comment.
   data_long <- data_sorted %>%
     select(HMDB.Name, gut_R2, oral_R2, r2_mean) %>%
     gather(key = "source", value = "value", c(gut_R2, oral_R2))
   
   # Translated comment.
   r2_mean_long <- data_sorted %>%
     select(HMDB.Name, r2_mean) %>%
     mutate(source = "r2_mean") %>%
     rename(value = r2_mean)
   
   # Translated comment.
   level_order <- data_sorted$HMDB.Name
   
   # Translated comment.
   data_long$HMDB.Name <- factor(data_long$HMDB.Name, levels = level_order)
   r2_mean_long$HMDB.Name <- factor(r2_mean_long$HMDB.Name, levels = level_order)
   
   # Translated comment.
   ggplot() +
     # Translated comment.
     geom_col(data = data_long, 
              aes(x = HMDB.Name, y = value, fill = source),
              position = "stack",
              width = 0.4) +
     # Translated comment.
     geom_col(data = r2_mean_long,
              aes(x = HMDB.Name, y = value, fill = source),
              width = 0.4,
              position = position_nudge(x = 0.4)) +  # translated comment
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
   
   # Translated comment.
   data_sorted$sum_R2 <- data_sorted$gut_R2 + data_sorted$oral_R2
   print("Summary of gut_R2 + oral_R2 vs r2_mean:")
   data_sorted$difference <- data_sorted$sum_R2 - data_sorted$r2_mean
   summary(data_sorted$difference)
   
   # Translated comment.
   
   
   
   
   
   # Translated comment.
   
   feature_importance_df_co_influence<-feature_importance_df
   
   # Translated comment.
   
   feature_importance_df_co_influence<-subset(feature_importance_df_co_influence,importance>1)
   
   
   feature_importance_df_co_influence<-merge(feature_importance_df_co_influence,metabolite_annotation[,c("variable_id","HMDB.Name","HMDB.Class")],by.x="metabolite",by.y="variable_id")
   # Translated comment.
   
   merge_tax<-data.frame(rbind(gut_temp_object@variable_info,oral_temp_object@variable_info))
   
   # Translated comment.
   feature_importance_df_co_influence <- separate(feature_importance_df_co_influence, "species", into = c("group", "variable_id"), sep = "_",extra = "merge")
  
   feature_importance_df_co_influence<-merge(feature_importance_df_co_influence,merge_tax[,c("Genus","variable_id")],by="variable_id",all=TRUE)
   
  
   
   
   feature_importance_df_co_influence<-subset(feature_importance_df_co_influence,metabolite%in%data_sorted$metabolite)
   
  
   
  # Translated comment.
  library(ComplexHeatmap)
  library(circlize)
  library(tidyr)
  library(dplyr)
  
  # Translated comment.
  # Translated comment.
  feature_importance_df_co_influence <- feature_importance_df_co_influence %>%
    mutate(
      group_genus = paste(group, Genus, sep = "_"),
      hmdb_metabolite = paste(HMDB.Name, metabolite, sep = "_")
    )
  
  # Translated comment.
  heatmap_matrix <- feature_importance_df_co_influence %>%
    select(hmdb_metabolite, group_genus, importance) %>%
    pivot_wider(names_from = group_genus,
                values_from = importance,
                values_fill = list(importance = 0)) %>%
    as.data.frame()
  
  # Translated comment.
  rownames(heatmap_matrix) <- heatmap_matrix$hmdb_metabolite
  heatmap_matrix$hmdb_metabolite <- NULL
  
  # Translated comment.
  matrix_data <- as.matrix(heatmap_matrix)
  
  
  cols_to_keep <- apply(matrix_data, 2, function(x) sum(x > 0) > 3)
  
  # Translated comment.
  matrix_data <- matrix_data[, cols_to_keep]
  # Translated comment.
  
  
  # Translated comment.
  # Translated comment.
  scale_by_row <- function(x) {
    # Translated comment.
    non_zero <- x != 0
    if(sum(non_zero) > 0) {  # translated comment
      # Translated comment.
      x[non_zero] <- (x[non_zero] - min(x[non_zero])) / (max(x[non_zero]) - min(x[non_zero]))
    }
    return(x)
  }
  
  # Translated comment.
  scaled_matrix <- t(apply(matrix_data, 1, scale_by_row))
  scaled_matrix[27,13]<-0
  # Translated comment.
  col_fun = colorRamp2(c(0, 0.5, 1), 
                       c("# translated comment
  
  # Translated comment.
  cellwidth = 0.7
  cellheight = 0.7
  cn = dim(scaled_matrix)[2]
  rn = dim(scaled_matrix)[1]
  w = cellwidth * cn
  h = cellheight * rn
  
  # Translated comment.
  bacteria_counts <- colSums(scaled_matrix > 0)
  metabolite_counts <- rowSums(scaled_matrix > 0)
  
  # Translated comment.
  column_ha = HeatmapAnnotation(
    "Non-zero metabolites" = anno_barplot(bacteria_counts,
                                          height = unit(2, "cm"),
                                          gp = gpar(fill = "#008bd0", col = "grey"),
                                          border = TRUE),
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(fontsize = 8)
  )
  
  # Translated comment.
  row_ha = rowAnnotation(
    "Non-zero bacteria" = anno_barplot(metabolite_counts,
                                       width = unit(2, "cm"),
                                       gp = gpar(fill = "#ffa61d", col = "grey"),
                                       border = TRUE),
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(fontsize = 8)
  )
  
  # Translated comment.
  heatmap <- Heatmap(scaled_matrix,
                     name = "Scaled Importance",
                     col = col_fun,
                     # Translated comment.
                     width = unit(w, "cm"),
                     height = unit(h, "cm"),
                     rect_gp = gpar(col = "white", lwd = 2),
                     
                     # Translated comment.
                     border = TRUE,
                     border_gp = gpar(col = "grey", lwd = 2),
                     
                     # Translated comment.
                     cluster_rows = TRUE,
                     cluster_columns = TRUE,
                     show_row_dend = FALSE,
                     show_column_dend = FALSE,
                     
                     # Translated comment.
                     top_annotation = column_ha,
                     right_annotation = row_ha,
                     
                     # Translated comment.
                     show_row_names = TRUE,
                     show_column_names = TRUE,
                     row_names_gp = gpar(fontsize = 8),
                     column_names_gp = gpar(fontsize = 8),
                     column_names_rot = 45,
                     
                     # Translated comment.
                     row_title = "HMDB_Metabolites",
                     column_title = "Group_Genus",
                     column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                     row_title_gp = gpar(fontsize = 12, fontface = "bold"),
                     
                     # Translated comment.
                     heatmap_legend_param = list(
                       title = "Scaled Importance",
                       title_gp = gpar(fontsize = 10, fontface = "bold"),
                       labels_gp = gpar(fontsize = 8)
                     ))
  
  draw(heatmap)