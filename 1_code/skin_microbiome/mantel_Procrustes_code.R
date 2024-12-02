# Load required packages
library(vegan)
library(ggplot2)
library(ggrepel)
library(tidyverse)

#' Select metabolites significantly associated with gut microbiome
#' @param gut_microbiome Preprocessed gut microbiome data (species × samples)
#' @param metabolome Preprocessed metabolome data (metabolites × samples)
#' @param cor_threshold Correlation threshold (absolute value)
#' @param p_threshold P-value threshold for significance
select_significant_metabolites <- function(gut_microbiome, metabolome, 
                                           cor_threshold = 0.2,
                                           p_threshold = 0.05) {
  # Match samples
  common_samples <- intersect(colnames(gut_microbiome), colnames(metabolome))
  gut_matched <- gut_microbiome[, common_samples]
  meta_matched <- metabolome[, common_samples]
  
  # Store results
  n_metabolites <- nrow(meta_matched)
  n_species <- nrow(gut_matched)
  
  results <- data.frame(
    metabolite = rownames(meta_matched),
    max_cor = numeric(n_metabolites),
    max_species = character(n_metabolites),
    n_sig_cors = numeric(n_metabolites)
  )
  
  # For each metabolite
  for(i in 1:n_metabolites) {
    # Calculate correlations with all species
    cors <- numeric(n_species)
    p_values <- numeric(n_species)
    
    for(j in 1:n_species) {
      cor_test <- cor.test(as.numeric(meta_matched[i,]), 
                           as.numeric(gut_matched[j,]), 
                           method = "spearman")
      cors[j] <- cor_test$estimate
      p_values[j] <- cor_test$p.value
    }
    
    # Count significant correlations
    sig_strong_cors <- abs(cors) >= cor_threshold & p_values < p_threshold
    results$n_sig_cors[i] <- sum(sig_strong_cors)
    
    # Record the strongest correlation
    max_idx <- which.max(abs(cors))
    results$max_cor[i] <- cors[max_idx]
    results$max_species[i] <- rownames(gut_matched)[max_idx]
  }
  
  # Select significant metabolites
  sig_idx <- which(results$n_sig_cors > 0)
  sig_metabolites <- results$metabolite[sig_idx]
  
  # Print summary
  cat("\n=== Metabolite Selection Results ===\n")
  cat(sprintf("Total metabolites analyzed: %d\n", n_metabolites))
  cat(sprintf("Significant metabolites found: %d\n", length(sig_metabolites)))
  cat(sprintf("Using correlation threshold: %.2f\n", cor_threshold))
  cat(sprintf("Using p-value threshold: %.3f\n", p_threshold))
  
  # Sort results by number of significant correlations
  results <- results[order(-results$n_sig_cors), ]
  
  # Print top metabolites
  cat("\nTop 10 metabolites by number of significant correlations:\n")
  print(head(results[, c("metabolite", "n_sig_cors", "max_cor")], 10))
  
  return(list(
    all_results = results,
    significant_metabolites = sig_metabolites,
    significant_indices = sig_idx
  ))
}

#' Plot metabolite selection results with compound names
#' @param selection_results Results from select_significant_metabolites function
#' @param metabolomics_class Data frame containing metabolite class information
plot_metabolite_selection <- function(selection_results, metabolomics_class) {
  # 创建代谢物ID和名称的对应关系
  id_to_name <- setNames(metabolomics_class$Compound.name, 
                         rownames(metabolomics_class))
  
  # 添加代谢物名称到结果数据框
  results_df <- selection_results$all_results
  results_df$compound_name <- id_to_name[results_df$metabolite]
  
  # Get number of selected metabolites
  n_selected <- sum(results_df$n_sig_cors > 0)
  
  # Create scatter plot
  p <- ggplot(results_df, 
              aes(x = abs(max_cor), 
                  y = n_sig_cors)) +
    geom_point(aes(color = n_sig_cors > 0)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
    geom_vline(xintercept = 0.3, linetype = "dashed", color = "red") +
    ggrepel::geom_text_repel(
      data = head(results_df[order(-results_df$n_sig_cors),], 10),
      aes(label = compound_name),  # 使用compound_name而不是metabolite
      size = 3,
      box.padding = 0.5,
      force = 10,
      max.overlaps = Inf
    ) +
    annotate("text",
             x = min(abs(results_df$max_cor)),
             y = max(results_df$n_sig_cors),
             label = sprintf("Selected metabolites: %d", n_selected),
             hjust = 0,
             vjust = 1,
             size = 4,
             fontface = "bold") +
    theme_minimal() +
    labs(title = "Metabolite Selection Results",
         x = "Maximum Absolute Correlation",
         y = "Number of Significant Correlations",
         color = "Selected") +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.3))) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))
  
  # Save plot
  ggsave("metabolite_selection.pdf", p, width = 12, height = 10)
  
  return(p)
}
analyze_associations <- function(gut_microbiome, metabolome, selected_metabolites, metadata) {
  # Extract significant metabolites
  meta_significant <- metabolome[selected_metabolites, ]
  
  # Match samples
  common_samples <- intersect(colnames(gut_microbiome), colnames(meta_significant))
  common_samples <- intersect(common_samples, rownames(metadata))
  
  gut_matched <- gut_microbiome[, common_samples]
  meta_matched <- meta_significant[, common_samples]
  meta_info <- metadata[common_samples, ]
  
  # Prepare data for co-inertia
  gut_df <- as.data.frame(t(gut_matched))
  meta_df <- as.data.frame(t(meta_matched))
  
  # 1. Overall Analysis
  ## 1.1 Mantel test
  gut_dist <- vegdist(t(gut_matched), method = "bray")
  meta_dist <- vegdist(t(meta_matched), method = "euclidean")
  mantel_result <- mantel(gut_dist, meta_dist, permutations = 999)
  
  ## 1.2 Procrustes analysis
  gut_pcoa <- cmdscale(gut_dist, k=2)
  meta_pcoa <- cmdscale(meta_dist, k=2)
  proc_result <- protest(gut_pcoa, meta_pcoa, permutations = 999)
  
  ## 1.3 Co-inertia analysis
  gut_pca <- dudi.pca(gut_df, scannf = FALSE, nf = 2)
  meta_pca <- dudi.pca(meta_df, scannf = FALSE, nf = 2)
  coin_result <- coinertia(gut_pca, meta_pca, scannf = FALSE, nf = 2)
  
  # 2. Individual Analysis
  ## Calculate individual co-inertia distances
  gut_coords <- coin_result$mX
  meta_coords <- coin_result$mY
  individual_distances <- sqrt(rowSums((gut_coords - meta_coords)^2))
  
  ## Create individual results data frame
  individual_results <- data.frame(
    sample = common_samples,
    IRIS = meta_info$IRIS,
    coinertia_distance = individual_distances,
    gut_x = gut_coords[,1],
    gut_y = gut_coords[,2],
    meta_x = meta_coords[,1],
    meta_y = meta_coords[,2]
  )
  
  # Print results
  cat("\n=== Overall Association Results ===\n")
  cat(sprintf("Number of metabolites used: %d\n", length(selected_metabolites)))
  
  cat("\n1. Mantel Test Results:\n")
  cat(sprintf("Correlation: r = %.3f\n", mantel_result$statistic))
  cat(sprintf("P-value: %.3f\n", mantel_result$signif))
  
  cat("\n2. Procrustes Analysis Results:\n")
  cat(sprintf("Correlation: %.3f\n", sqrt(1 - proc_result$ss)))
  cat(sprintf("P-value: %.3f\n", proc_result$signif))
  
  cat("\n3. Co-inertia Analysis Results:\n")
  cat(sprintf("RV coefficient: %.3f\n", coin_result$RV))
  
  cat("\n=== Individual Analysis by IRIS Groups ===\n")
  print(tapply(individual_results$coinertia_distance, 
               individual_results$IRIS, 
               function(x) c(mean = mean(x), sd = sd(x))))
  
  return(list(
    mantel = mantel_result,
    procrustes = proc_result,
    coinertia = coin_result,
    individual = individual_results
  ))
}

#' Plot association results with IRIS grouping
#' @param results Results from analyze_associations
plot_associations <- function(results) {
  # Co-inertia ordination plot
  p1 <- ggplot(results$individual) +
    geom_segment(aes(x = gut_x, y = gut_y,
                     xend = meta_x, yend = meta_y,
                     color = IRIS),
                 arrow = arrow(length = unit(0.2, "cm"))) +
    geom_point(aes(x = gut_x, y = gut_y), color = "blue", size = 2) +
    geom_point(aes(x = meta_x, y = meta_y), color = "red", size = 2) +
    ggrepel::geom_text_repel(
      aes(x = gut_x, y = gut_y, label = sample),
      size = 3,
      max.overlaps = 10
    ) +scale_color_manual(values = c(iris_color,"grey50") )+
    theme_minimal() +
    labs(title = "Overall Co-inertia Analysis",
         subtitle = sprintf("RV coefficient = %.3f", results$coinertia$RV),
         caption = "Blue: Microbiome positions\nRed: Metabolome positions",
         color = "IRIS Group")
  
  # Individual distance distribution by IRIS
  p2 <- ggplot(results$individual, 
               aes(x = IRIS, y = coinertia_distance, fill = IRIS)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    theme_minimal() +
    labs(title = "Individual Association Strengths by IRIS Group",
         x = "IRIS Group",
         y = "Co-inertia Distance",
         fill = "IRIS Group")+scale_fill_manual(values = c(iris_color,"grey50"))
  
  # Save plots
  ggsave("overall_coinertia.pdf", p1, width = 12, height = 10)
  ggsave("individual_distances_by_IRIS.pdf", p2, width = 10, height = 8)
  
  return(list(
    ordination = p1,
    distances = p2
  ))
}




# Load required packages
library(ggplot2)
library(ggrepel)

#' Plot metabolite selection results in four quadrants of a single plot
#' @param selection_results_list List of selection results for each body site
#' @param metabolomics_class Data frame containing metabolite information
#' @param sites Names of body sites in order: (+,+), (-,+), (-,-), (+,-)
plot_quadrant_metabolite_selection <- function(selection_results_list, metabolomics_class, sites) {
  # 创建代谢物ID和名称的对应关系
  id_to_name <- setNames(metabolomics_class$Compound.name, 
                         rownames(metabolomics_class))
  
  # 创建空的数据框来存储所有结果
  plot_data <- data.frame()
  
  # 处理每个部位的数据
  for(i in seq_along(sites)) {
    site <- sites[i]
    temp_data <- selection_results_list[[site]]$all_results
    temp_data$compound_name <- id_to_name[temp_data$metabolite]
    temp_data$site <- site
    temp_data$selected <- temp_data$n_sig_cors > 0
    
    # 根据象限设置坐标
    if(i == 1) {  # 第一象限
      temp_data$plot_x <- abs(temp_data$max_cor)
      temp_data$plot_y <- temp_data$n_sig_cors
    } else if(i == 2) {  # 第二象限
      temp_data$plot_x <- -abs(temp_data$max_cor)
      temp_data$plot_y <- temp_data$n_sig_cors
    } else if(i == 3) {  # 第三象限
      temp_data$plot_x <- -abs(temp_data$max_cor)
      temp_data$plot_y <- -temp_data$n_sig_cors
    } else {  # 第四象限
      temp_data$plot_x <- abs(temp_data$max_cor)
      temp_data$plot_y <- -temp_data$n_sig_cors
    }
    
    # 标记每个部位的前10个代谢物
    temp_data$is_top10 <- FALSE
    temp_data$is_top10[order(-temp_data$n_sig_cors)[1:10]] <- TRUE
    
    plot_data <- rbind(plot_data, temp_data)
  }
  
  # Create plot
  p <- ggplot(plot_data, 
              aes(x = plot_x, y = plot_y)) +
    # Add grey points for unselected metabolites
    geom_point(data = subset(plot_data, !selected), 
               color = "grey80", 
               alpha = 0.5) +
    # Add colored points for selected metabolites
    geom_point(data = subset(plot_data, selected),
               aes(color = site)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    # Add site labels
    annotate("text", 
             x = c(max(plot_data$plot_x), -max(plot_data$plot_x), 
                   -max(plot_data$plot_x), max(plot_data$plot_x)) * 0.8,
             y = c(max(plot_data$plot_y), max(plot_data$plot_y), 
                   -max(plot_data$plot_y), -max(plot_data$plot_y)) * 0.8,
             label = sites,
             fontface = "bold",
             size = 5) +
    # Add labels for top 10 metabolites
    geom_text_repel(
      data = subset(plot_data, is_top10 & selected),
      aes(label = compound_name),
      size = 2,
      box.padding = 0.5,
      point.padding = 0.1,
      force = 10,
      max.overlaps = Inf,
      min.segment.length = 0,
      segment.color = "grey50",
      direction = "both"
    ) +
    theme_minimal() +
    labs(title = "Metabolite Selection Results Across Body Sites",
         x = "Maximum Correlation",
         y = "Number of Significant Correlations",
         color = "Body Site") +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Save plot
  ggsave("quadrant_metabolite_selection.pdf", p, width = 15, height = 15)
  
  return(p)
}