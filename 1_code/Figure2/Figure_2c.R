setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")

library(tidyverse)
library(tidymass)
load("3_data_analysis/oral_microbiome/data_preparation/object_cross_section")

gut_object <- object_cross_section

load(
  "3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section"
)

metabolomics_object <- object_cross_section
metabolite_annotation <- readxl::read_excel(
  "3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx"
)
####only remain the genus level
library(microbiomedataset)

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

library(plyr)

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


load("3_data_analysis/oral_microbiome/data_preparation/object_cross_section")

oral_object <- object_cross_section

load(
  "3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section"
)

metabolomics_object <- object_cross_section

####only remain the genus level
library(microbiomedataset)

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

library(plyr)

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

load("3_data_analysis/skin_microbiome/data_preparation/object_cross_section")

skin_object <- object_cross_section

load(
  "3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section"
)

metabolomics_object <- object_cross_section

####only remain the genus level
library(microbiomedataset)

skin_object <-
  skin_object %>%
  summarize_variables(what = "sum_intensity", group_by = "Genus")

##only remain the genus at least 10% samples
dim(skin_object)

non_zero_per <-
  apply(skin_object, 1, function(x) {
    sum(x != 0) / ncol(skin_object)
  })

idx <-
  which(non_zero_per > 0.1)

skin_object <-
  skin_object[idx, ]


skin_object <-
  skin_object %>%
  transform2relative_intensity()



##
##adjust BMI, sex, and IRIS, ethnicity
library(tidyverse)
library(ggpubr)
library(rstatix)

skin_expression_data <-
  extract_expression_data(skin_object) %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

library(plyr)

skin_sample_info <-
  skin_object@sample_info

#######adjust BMI, sex, and IRIS, ethnicity
skin_expression_data <-
  lm_adjust(expression_data = skin_expression_data,
            sample_info = skin_sample_info,
            threads = 3)

skin_temp_object <- skin_object
skin_temp_object@expression_data <- skin_expression_data



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


load("3_data_analysis/nasal_microbiome/data_preparation/object_cross_section")

nasal_object <- object_cross_section

load(
  "3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section"
)

metabolomics_object <- object_cross_section

####only remain the genus level
library(microbiomedataset)

nasal_object <-
  nasal_object %>%
  summarize_variables(what = "sum_intensity", group_by = "Genus")

##only remain the genus at least 10% samples
dim(nasal_object)

non_zero_per <-
  apply(nasal_object, 1, function(x) {
    sum(x != 0) / ncol(nasal_object)
  })

idx <-
  which(non_zero_per > 0.1)

nasal_object <-
  nasal_object[idx, ]


nasal_object <-
  nasal_object %>%
  transform2relative_intensity()



##
##adjust BMI, sex, and IRIS, ethnicity
library(tidyverse)
library(ggpubr)
library(rstatix)

nasal_expression_data <-
  extract_expression_data(nasal_object) %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

library(plyr)

nasal_sample_info <-
  nasal_object@sample_info

#######adjust BMI, sex, and IRIS, ethnicity
nasal_expression_data <-
  lm_adjust(expression_data = nasal_expression_data,
            sample_info = nasal_sample_info,
            threads = 3)

nasal_temp_object <- nasal_object
nasal_temp_object@expression_data <- nasal_expression_data



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
select_significant_metabolites <-
  function(gut_microbiome,
           metabolome,
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
    for (i in 1:n_metabolites) {
      # Calculate correlations with all species
      cors <- numeric(n_species)
      p_values <- numeric(n_species)
      
      for (j in 1:n_species) {
        cor_test <- cor.test(as.numeric(meta_matched[i, ]),
                             as.numeric(gut_matched[j, ]),
                             method = "spearman")
        cors[j] <- cor_test$estimate
        p_values[j] <- cor_test$p.value
      }
      
      # Count significant correlations
      sig_strong_cors <- abs(cors) >= cor_threshold &
        p_values < p_threshold
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
    cat(sprintf(
      "Significant metabolites found: %d\n",
      length(sig_metabolites)
    ))
    cat(sprintf("Using correlation threshold: %.2f\n", cor_threshold))
    cat(sprintf("Using p-value threshold: %.3f\n", p_threshold))
    
    # Sort results by number of significant correlations
    results <- results[order(-results$n_sig_cors), ]
    
    # Print top metabolites
    cat("\nTop 10 metabolites by number of significant correlations:\n")
    print(head(results[, c("metabolite", "n_sig_cors", "max_cor")], 10))
    
    return(
      list(
        all_results = results,
        significant_metabolites = sig_metabolites,
        significant_indices = sig_idx
      )
    )
  }

#' Plot metabolite selection results with compound names
#' @param selection_results Results from select_significant_metabolites function
#' @param metabolomics_class Data frame containing metabolite class information






# Load required packages
library(ggplot2)
library(ggrepel)

#' Plot metabolite selection results in four quadrants of a single plot
#' @param selection_results_list List of selection results for each body site
#' @param metabolomics_class Data frame containing metabolite information
#' @param sites Names of body sites in order: (+,+), (-,+), (-,-), (+,-)
plot_quadrant_metabolite_selection <-
  function(selection_results_list,
           metabolomics_class,
           sites) {
    # 创建代谢物ID和名称的对应关系
    id_to_name <- setNames(metabolomics_class$Compound.name,
                           metabolomics_class$variable_id)
    
    # 创建空的数据框来存储所有结果
    plot_data <- data.frame()
    
    # 处理每个部位的数据
    for (i in seq_along(sites)) {
      site <- sites[i]
      temp_data <-
        selection_results_list[[site]]$all_results
      temp_data$compound_name <- id_to_name[temp_data$metabolite]
      temp_data$site <- site
      temp_data$selected <- temp_data$n_sig_cors > 0
      
      # 根据象限设置坐标
      if (i == 1) {
        # 第一象限
        temp_data$plot_x <- abs(temp_data$max_cor)
        temp_data$plot_y <- temp_data$n_sig_cors
      } else if (i == 2) {
        # 第二象限
        temp_data$plot_x <- -abs(temp_data$max_cor)
        temp_data$plot_y <- temp_data$n_sig_cors
      } else if (i == 3) {
        # 第三象限
        temp_data$plot_x <- -abs(temp_data$max_cor)
        temp_data$plot_y <- -temp_data$n_sig_cors
      } else {
        # 第四象限
        temp_data$plot_x <- abs(temp_data$max_cor)
        temp_data$plot_y <- -temp_data$n_sig_cors
      }
      
      # 标记每个部位的前10个代谢物
      temp_data$is_top10 <- FALSE
      temp_data$is_top10[order(-temp_data$n_sig_cors)[1:10]] <- TRUE
      
      plot_data <- rbind(plot_data, temp_data)
    }
    
    # Create plot
    p <- ggplot(plot_data, aes(x = plot_x, y = plot_y)) +
      # Add grey points for unselected metabolites
      geom_point(
        data = subset(plot_data, !selected),
        color = "grey80",
        alpha = 0.5
      ) +
      # Add colored points for selected metabolites
      geom_point(data = subset(plot_data, selected), 
                 aes(color = site),
                 alpha = 0.6,
                 size = 5) +
      scale_color_manual(values = body_site_color) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      # Add site labels
      annotate(
        "text",
        x = c(
          max(plot_data$plot_x),-max(plot_data$plot_x),
          -max(plot_data$plot_x),
          max(plot_data$plot_x)
        ) * 0.8,
        y = c(
          max(plot_data$plot_y),
          max(plot_data$plot_y),
          -max(plot_data$plot_y),-max(plot_data$plot_y)
        ) * 0.8,
        label = sites,
        fontface = "bold",
        size = 5
      ) +
      # Add labels for top 10 metabolites
      geom_text_repel(
        data = subset(plot_data, is_top10 & selected),
        aes(plot_x, plot_y, label = compound_name),
        size = 2,
        box.padding = 0.5,
        point.padding = 0.1,
        force = 10,
        max.overlaps = Inf,
        min.segment.length = 0,
        segment.color = "grey50",
        direction = "both"
      ) +
      theme_bw() +
      labs(# title = "Metabolite Selection Results Across Body Sites",
        x = "Maximum Correlation", y = "Number of Significant Correlations", color = "Site") +
      theme(plot.title = element_text(hjust = 0.5))
    
    # # Save plot
    # ggsave("quadrant_metabolite_selection.pdf",
    #        p,
    #        width = 15,
    #        height = 15)
    
    return(p)
  }

gut_result <- select_significant_metabolites(gut_temp_object@expression_data,
                                             metabolome = metabolomics_temp_object@expression_data)
oral_result <- select_significant_metabolites(oral_temp_object@expression_data,
                                              metabolome = metabolomics_temp_object@expression_data)
skin_result <- select_significant_metabolites(skin_temp_object@expression_data,
                                              metabolome = metabolomics_temp_object@expression_data)
nasal_result <- select_significant_metabolites(nasal_temp_object@expression_data,
                                               metabolome = metabolomics_temp_object@expression_data)

# 假设我们有四个部位的选择结果
selection_results_list <- list(
  "gut" = gut_result,
  "oral" = oral_result,
  "skin" = skin_result,
  "nasal" = nasal_result
)

# 定义部位名称
sites <- c("gut", "oral", "skin", "nasal")
rownames(metabolite_annotation) <- metabolite_annotation$variable_id
# 绘制多部位比较图
p <- plot_quadrant_metabolite_selection(
  selection_results_list = selection_results_list,
  metabolomics_class = metabolite_annotation,
  sites = sites
)
p

ggsave(
  p,
  filename = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_2/figure_2c.pdf"
  ),
  width = 8,
  height = 6
)
