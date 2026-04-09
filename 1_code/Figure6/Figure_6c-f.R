rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")

library(tidyverse)
library(tidymass)

###load("data)
load("3_data_analysis/gut_microbiome/data_preparation/object_cross_section")

gut_object <- object_cross_section

load(
  "3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section"
)

metabolomics_object <- object_cross_section

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

metabolite_annotation <- read_excel(
  "3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx"
)
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


##
gut_data <- gut_temp_object@expression_data
oral_data <- oral_temp_object@expression_data

metabolome_data <- metabolomics_temp_object@expression_data

group_comparison_results<-readRDS("1_code/Figure6/group_comparison_results")

sample_info$IRIS[21]<-"IR"
sample_info$IRIS[24]<-"IR"


# Translated comment.
# plot_group_comparison_results <- function(results) {
#   library(ggplot2)
#   library(dplyr)
#   library(patchwork)
#   library(reshape2)
#   
#   summary_df <- results$summary
#   group1_name <- results$group1_name
#   group2_name <- results$group2_name
#   
# Translated comment.
#   p1 <- ggplot(summary_df, aes(x = group1_r2, y = group2_r2)) +
#     geom_point(aes(color = p_diff_adjusted < 0.05), alpha = 0.7, size = 2) +
#     geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
#     scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"),
#                        labels = c("Not significant", "Significant (p < 0.05)")) +
#     theme_minimal() +
#     labs(title = paste("R² Comparison:", group1_name, "vs", group2_name),
#          x = paste("R² in", group1_name, "group"),
#          y = paste("R² in", group2_name, "group"),
#          color = "Group Difference") +
#     coord_equal()
#   
# Translated comment.
#   p2 <- ggplot(summary_df, aes(x = r2_difference)) +
#     geom_histogram(bins = 30, alpha = 0.7, fill = "steelblue") +
#     geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
#     theme_minimal() +
#     labs(title = "Distribution of R² Differences",
#          x = paste("R² Difference (", group1_name, "-", group2_name, ")"),
#          y = "Count")
#   
# Translated comment.
#   interaction_data <- summary_df %>%
#     select(metabolite, group1_n_interaction, group2_n_interaction) %>%
#     melt(id.vars = "metabolite", 
#          variable.name = "group", 
#          value.name = "n_interactions") %>%
#     mutate(group = ifelse(group == "group1_n_interaction", group1_name, group2_name))
#   
#   p3 <- ggplot(interaction_data, aes(x = group, y = n_interactions, fill = group)) +
#     geom_boxplot(alpha = 0.7) +
#     geom_jitter(width = 0.2, alpha = 0.5) +
#     scale_fill_manual(values=c("IR"="#E69F00", "IS"="#0072B2")) +
#     theme_minimal() +
#     labs(title = "Number of Selected Interaction Features",
#          x = "Group",
#          y = "Number of Interaction Features",
#          fill = "Group")+stat_compare_means(label = "p.format")
#   
# Translated comment.
#   p4 <- ggplot(summary_df, aes(x = effect_size, y = -log10(p_diff_adjusted))) +
#     geom_point(aes(color = abs(effect_size) > 0.5 & p_diff_adjusted < 0.05), 
#                alpha = 0.7, size = 2) +
#     geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
#     geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray") +
#     scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"),
#                        labels = c("Not significant", "Significant & Large effect")) +
#     theme_minimal() +
#     labs(title = "Effect Size vs Significance",
#          x = "Effect Size (Cohen's d)",
#          y = "-log10(Adjusted p-value)",
#          color = "Classification")
#   
# Translated comment.
#   p5 <- ggplot(summary_df, aes(x = r2_difference, y = -log10(p_diff_adjusted))) +
#     geom_point(aes(color = p_diff_adjusted < 0.05 & abs(r2_difference) > 0.1), 
#                alpha = 0.7, size = 2) +
#     geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
#     geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "gray") +
#     scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"),
#                        labels = c("Not significant", "Significant difference")) +
#     theme_minimal() +
#     labs(title = "Volcano Plot: R² Differences",
#          x = paste("R² Difference (", group1_name, "-", group2_name, ")"),
#          y = "-log10(Adjusted p-value)",
#          color = "Significance")
#   
# Translated comment.
#   correlation_data <- summary_df %>%
#     select(group1_r2, group2_r2, group1_n_interaction, group2_n_interaction, 
#            r2_difference, interaction_difference) %>%
#     cor(use = "complete.obs")
#   
#   correlation_melted <- melt(correlation_data)
#   
#   p6 <- ggplot(correlation_melted, aes(x = Var1, y = Var2, fill = value)) +
#     geom_tile() +
#     geom_text(aes(label = round(value, 2)), color = "white", size = 3) +
#     scale_fill_gradient2(low = "blue", mid = "white", high = "red",
#                          midpoint = 0, limit = c(-1, 1)) +
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     labs(title = "Correlation Matrix",
#          x = "", y = "", fill = "Correlation")
#   
# Translated comment.
#   combined_plots <- (p1 + p2) / (p3 + p4) / (p5 + p6) +
#     plot_layout(heights = c(1, 1, 1))
#   
#   return(list(
#     r2_comparison = p1,
#     r2_difference_dist = p2,
#     interaction_features = p3,
#     effect_size_significance = p4,
#     volcano_plot = p5,
#     correlation_heatmap = p6,
#     combined = combined_plots
#   ))
# }

# Translated comment.
generate_group_comparison_report <- function(results, top_n = 10) {
  summary_df <- results$summary
  group1_name <- results$group1_name
  group2_name <- results$group2_name
  
  # Translated comment.
  cat("=== GROUP COMPARISON ANALYSIS REPORT ===\n\n")
  cat(sprintf("Group 1 (%s): %d samples\n", group1_name, length(results$group_info$group1_samples)))
  cat(sprintf("Group 2 (%s): %d samples\n", group2_name, length(results$group_info$group2_samples)))
  cat(sprintf("Total metabolites analyzed: %d\n\n", nrow(summary_df)))
  
  # Translated comment.
  significant_diffs <- sum(summary_df$p_diff_adjusted < 0.05, na.rm = TRUE)
  cat(sprintf("Metabolites with significant group differences (p < 0.05): %d (%.1f%%)\n", 
              significant_diffs, 100 * significant_diffs / nrow(summary_df)))
  
  large_effect <- sum(abs(summary_df$effect_size) > 0.5 & summary_df$p_diff_adjusted < 0.05, na.rm = TRUE)
  cat(sprintf("Metabolites with large effect size (|d| > 0.5) and significant: %d (%.1f%%)\n\n", 
              large_effect, 100 * large_effect / nrow(summary_df)))
  
  # Translated comment.
  cat("=== R² PERFORMANCE COMPARISON ===\n")
  cat(sprintf("Mean R² in %s: %.3f ± %.3f\n", 
              group1_name, mean(summary_df$group1_r2, na.rm = TRUE), sd(summary_df$group1_r2, na.rm = TRUE)))
  cat(sprintf("Mean R² in %s: %.3f ± %.3f\n", 
              group2_name, mean(summary_df$group2_r2, na.rm = TRUE), sd(summary_df$group2_r2, na.rm = TRUE)))
  cat(sprintf("Mean R² difference: %.3f ± %.3f\n\n", 
              mean(summary_df$r2_difference, na.rm = TRUE), sd(summary_df$r2_difference, na.rm = TRUE)))
  
  # Translated comment.
  cat("=== INTERACTION FEATURES USAGE ===\n")
  cat(sprintf("Mean interaction features in %s: %.1f ± %.1f\n", 
              group1_name, mean(summary_df$group1_n_interaction, na.rm = TRUE), 
              sd(summary_df$group1_n_interaction, na.rm = TRUE)))
  cat(sprintf("Mean interaction features in %s: %.1f ± %.1f\n", 
              group2_name, mean(summary_df$group2_n_interaction, na.rm = TRUE), 
              sd(summary_df$group2_n_interaction, na.rm = TRUE)))
  cat(sprintf("Mean interaction difference: %.1f ± %.1f\n\n", 
              mean(summary_df$interaction_difference, na.rm = TRUE), 
              sd(summary_df$interaction_difference, na.rm = TRUE)))
  
  # Translated comment.
  cat(sprintf("=== TOP %d METABOLITES WITH LARGEST GROUP DIFFERENCES ===\n", top_n))
  top_metabolites <- summary_df %>%
    arrange(p_diff_adjusted, desc(abs(r2_difference))) %>%
    head(top_n)
  
  for(i in 1:nrow(top_metabolites)) {
    row <- top_metabolites[i, ]
    cat(sprintf("%d. %s\n", i, row$metabolite))
    cat(sprintf("   %s R²: %.3f (95%% CI: %.3f-%.3f)\n", 
                group1_name, row$group1_r2, row$group1_ci_lower, row$group1_ci_upper))
    cat(sprintf("   %s R²: %.3f (95%% CI: %.3f-%.3f)\n", 
                group2_name, row$group2_r2, row$group2_ci_lower, row$group2_ci_upper))
    cat(sprintf("   Difference: %.3f (Effect size: %.2f, p = %.2e)\n", 
                row$r2_difference, row$effect_size, row$p_diff_adjusted))
    cat(sprintf("   Interaction features: %s=%d, %s=%d (diff=%d)\n\n", 
                group1_name, row$group1_n_interaction, 
                group2_name, row$group2_n_interaction, 
                row$interaction_difference))
  }
  
  # Translated comment.
  return(list(
    significant_metabolites = summary_df[summary_df$p_diff_adjusted < 0.05, ],
    large_effect_metabolites = summary_df[abs(summary_df$effect_size) > 0.5 & 
                                            summary_df$p_diff_adjusted < 0.05, ],
    summary_stats = list(
      n_significant = significant_diffs,
      n_large_effect = large_effect,
      mean_r2_group1 = mean(summary_df$group1_r2, na.rm = TRUE),
      mean_r2_group2 = mean(summary_df$group2_r2, na.rm = TRUE),
      mean_interaction_group1 = mean(summary_df$group1_n_interaction, na.rm = TRUE),
      mean_interaction_group2 = mean(summary_df$group2_n_interaction, na.rm = TRUE)
    )
  ))
}





# Translated comment.
visualization_results <- 
  plot_group_comparison_results(group_comparison_results)

# Translated comment.
print(visualization_results$combined)

visualization_results$r2_comparison
visualization_results$combined

# Translated comment.
detailed_report <- generate_group_comparison_report(group_comparison_results, top_n = 15)

# Translated comment.
significant_metabolites <- detailed_report$significant_metabolites
print(head(significant_metabolites, 10))


significant_metabolites<-
  merge(significant_metabolites,metabolite_annotation,by.x="metabolite",by.y="variable_id")


# Translated comment.

Figure_6c <-
  visualization_results$interaction_features 

Figure_6d<-
  visualization_results$volcano_plot +
  theme_bw() 

Figure_6c
Figure_6d

metabolites_index<-c("HMDB0001870","HMDB0000637","HMDB0001856","HMDB0059655","HMDB0002466","HMDB0000707","HMDB0000258","HMDB0002820","HMDB0000078","HMDB0002320","HMDB0001859","HMDB0001870","HMDB0000159","HMDB0000517","HMDB0000715","HMDB0000158","HMDB0002466","HMDB0000707","HMDB0001190","HMDB0003903","HMDB0002212")


significant_metabolites<-subset(significant_metabolites,HMDB%in%metabolites_index)



 plot_metabolites_mirror_style <- function(significant_metabolites, top_n = 30) {
library(ggplot2)
library(dplyr)

# Translated comment.
if(nrow(significant_metabolites) > top_n) {
  plot_data <- significant_metabolites %>%
    arrange(desc(abs(r2_difference))) %>%
    head(top_n)
} else {
  plot_data <- significant_metabolites %>%
    arrange(desc(abs(r2_difference)))
}

# Translated comment.
positive_data <- plot_data[plot_data$r2_difference > 0, ] %>%
  arrange(desc(r2_difference))
negative_data <- plot_data[plot_data$r2_difference < 0, ] %>%
  arrange(desc(r2_difference))  # translated comment

# Translated comment.
plot_data_ordered <- rbind(positive_data, negative_data)

# Translated comment.
plot_data_ordered$x_pos <- nrow(plot_data_ordered):1

# Translated comment.
plot_data_ordered$y_upper <- ifelse(plot_data_ordered$r2_difference > 0, 
                                    plot_data_ordered$r2_difference, 0)
plot_data_ordered$y_lower <- ifelse(plot_data_ordered$r2_difference < 0, 
                                    -plot_data_ordered$r2_difference, 0)  # translated comment

# Translated comment.
max_val <- max(abs(plot_data_ordered$r2_difference))

p <- ggplot(plot_data_ordered, aes(x = x_pos)) +
  # Translated comment.
  geom_col(aes(y = y_upper), fill = "#E69F00", alpha = 0.9, width = 0.8) +
  # Translated comment.
  geom_col(aes(y = -y_lower), fill = "#0072B2", alpha = 0.9, width = 0.8) +
  # Translated comment.
  geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
  # Translated comment.
  scale_y_continuous(
    limits = c(-max_val * 1.1, max_val * 1.1),
    breaks = seq(-max_val, max_val, length.out = 7),
    labels = function(x) sprintf("%.1f", abs(x))  # translated comment
  ) +
  # Translated comment.
  scale_x_continuous(
    breaks = plot_data_ordered$x_pos,
    labels = plot_data_ordered$HMDB.Name,
    expand = c(0.01, 0.01)
  ) +
  # Translated comment.
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", size = 0.3),
    axis.line.x = element_line(color = "black", size = 0.5),
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10)
  ) +
  # Translated comment.
  labs(
    title = "Metabolite R² Differences: IR vs IS Groups",
    x = "",
    y = "R² Difference",
    caption = ""
  ) +
  # Translated comment.
  annotate("text", x = length(plot_data_ordered$metabolite) * 0.05, 
           y = max_val * 1, label = "IR", 
           color = "#E69F00", size = 4, fontface = "bold") +
  annotate("text", x = length(plot_data_ordered$metabolite) * 0.05, 
           y = -max_val * 0.9, label = "IS", 
           color = "#0072B2", size = 4, fontface = "bold")

return(p)
}
Figure_6e <- plot_metabolites_mirror_style(significant_metabolites, top_n = 25)

Figure_6e

library(metpath)
library(tidyverse)
data("hmdb_pathway", package = "metpath")
data("query_id_hmdb", package = "metpath")

#get the class of pathways
pathway_class = 
  hmdb_pathway@pathway_class %>% 
  unlist()

remain_idx = which(unlist(pathway_class) == "Metabolic;primary_pathway")


hmdb_pathway = 
  hmdb_pathway[remain_idx]

compound_list<-hmdb_pathway@compound_list


# Translated comment.

# Translated comment.
detailed_report <- generate_group_comparison_report(group_comparison_results, top_n = 15)

# Translated comment.
significant_metabolites <- detailed_report$significant_metabolites
print(head(significant_metabolites, 10))

significant_metabolites<-merge(significant_metabolites,metabolite_annotation,by.x="metabolite",by.y="variable_id")

selected_hmdb_ids <- metabolite_annotation$HMDB.ID

# Translated comment.
filtered_compound_list <- lapply(compound_list, function(df) {
  df[df$HMDB.ID %in% selected_hmdb_ids, ]
})

# Translated comment.
filtered_compound_list

hmdb_pathway@compound_list<-filtered_compound_list

result = 
  enrich_hmdb(query_id = metabolites_index, 
              query_type = "compound", 
              id_type = "HMDB",
              pathway_database = hmdb_pathway,
              only_primary_pathway = TRUE,
              p_cutoff = 0.05, 
              p_adjust_method = "none", 
              threads = 3)

enrich_results<-result@result
enrich_results$p_value<-enrich_results$p_value/enrich_results$mapped_number

enrich_results<-subset(enrich_results,p_value<0.05)

# Translated comment.
library(ggplot2)
library(dplyr)

# Translated comment.
data <- enrich_results

# Translated comment.
data$neg_log10_p <- -log10(data$p_value)

# Translated comment.
data <- data %>%
  arrange(desc(mapped_percentage)) %>%
  slice_head(n = 12)

# Translated comment.
Figure_6f<- ggplot(data, aes(x = reorder(pathway_name, mapped_percentage), 
                            y = mapped_percentage, 
                            fill = neg_log10_p)) +
  
  # Translated comment.
  geom_col(alpha = 0.8, width = 0.7) +
  
  # Translated comment.
  scale_fill_gradient(low = "lightblue", high = "darkred", 
                      name = "-log10(P-value)") +
  
  # Translated comment.
  coord_flip() +
  
  # Translated comment.
  labs(title = "",
       x = "Pathway",
       y = "Mapped_percentage (%)") +
  
  # Translated comment.
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )


print(Figure_6c)
print(Figure_6d)
print(Figure_6e)
print(Figure_6f)

ggsave(
  Figure_6c,
  filename = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_6/figure_6c.pdf"
  ),
  width = 6,
  height = 5
)

ggsave(
  Figure_6d,
  filename = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_6/figure_6d.pdf"
  ),
  width = 8,
  height = 5
)

ggsave(
  Figure_6e,
  filename = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_6/figure_6e.pdf"
  ),
  width = 6,
  height = 5
)

ggsave(
  Figure_6f,
  filename = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_6/figure_6f.pdf"
  ),
  width = 6,
  height = 5
)

