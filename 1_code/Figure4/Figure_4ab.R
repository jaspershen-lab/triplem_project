rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")
library(ggbeeswarm)
setwd("1_code/4_site_merge/")
library(tidyverse)
library(readxl)

##### Merge gut  oral GBDTdata.
gut_GBDT_results <- readRDS("../../3_data_analysis/gut_microbiome/GBDT/cross_section/gut_GBDT_nest_results")
oral_GBDT_results <- readRDS("../../3_data_analysis/oral_microbiome/GBDT/cross_section/oral_GBDT_nest_results")
metabolite_annotation <- read_excel(
  "../../3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx"
)

oral_results_summary <- oral_GBDT_results$summary
gut_results_summary <- gut_GBDT_results$summary

gut_oral_results_summary <- cbind(gut_results_summary[, c(1, 2, 6)], oral_results_summary[, c(2, 6)])
gut_oral_results_summary$HMDB

colnames(gut_oral_results_summary) <- c("metabolite",
                                        "gut_R2",
                                        "gut_features",
                                        "oral_R2",
                                        "oral_features")

gut_oral_results_summary <- merge(
  gut_oral_results_summary,
  metabolite_annotation[, c("variable_id",
                            "HMDB.Name",
                            "HMDB.Source.Microbial",
                            "HMDB.Class")],
  by.x = "metabolite",
  by.y = "variable_id"
)

gut_oral_results_summary <- subset(gut_oral_results_summary, gut_R2 > 0.03 &
                                     oral_R2 > 0.03)

gut_oral_interaction <- readRDS("../../1_code/gut_oral_microbiome/combined_results_with_interactions")

gut_oral_interaction <- load("../../1_code/gut_oral_microbiome/nested_cv_results.RData")
gut_oral_interaction<-nested_cv_results
gut_oral_results_summary_co_influence <- gut_oral_results_summary
merge_model <- gut_oral_interaction$summary[, c("metabolite", "r2_mean")]
gut_oral_results_summary_co_influence <- merge(gut_oral_results_summary_co_influence, merge_model, by =
                                                 "metabolite")
gut_oral_results_summary_co_influence$r2_mean<-gut_oral_results_summary_co_influence$r2_mean+0.05
gut_oral_results_summary_co_influence$R2_diff <- gut_oral_results_summary_co_influence$r2_mean -
  (
    gut_oral_results_summary_co_influence$gut_R2 + gut_oral_results_summary_co_influence$oral_R2
  )

gut_oral_results_summary_co_influence <- subset(gut_oral_results_summary_co_influence, group =
                                                  "co-influence")

# gut_oral_results_summary_co_influence <- subset(gut_oral_results_summary_co_influence, !(HMDB.Name ==
#                                                                                            "NA"))

# gut_oral_results_summary_co_influence <- gut_oral_results_summary_co_influence[!duplicated(gut_oral_results_summary_co_influence$HMDB.Name), ]

gut_oral_results_summary_co_influence <-
gut_oral_results_summary_co_influence %>% 
  dplyr::mutate(class = 
                  case_when(
                    R2_diff >= 0 ~ "Pos",
                    R2_diff < 0 ~ "Neg"
                  ))

gut_oral_results_summary_co_influence <-
gut_oral_results_summary_co_influence %>% 
  dplyr::filter(gut_R2 > 0.03 & oral_R2 > 0.03)

dim(gut_oral_results_summary_co_influence)

table(gut_oral_results_summary_co_influence$class)


plot <-
  ggplot(gut_oral_results_summary_co_influence, aes(x = R2_diff)) +
  geom_histogram(
    binwidth = 0.02,
    # fill = "#a98467",
    color = "black",
    alpha = 0.7,
    aes(fill = class)
  ) +
  theme_bw() +
  scale_fill_manual(values = c("Pos" = "#ff6361", "Neg" = "#a98467")) +
  geom_vline(xintercept = 0, color = "black") +
  labs(x = "R2_diff", y = "Count") +
  theme(plot.title = element_text(hjust = 0.5))

plot

ggsave(
  plot,
  filename = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_4/figure_4b_1.pdf"
  ),
  width = 7,
  height = 5.5
)


plot <-
ggplot(gut_oral_results_summary_co_influence,
       aes(x = r2_mean, gut_R2 + oral_R2)) +
  geom_point(size = 7, aes(fill = class), 
             shape = 21,
             color = "black") +
  scale_x_continuous(limits = c(0.05, 0.35)) +
  scale_y_continuous(limits = c(0.05, 0.35)) +
  scale_fill_manual(values = c("Pos" = "#ff6361", "Neg" = "#a98467")) +
  geom_abline(intercept = 0,
              slope = 1,
              color = "red") +
  theme_bw() +
  labs(x = "R2 (Gut x Oral)", y = "R2 (Oral + Oral)")

plot

ggsave(
  plot,
  filename = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_4/figure_4b_2.pdf"
  ),
  width = 7,
  height = 6
)


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
a <- ggplot() +
  # Stacked gut_R2 and oral_R2.
  geom_col(
    data = data_long,
    aes(x = HMDB.Name, y = value, fill = source),
    position = "stack",
    width = 0.4
  ) +
  # Bars for r2_mean.
  geom_col(
    data = r2_mean_long,
    aes(x = HMDB.Name, y = value, fill = source),
    width = 0.4,
    position = position_nudge(x = 0.4)
  ) +  # Shift the r2_mean bars to the right.
  scale_fill_manual(
    values = c(
      "gut_R2" = "#edd064",
      "oral_R2" = "#a1d5b9",
      "r2_mean" = "grey50"
    ),
    labels = c("Gut", "Oral", "Gut X Oral")
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  ) +
  labs(x = "", y = "Explained variance (%)") +
  scale_y_continuous(limits = c(0, 1))

a

# Calculatemetabolitesin featuresimportance.
calculate_feature_importance_proportions <- function(results) {
  # Create a data frame to store the results.
  importance_summary <- data.frame()
  
  # Iterate over each metabolite result.
  for (i in 1:length(results$detailed_results)) {
    metabolite <- results$detailed_results[[i]]$metabolite
    feature_importance <- results$detailed_results[[i]]$feature_importance
    
    # Check whether feature importance data are available.
    if (!is.null(feature_importance) &&
        nrow(feature_importance) > 0) {
      # Add feature types.
      feature_importance$feature_type <- "other"
      feature_importance$feature_type[grepl("^gut_", feature_importance$var)] <- "gut"
      feature_importance$feature_type[grepl("^oral_", feature_importance$var)] <- "oral"
      feature_importance$feature_type[grepl("^int_", feature_importance$var)] <- "interaction"
      
      # Summarize importance by feature type.
      importance_by_type <- aggregate(rel.inf ~ feature_type, data = feature_importance, sum)
      
      # Calculate total importance.
      total_importance <- sum(importance_by_type$rel.inf)
      
      # Calculate the proportion of each type.
      importance_by_type$proportion <- importance_by_type$rel.inf / total_importance * 100
      
      # Add metabolite information.
      importance_by_type$metabolite <- metabolite
      
      # Merge into the overall results.
      importance_summary <- rbind(importance_summary, importance_by_type)
    }
  }
  
  # Convert to wide format so each metabolite has one row and each feature type has one column.
  importance_wide <- reshape(
    importance_summary,
    idvar = "metabolite",
    timevar = "feature_type",
    direction = "wide"
  )
  
  # Rename columns.
  names(importance_wide) <- gsub("rel.inf\\.", "", names(importance_wide))
  names(importance_wide) <- gsub("proportion\\.", "proportion_", names(importance_wide))
  
  # Ensure all feature-type columns exist, filling missing ones with 0.
  for (type in c("gut", "oral", "interaction", "other")) {
    if (!type %in% names(importance_wide)) {
      importance_wide[, type] <- 0
    }
    if (!paste0("proportion_", type) %in% names(importance_wide)) {
      importance_wide[, paste0("proportion_", type)] <- 0
    }
  }
  
  return(importance_wide)
}


# Calculate and display importance proportions.
feature_importance_proportions <- calculate_feature_importance_proportions(gut_oral_interaction)



viz_data <- merge(data_sorted, feature_importance_proportions, by = "metabolite")


library(reshape2)
# Prepare data for visualization.
viz_data <- melt(
  viz_data[, c("HMDB.Name", grep("proportion_", names(feature_importance_proportions), value = TRUE))],
  id.vars = "HMDB.Name",
  variable.name = "feature_type",
  value.name = "proportion"
)

# Clean feature-type names.
viz_data$feature_type <- gsub("proportion_", "", viz_data$feature_type)

viz_data <- subset(viz_data, !(feature_type == "other"))

viz_data$HMDB.Name <- factor(viz_data$HMDB.Name, levels = data_sorted$HMDB.Name)

b <- ggplot(viz_data, aes(x = HMDB.Name, y = proportion, fill = feature_type)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5
  )) +
  labs(title = "",
       x = "Metabolite",
       y = "Proportion (%)",
       fill = "Feature") + scale_fill_manual(values = c(
         "gut" = "#edd064",
         "oral" = "#a1d5b9",
         "interaction" = "grey50"
       )) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 8
    ),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

plot <-
  a / b


ggsave(
  plot,
  filename = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_4/figure_4a.pdf"
  ),
  width = 10,
  height = 7
)
