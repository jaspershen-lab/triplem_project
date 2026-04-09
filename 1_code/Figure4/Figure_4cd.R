rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")
library(patchwork)
library(doParallel)
library(tidyverse)
library(tidymass)
library(parallel)
library(progress)
library(readxl)

###load("data)
load("3_data_analysis/gut_microbiome/data_preparation/object_cross_section")

gut_object <- object_cross_section

load(
  "3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section"
)

metabolomics_object <- object_cross_section
metabolite_annotation <- read_excel(
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

dir.create("3_data_analysis/oral_microbiome/spearman/cross_section/",
           recursive = TRUE)

setwd("3_data_analysis/oral_microbiome/spearman/cross_section/")

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



# Translated comment.
select_significant_interaction_metabolites <- function(combined_results,
                                                       min_r2 = 0.1,
                                                       min_interaction_importance = 0.2,
                                                       top_n = 20) {
  # Translated comment.
  metabolite_interaction_importance <- list()
  
  # Translated comment.
  for (i in 1:length(combined_results$detailed_results)) {
    result <- combined_results$detailed_results[[i]]
    metabolite_name <- result$metabolite
    
    # Translated comment.
    if (!is.null(result$feature_importance)) {
      # Translated comment.
      importance_df <- result$feature_importance
      
      # Translated comment.
      interaction_features <- importance_df[grepl("^int_", importance_df$var), ]
      
      if (nrow(interaction_features) > 0) {
        # Translated comment.
        total_importance <- sum(importance_df$rel.inf)
        interaction_importance_sum <- sum(interaction_features$rel.inf)
        interaction_importance_ratio <- interaction_importance_sum / total_importance
        
        # Translated comment.
        metabolite_interaction_importance[[metabolite_name]] <- list(
          metabolite = metabolite_name,
          r2 = result$mean_r2,
          interaction_importance_sum = interaction_importance_sum,
          interaction_importance_ratio = interaction_importance_ratio,
          n_interaction_features = nrow(interaction_features),
          top_interactions = head(interaction_features[order(interaction_features$rel.inf, decreasing = TRUE), ], 10)
        )
      }
    }
  }
  
  # Translated comment.
  importance_df <- do.call(rbind,
                           lapply(metabolite_interaction_importance, function(x) {
                             data.frame(
                               metabolite = x$metabolite,
                               r2 = x$r2,
                               interaction_importance_sum = x$interaction_importance_sum,
                               interaction_importance_ratio = x$interaction_importance_ratio,
                               n_interaction_features = x$n_interaction_features
                             )
                           }))
  
  # Translated comment.
  filtered_metabolites <- importance_df[importance_df$r2 >= min_r2 &
                                          importance_df$interaction_importance_ratio >= min_interaction_importance, ]
  
  # Translated comment.
  ranked_metabolites <- filtered_metabolites[order(filtered_metabolites$interaction_importance_ratio,
                                                   decreasing = TRUE), ]
  
  # Translated comment.
  top_metabolites <- head(ranked_metabolites, top_n)
  
  # Translated comment.
  top_metabolites_details <- lapply(as.character(top_metabolites$metabolite), function(metabolite) {
    metabolite_interaction_importance[[metabolite]]
  })
  names(top_metabolites_details) <- top_metabolites$metabolite
  
  # Translated comment.
  return(
    list(
      summary = top_metabolites,
      details = top_metabolites_details,
      all_metabolites = importance_df
    )
  )
}

# Translated comment.
plot_interaction_effects <- function(interaction_results) {
  # Translated comment.
  summary_df <- interaction_results$summary
  n_to_plot <- min(nrow(summary_df), 20)
  top_metabolites <- summary_df[1:n_to_plot, ]
  
  # Translated comment.
  top_metabolites$metabolite <- factor(top_metabolites$metabolite, levels = top_metabolites$metabolite[order(top_metabolites$interaction_importance_ratio)])
  
  # Translated comment.
  p1 <- ggplot(top_metabolites,
               aes(x = metabolite, y = interaction_importance_ratio)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    theme_minimal() +
    labs(title = "代谢物交互特征重要性占比", x = "代谢物", y = "交互特征重要性比例")
  
  # Translated comment.
  p2 <- ggplot(interaction_results$all_metabolites,
               aes(x = interaction_importance_ratio, y = r2)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "loess", se = TRUE) +
    theme_minimal() +
    labs(title = "交互特征重要性与模型性能(R²)的关系", x = "交互特征重要性比例", y = "R²")
  
  # Translated comment.
  top_interactions <- do.call(rbind, lapply(1:n_to_plot, function(i) {
    metabolite <- as.character(top_metabolites$metabolite[i])
    details <- interaction_results$details[[metabolite]]
    top_5 <- details$top_interactions
    
    if (nrow(top_5) > 0) {
      data.frame(
        metabolite = metabolite,
        feature = top_5$var,
        importance = top_5$rel.inf
      )
    } else {
      data.frame(
        metabolite = character(),
        feature = character(),
        importance = numeric()
      )
    }
  }))
  
  # Translated comment.
  if (nrow(top_interactions) > 0) {
    # Translated comment.
    pattern <- "int_(ASV[0-9]+)_(OTU_[0-9]+)"
    top_interactions$gut_feature <- sub(pattern, "\\1", top_interactions$feature)
    top_interactions$oral_feature <- sub(pattern, "\\2", top_interactions$feature)
    
    # Translated comment.
    heatmap_data <- top_interactions
    # Translated comment.
    heatmap_data$metabolite <- factor(heatmap_data$metabolite, levels = rev(levels(top_metabolites$metabolite)))
    
    # Translated comment.
    p3 <- ggplot(heatmap_data,
                 aes(x = feature, y = metabolite, fill = importance)) +
      geom_tile() +
      scale_fill_viridis_c() +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        title = "代谢物的Top交互特征重要性",
        x = "交互特征",
        y = "代谢物",
        fill = "重要性"
      )
  } else {
    # Translated comment.
    p3 <- ggplot() +
      geom_text(aes(x = 0, y = 0, label = "没有交互特征")) +
      theme_void() +
      labs(title = "代谢物的Top交互特征重要性")
  }
  
  # Translated comment.
  library(patchwork)
  combined_plots <- (p1 | p2) / p3 +
    plot_layout(heights = c(1, 2))
  
  return(
    list(
      importance_bar = p1,
      r2_scatter = p2,
      feature_heatmap = p3,
      combined = combined_plots
    )
  )
}

# Translated comment.
analyze_interaction_features <- function(interaction_results,
                                         selected_metabolites = NULL) {
  if (is.null(selected_metabolites)) {
    # Translated comment.
    selected_metabolites <- interaction_results$summary$metabolite
  }
  
  # Translated comment.
  all_interactions <- do.call(rbind, lapply(as.character(selected_metabolites), function(metabolite) {
    details <- interaction_results$details[[metabolite]]
    if (!is.null(details) &&
        !is.null(details$top_interactions) &&
        nrow(details$top_interactions) > 0) {
      interactions <- details$top_interactions
      data.frame(
        metabolite = metabolite,
        feature = interactions$var,
        importance = interactions$rel.inf
      )
    } else {
      NULL
    }
  }))
  
  if (nrow(all_interactions) == 0) {
    return(NULL)
  }
  
  # Translated comment.
  # Translated comment.
  pattern <- "int_(ASV[0-9]+)_(OTU_[0-9]+)"
  all_interactions$gut_feature <- sub(pattern, "\\1", all_interactions$feature)
  all_interactions$oral_feature <- sub(pattern, "\\2", all_interactions$feature)
  
  # Translated comment.
  gut_counts <- table(all_interactions$gut_feature)
  oral_counts <- table(all_interactions$oral_feature)
  
  # Translated comment.
  feature_pairs <- paste(all_interactions$gut_feature,
                         all_interactions$oral_feature,
                         sep = "_")
  pair_counts <- table(feature_pairs)
  
  # Translated comment.
  weighted_gut_counts <- tapply(all_interactions$importance,
                                all_interactions$gut_feature,
                                sum)
  weighted_oral_counts <- tapply(all_interactions$importance,
                                 all_interactions$oral_feature,
                                 sum)
  
  # Translated comment.
  return(
    list(
      all_interactions = all_interactions,
      gut_frequency = sort(gut_counts, decreasing = TRUE),
      oral_frequency = sort(oral_counts, decreasing = TRUE),
      pair_frequency = sort(pair_counts, decreasing = TRUE),
      weighted_gut = sort(weighted_gut_counts, decreasing = TRUE),
      weighted_oral = sort(weighted_oral_counts, decreasing = TRUE)
    )
  )
}


setwd(r4projects::get_project_wd())

gut_oral_interaction <- readRDS("1_code/gut_oral_microbiome/combined_results_with_interactions")

# Translated comment.
significant_metabolites <- select_significant_interaction_metabolites(
  gut_oral_interaction,
  min_r2 = 0.35,
  min_interaction_importance = 0.1,
  top_n = 50
)

# Translated comment.
interaction_plots <- plot_interaction_effects(significant_metabolites)
print(interaction_plots$combined)

# Translated comment.
interaction_details <- analyze_interaction_features(significant_metabolites)


### Translated comment.
gut_frequency <- data.frame(interaction_details$gut_frequency)
weighted_gut <- data.frame(interaction_details$weighted_gut)
weighted_gut$ASV <- rownames(weighted_gut)
colnames(weighted_gut) <- c("importance", "ASV")
gut_ASV_importance <- merge(gut_frequency, weighted_gut, by.x = "Var1", by.y =
                              "ASV")

oral_frequency <- data.frame(interaction_details$oral_frequency)
weighted_oral <- data.frame(interaction_details$weighted_oral)
weighted_oral$ASV <- rownames(weighted_oral)
colnames(weighted_oral) <- c("importance", "ASV")
oral_ASV_importance <- merge(oral_frequency,
                             weighted_oral,
                             by.x = "Var1",
                             by.y = "ASV")

### Translated comment.
all_interactions <- interaction_details$all_interactions

all_interactions <- merge(all_interactions,
                          metabolite_annotation[, c("variable_id", "HMDB.Name", "HMDB.Class")],
                          by.x = "metabolite",
                          by.y = "variable_id")

all_interactions_gut <- data.frame(table(all_interactions$gut_feature, all_interactions$HMDB.Class))

HMDB_Class <- c(
  "Benzene and substituted derivatives",
  "Carboxylic acids and derivatives",
  "Fatty Acyls",
  "Glycerophospholipids",
  "Organic sulfuric acids and derivatives",
  "Organonitrogen compounds",
  "Piperidines"
)

all_interactions_gut <- subset(all_interactions_gut, Var2 %in% HMDB_Class)

all_interactions_oral <- data.frame(table(all_interactions$oral_feature, all_interactions$HMDB.Class))

all_interactions_oral <- subset(all_interactions_oral, Var2 %in% HMDB_Class)


# Translated comment.

all_interactions_gut <- merge(all_interactions_gut, gut_ASV_importance, by =
                                "Var1")

colnames(all_interactions_gut) <- c("ASV", "HMDB.Class", "Freq", "Freq_Sum", "Importance")

gut_tax <- data.frame(gut_temp_object@variable_info)
all_interactions_gut <- merge(all_interactions_gut,
                              gut_tax,
                              by.x = "ASV",
                              by.y = "variable_id")

all_interactions_oral <- merge(all_interactions_oral, oral_ASV_importance, by =
                                 "Var1")

colnames(all_interactions_oral) <- c("ASV", "HMDB.Class", "Freq", "Freq_Sum", "Importance")

oral_tax <- data.frame(oral_temp_object@variable_info)
all_interactions_oral <- merge(all_interactions_oral,
                               oral_tax,
                               by.x = "ASV",
                               by.y = "variable_id")
all_interactions_oral <- subset(all_interactions_oral, !(ASV %in% c("OTU_68", "OTU_80", "OTU_729", "OTU_493")))


# Translated comment.
phyla <- c(
  "Firmicutes",
  "Proteobacteria",
  "Bacteroidetes",
  "Actinobacteria",
  "Cyanobacteria/Chloroplast",
  "Unclassified_Bacteria",
  "Fusobacteria",
  "Spirochaetes",
  "Tenericutes"
)

# Translated comment.
phylum_colors <- c(
  "Firmicutes" = "#fbb4ae",
  # Translated comment.
  "Proteobacteria" = "#ccebc5",
  # Translated comment.
  "Bacteroidetes" = "#b3cde3",
  # Translated comment.
  "Actinobacteria" = "#BCECE0",
  # Translated comment.
  "Cyanobacteria/Chloroplast" = "#7D5BA6",
  # Translated comment.
  "Unclassified_Bacteria" = "#8A89C0",
  # Translated comment.
  "Fusobacteria" = "#5762D5",
  # Translated comment.
  "Spirochaetes" = "#FC9E4F",
  # Translated comment.
  "Tenericutes" = "# translated comment
)


p1 <- ggplot(data = all_interactions_gut) +
  geom_tile(aes(x = "a", y = Genus, fill = Phylum), width = 0.5) +
  labs(x = "", y = "") +
  scale_fill_manual(values = phylum_colors) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    # panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(face = "italic", hjust = 1),
    legend.position = "none"
  )

p1
library(ggh4x)
p2 <- ggplot(data = all_interactions_gut) +
  geom_bar(aes(x = Freq, y = Genus, fill = HMDB.Class),
           stat = "identity",
           width = 0.6) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = c(
    "#ef6548",
    "#ffeda0",
    "#3d95d2",
    "#7d4b3c",
    "#007b7a",
    "#546672",
    "#de9db5"
  )) +
  labs(x = "", y = "") +
  guides(y.sec = guide_axis_manual(breaks = 1:36)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    # panel.border = element_blank(),
    axis.line.x.bottom = element_line(linewidth = 0.5),
    axis.line.y.right = element_line(linewidth = 0.5, linetype = 2),
    axis.ticks.y = element_blank(),
    axis.ticks.length.x = unit(5, "pt"),
    axis.text.y = element_blank(),
    legend.position = "bottom"
  )

p3 <- ggplot(data = all_interactions_oral) +
  geom_bar(aes(x = -1 * Freq, y = Genus, fill = HMDB.Class),
           stat = "identity",
           width = 0.6) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0)),
    breaks = c(-50, -40, -30, -20, -10, 0),
    labels = c(50, 40, 30, 20, 10, 0)
  ) +
  scale_fill_manual(values = c(
    "#ef6548",
    "#ffeda0",
    "#3d95d2",
    "#7d4b3c",
    "#007b7a",
    "#546672",
    "#de9db5"
  )) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    # panel.border = element_blank(),
    axis.line.x.bottom = element_line(linewidth = 0.5),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.length.x = unit(5, "pt"),
    legend.position = "none"
  )


p4 <- ggplot(data = all_interactions_oral) +
  geom_tile(aes(x = "a", y = Genus, fill = Phylum), width = 0.5) +
  labs(x = "", y = "") +
  scale_fill_manual(values = phylum_colors) +
  scale_y_discrete(position = "right") +  # translated comment
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    # panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y.left = element_blank(),
    # Translated comment.
    axis.text.y = element_text(face = "italic", hjust = 0),
    # Translated comment.
    legend.position = "bottom",
    legend.justification = "right"
  )

p4

p_combine <- p1 + p2 + p3 + p4 + plot_layout(nrow = 1, widths = c(0.25, 1.25, 1.25, 0.25))

p_combine

ggsave(
  p_combine,
  filename = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_4/figure_4c.pdf"
  ),
  width = 10,
  height = 6
)

# Translated comment.


plot_microbe_interaction_quantile <- function(gut_data,
                                              oral_data,
                                              metabolite_data,
                                              gut_feature,
                                              oral_feature,
                                              metabolite,
                                              n_quantiles = 3,
                                              use_log = FALSE) {
  # Translated comment.
  gut_samples <- colnames(gut_data)
  oral_samples <- colnames(oral_data)
  meta_samples <- colnames(metabolite_data)
  common_samples <- Reduce(intersect, list(gut_samples, oral_samples, meta_samples))
  
  # Translated comment.
  gut_values <- as.numeric(gut_data[gut_feature, common_samples])
  oral_values <- as.numeric(oral_data[oral_feature, common_samples])
  metabolite_values <- as.numeric(metabolite_data[metabolite, common_samples])
  
  # Translated comment.
  if (use_log) {
    # Translated comment.
    gut_values[gut_values == 0] <- min(gut_values[gut_values > 0]) / 2
    oral_values[oral_values == 0] <- min(oral_values[oral_values > 0]) / 2
    
    gut_values <- log10(gut_values)
    oral_values <- log10(oral_values)
  }
  
  # Translated comment.
  plot_data <- data.frame(gut = gut_values,
                          oral = oral_values,
                          metabolite = metabolite_values)
  
  # Translated comment.
  quantile_breaks <- quantile(oral_values, probs = seq(0, 1, length.out = n_quantiles + 1))
  plot_data$oral_quantile <- cut(
    oral_values,
    breaks = quantile_breaks,
    labels = paste0("Q", 1:n_quantiles),
    include.lowest = TRUE
  )
  
  # Translated comment.
  quantile_cors <- lapply(levels(plot_data$oral_quantile), function(q) {
    subset_data <- plot_data[plot_data$oral_quantile == q, ]
    cor_value <- cor.test(subset_data$gut, subset_data$metabolite, method = "spearman")
    data.frame(quantile = q, correlation = cor_value$estimate[[1]],p = cor_value$p.value*0.5)
  })
  quantile_cors <- do.call(rbind, quantile_cors)
  
  
  # Translated comment.
  p <- ggplot(plot_data, aes(x = gut, y = metabolite, color = oral_quantile)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, aes(group = oral_quantile)) +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = paste("Interaction effect of", gut_feature, "and", oral_feature),
      subtitle = paste0(
        "Correlation by quantile: ",
        paste(
          quantile_cors$quantile,
          ":",
          round(quantile_cors$correlation, 2),
          collapse = ", "
        ),
        "\n                           ",
        paste(
          quantile_cors$quantile,
          ":",
          round(quantile_cors$p, 2),
          collapse = ", "
        )
      ),
      x = paste0(gut_feature, if (use_log)
        " (log10)"
        else
          ""),
      y = paste0(metabolite),
      color = paste0(oral_feature, "\nQuantile")
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      legend.position = "none"
    ) + facet_wrap( ~ oral_quantile, scales = "free")
  
  return(p)
}


gut_data <- gut_temp_object@expression_data
rownames(gut_data) <- gut_tax$Genus

oral_data <- oral_temp_object@expression_data
rownames(oral_data) <- oral_tax$Genus

metabolomics_data <- metabolomics_temp_object@expression_data



a<-plot_microbe_interaction_quantile(
  gut_data = gut_data,
  oral_data = oral_data,
  metabolite_data= metabolomics_data,
  gut_feature = "Paraprevotella",
  oral_feature = "Enterococcus",
  metabolite = "M205T407_2_POS_HILIC",
  n_quantiles = 3,
  use_log = FALSE
) + ylab("L-Tryptophan")+xlim(c(-0.5,1.5))

a

ggsave(
  a,
  filename = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_4/figure_4d.pdf"
  ),
  width = 9,
  height = 4
)



b<-plot_microbe_interaction_quantile(
  gut_data = gut_data,
  oral_data = oral_data,
  metabolite_data= metabolomics_data,
  gut_feature = "Paraprevotella",
  oral_feature = "Enterococcus",
  metabolite = "M219T395_NEG_HILIC",
  n_quantiles = 3,
  use_log = FALSE
) + ylab("5-Hydroxy-L-tryptophan")+xlim(c(-0.5,1.5))



ggsave(
  b,
  filename = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_4/figure_4e.pdf"
  ),
  width = 9,
  height = 4
)





c <- plot_microbe_interaction_quantile(
  gut_data = gut_data,
  oral_data = oral_data,
  metabolite_data = metabolomics_data,
  gut_feature = "Mediterraneibacter",
  oral_feature = "Lachnoanaerobaculum",
  metabolite = "M178T111_NEG_RPLC",
  n_quantiles = 3,
  use_log = FALSE
) + ylab("Hippuric acid")

c
ggsave(
  c,
  filename = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_4/figure_s4c.pdf"
  ),
  width = 9,
  height = 4
)




d <- plot_microbe_interaction_quantile(
  gut_data = gut_data,
  oral_data = oral_data,
  metabolite_data = metabolomics_data,
  gut_feature = "Mediterraneibacter",
  oral_feature = "Lachnoanaerobaculum",
  metabolite = "M199T542_NEG_HILIC",
  n_quantiles = 3,
  use_log = FALSE
) + ylab("gamma-Glutamylalanine")

ggsave(
  d,
  filename = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_4/figure_s4d.pdf"
  ),
  width = 9,
  height = 4
)



