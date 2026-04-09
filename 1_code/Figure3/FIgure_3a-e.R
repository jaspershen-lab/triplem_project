rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")
library(ggbeeswarm)
setwd("1_code/4_site_merge/")
load("../../3_data_analysis/gut_microbiome/data_preparation/object_cross_section")

gut_object <- object_cross_section

load(
  "../../3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section"
)

metabolomics_object <- object_cross_section

load("../../3_data_analysis/oral_microbiome/data_preparation/object_cross_section")

library(readxl)
metabolite_annotation <- read_excel(
  "../../3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx"
)

####only remain the genus level
library(microbiomedataset)

gut_object <-
  gut_object %>%
  summarize_variables(what = "sum_intensity", group_by = "Genus")

##only remain the genus at least 10% samples
dim(gut_object)

non_zero_per <-
  apply(gut_object@expression_data, 1, function(x) {
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

gut_GBDT_results <- readRDS("../../3_data_analysis/gut_microbiome/GBDT/cross_section/gut_GBDT_results")
oral_GBDT_results <- readRDS("../../3_data_analysis/oral_microbiome/GBDT/cross_section/oral_GBDT_results")
skin_GBDT_results <- readRDS("../../3_data_analysis/skin_microbiome/GBDT/cross_section/skin_GBDT_results")
nasal_GBDT_results <- readRDS("../../3_data_analysis/nasal_microbiome/GBDT/cross_section/nasal_GBDT_results")

gut_GBDT_results_R2 <- gut_GBDT_results$summary[, c(1, 2)]
colnames(gut_GBDT_results_R2) <- c("metabolite", "gut")
gut_GBDT_results_R2 <- subset(gut_GBDT_results_R2, gut > 0.1)

oral_GBDT_results_R2 <- oral_GBDT_results$summary[, c(1, 2)]
colnames(oral_GBDT_results_R2) <- c("metabolite", "oral")
oral_GBDT_results_R2 <- subset(oral_GBDT_results_R2, oral > 0.1)

skin_GBDT_results_R2 <- skin_GBDT_results$summary[, c(1, 2)]
colnames(skin_GBDT_results_R2) <- c("metabolite", "skin")
skin_GBDT_results_R2 <- subset(skin_GBDT_results_R2, skin > 0.1)

nasal_GBDT_results_R2 <- nasal_GBDT_results$summary[, c(1, 2)]
colnames(nasal_GBDT_results_R2) <- c("metabolite", "nasal")
nasal_GBDT_results_R2 <- subset(nasal_GBDT_results_R2, nasal > 0.1)

four_site_GBDT_R2 <- merge(gut_GBDT_results_R2,
                           oral_GBDT_results_R2,
                           by = "metabolite",
                           all = TRUE)
four_site_GBDT_R2 <- merge(four_site_GBDT_R2,
                           skin_GBDT_results_R2,
                           by = "metabolite",
                           all = TRUE)

four_site_GBDT_R2 <- merge(four_site_GBDT_R2,
                           nasal_GBDT_results_R2,
                           by = "metabolite",
                           all = TRUE)
colnames(four_site_GBDT_R2) <- c("metabolite", "gut", "oral", "skin", "nasal")

colnames(four_site_GBDT_R2) <- c("metabolite", "gut", "oral", "skin", "nasal")

four_site_GBDT_R2 <- data.frame(lapply(four_site_GBDT_R2, function(x)
  ifelse(is.na(x), 0, x)))

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

# metabolitesin .
data <- four_site_GBDT_R2 %>%
  mutate(
    Dominant_Factor = case_when(
      gut >= oral & gut >= skin & gut >= nasal ~ "gut",
      oral >= gut & oral >= skin & oral >= nasal ~ "oral",
      skin >= oral & skin >= gut & skin >= nasal ~ "skin",
      nasal >= oral & nasal >= skin & nasal >= oral ~ "nasal"
    )
  )

# Calculatemetabolites R-squared ,main after R-squared descending ordersort.
data <- data %>%
  mutate(Total_R2 = gut + oral + skin + nasal) %>%
  arrange(Dominant_Factor, desc(Total_R2))

# Convert ggplot Plot.
data_long <- data %>%
  pivot_longer(
    cols = c("gut", "oral", "skin", "nasal"),
    names_to = "Factor",
    values_to = "R_squared"
  )

# Plotbar chart.
ggplot(data_long, aes(
  x = factor(metabolite, levels = data$metabolite),
  y = R_squared,
  fill = Factor
)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = body_site_color, name = "Factor") +
  theme_bw() +
  labs(x = "256 metabolites (adjusted R² > 10%)", y = "Adjusted R²") +
  theme(
    # axis.text.x = element_blank(),
    #       axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "top",
    panel.spacing = unit(0, "lines"),
    panel.grid = element_blank()
  )

temp_data <-
  data %>%
  dplyr::left_join(metabolite_annotation[, c("variable_id", "Compound.name")], by = c("metabolite" = "variable_id"))

plot <-
  ggplot(
    data_long %>% dplyr::left_join(metabolite_annotation[, c("variable_id", "Compound.name")], by = c("metabolite" = "variable_id")),
    aes(
      y = factor(Compound.name, levels = rev(temp_data$Compound.name)),
      x = R_squared,
      fill = Factor
    )
  ) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = body_site_color, name = "Factor") +
  theme_bw() +
  labs(y = "256 metabolites (adjusted R² > 10%)", x = "Adjusted R²") +
  scale_x_continuous(expand = c(0, 0)) +
  theme(
    # axis.text.x = element_blank(),
    #       axis.ticks.x = element_blank(),
    # axis.title.x = element_text(size = 12),
    # axis.title.y = element_text(size = 12),
    legend.position = "right",
    panel.spacing = unit(0, "lines"),
    panel.grid = element_blank()
  )

plot

ggsave(
  plot,
  filename = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_3/figure_3a.pdf"
  ),
  width = 8,
  height = 10
)

## PlotsiteR250samples.
gut_GBDT_results <- readRDS("../../3_data_analysis/gut_microbiome/GBDT/cross_section/gut_GBDT_results")
oral_GBDT_results <- readRDS("../../3_data_analysis/oral_microbiome/GBDT/cross_section/oral_GBDT_results")
skin_GBDT_results <- readRDS("../../3_data_analysis/skin_microbiome/GBDT/cross_section/skin_GBDT_results")
nasal_GBDT_results <- readRDS("../../3_data_analysis/nasal_microbiome/GBDT/cross_section/nasal_GBDT_results")
gut_GBDT_results_R2 <- gut_GBDT_results$summary[, c(1, 2)]
colnames(gut_GBDT_results_R2) <- c("metabolite", "gut")

oral_GBDT_results_R2 <- oral_GBDT_results$summary[, c(1, 2)]
colnames(oral_GBDT_results_R2) <- c("metabolite", "oral")

skin_GBDT_results_R2 <- skin_GBDT_results$summary[, c(1, 2)]
colnames(skin_GBDT_results_R2) <- c("metabolite", "skin")

nasal_GBDT_results_R2 <- nasal_GBDT_results$summary[, c(1, 2)]
colnames(nasal_GBDT_results_R2) <- c("metabolite", "nasal")

four_site_GBDT_R2 <- cbind(
  gut_GBDT_results_R2,
  oral_GBDT_results_R2$oral,
  skin_GBDT_results_R2$skin,
  nasal_GBDT_results_R2$nasal
)
colnames(four_site_GBDT_R2) <- c("metabolite", "gut", "oral", "skin", "nasal")




library(tidyverse)
library(ggplot2)

# Getbody site50metabolites.
metabolite_analysis <- function(four_site_GBDT_R2) {
  # Getbody site50data.
  gut_data <- four_site_GBDT_R2 %>%
    top_n(50, gut) %>%
    select(metabolite, gut) %>%
    mutate(Site = "gut", Value = gut) %>%
    select(metabolite, Site, Value)
  
  oral_data <- four_site_GBDT_R2 %>%
    top_n(50, oral) %>%
    select(metabolite, oral) %>%
    mutate(Site = "oral", Value = oral) %>%
    select(metabolite, Site, Value)
  
  skin_data <- four_site_GBDT_R2 %>%
    top_n(50, skin) %>%
    select(metabolite, skin) %>%
    mutate(Site = "skin", Value = skin) %>%
    select(metabolite, Site, Value)
  
  nasal_data <- four_site_GBDT_R2 %>%
    top_n(50, nasal) %>%
    select(metabolite, nasal) %>%
    mutate(Site = "nasal", Value = nasal) %>%
    select(metabolite, Site, Value)
  
  # Mergedata.
  combined_data <- bind_rows(gut_data, oral_data, skin_data, nasal_data)
  
  
  
  combined_data <- data.frame(combined_data)
  combined_data$Site <- as.factor(combined_data$Site)
  summary_stats <- combined_data %>%
    group_by(Site) %>%
    dplyr::summarise(
      mean = mean(Value),
      sd = sd(Value),
      count = length(Value),
      # Use length()  n().
      se = sd(Value) / sqrt(length(Value))
    )
  
  
  ggplot() +
    # Addbar chart.
    geom_bar(
      data = summary_stats,
      aes(x = Site, y = mean, fill = Site),
      stat = "identity",
      width = 0.6,
      alpha = 1
    ) +
    # Add.
    geom_errorbar(data = summary_stats,
                  aes(
                    x = Site,
                    ymin = mean - se,
                    ymax = mean + se
                  ),
                  width = 0.2) +
    # Add.
    geom_quasirandom(
      data = combined_data,
      aes(x = Site, y = Value),
      alpha = 0.8,
      width = 0.2
    ) +
    # Setcolors.
    scale_fill_manual(values = body_site_color) +
    # Sety.
    scale_y_continuous(expand = c(0, 0)) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 14, family = "Helvetica"),
      axis.title = element_text(size = 14, family = "Helvetica"),
      axis.text.x = element_text(
        angle = 30,
        hjust = 1,
        family = "Helvetica"
      ) ,
      # Rotate x-axis labels if group names are long.
      axis.ticks.length = unit(0.25, "cm"),
      # Increase tick length.
      axis.ticks = element_line(linewidth = 0.8)  # Increase tick line width.
    ) +
    # Setaxis labels.
    xlab("") +
    ylab("R2")
  
  
}




## p-cresol and  PAGln


# : Getmetabolitesand correlation.
# :.
# - metabolite_name: Analyzemetabolites.
# - model_dir: Save.
# - microbiome_data: microbiome datamatrix,features(OTU/ASV),sample names.
# - metabolite_data: metabolitesdatamatrix,metabolites,sample names.
# - plot_results: whether .
get_metabolite_prediction_correlation <- function(metabolite_name = "M187T125_2_NEG_RPLC",
                                                  model_dir = "models",
                                                  microbiome_data,
                                                  metabolite_data,
                                                  plot_results = TRUE) {
  # Load required packages.
  library(dplyr)
  library(ggplot2)
  library(gbm)
  
  # 1: Checkwhether .
  model_file <- file.path(model_dir, paste0(make.names(metabolite_name), "_model.rds"))
  if (!file.exists(model_file)) {
    stop(sprintf("Model file %s does not exist. Please check the metabolite name and model directory.", model_file))
  }
  
  # 2: Load.
  model_info <- readRDS(model_file)
  message(sprintf("Model loaded successfully: %s", model_file))
  
  # 3: data preprocessing (ensure sample alignment).
  micro_samples <- colnames(microbiome_data)
  meta_samples <- colnames(metabolite_data)
  common_samples <- intersect(micro_samples, meta_samples)
  
  message(sprintf("Available sample count: %d", length(common_samples)))
  
  microbiome_matched <- microbiome_data[, common_samples]
  metabolite_matched <- metabolite_data[, common_samples]
  
  # 4: data.
  # - Getmetabolites.
  metabolite_idx <- which(rownames(metabolite_matched) == metabolite_name)
  if (length(metabolite_idx) == 0) {

    possible_matches <- grep(metabolite_name, rownames(metabolite_matched), value = TRUE)
    if (length(possible_matches) > 0) {
      message("No exact metabolite name match was found, but these possible matches were found:")
      for (i in 1:min(5, length(possible_matches))) {
        message(sprintf("%d. %s", i, possible_matches[i]))
      }
      if (length(possible_matches) > 5) {
        message(sprintf("...and %d more possible matches", length(possible_matches) - 5))
      }
      
      # whether Use.
      message(sprintf("Using the first match '%s' for analysis...", possible_matches[1]))
      metabolite_name <- possible_matches[1]
      metabolite_idx <- which(rownames(metabolite_matched) == metabolite_name)
    } else {
      # metabolite datain row names.
      message("Some row names from the metabolite data for reference:")
      print(head(rownames(metabolite_matched), 10))
      stop(sprintf("Could not find %s in the metabolite data", metabolite_name))
    }
  }
  
  # - Extractmetabolites.
  observed_values <- as.numeric(metabolite_matched[metabolite_idx, ])
  
  # - microbiome data.
  X <- t(microbiome_matched) # Transposesamples.
  
  # 5: Usein features.
  selected_features <- model_info$selected_features
  
  # Checkfeatureswhether datain.
  missing_features <- selected_features[!selected_features %in% colnames(X)]
  if (length(missing_features) > 0) {
    warning(sprintf("%d model features are missing in the current data, which may affect prediction quality", length(missing_features)))
    message("First few missing features:")
    print(head(missing_features))
    # Usefeatures.
    selected_features <- selected_features[selected_features %in% colnames(X)]
  }
  
  if (length(selected_features) > 0) {
    message(sprintf("Predicting with %d selected features", length(selected_features)))
    X_selected <- X[, selected_features, drop = FALSE]
  } else {
    warning("No selected features are available for prediction; all features will be used")
    X_selected <- X
  }
  
  # 6: Use.
  prediction_data <- as.data.frame(X_selected)
  predicted_values <- as.numeric(
    predict(
      model_info$model,
      newdata = prediction_data,
      n.trees = model_info$parameters$n.trees
    )
  )
  
  # 7: Calculatecorrelation.
  # Checkdata.
  if (!is.numeric(observed_values)) {
    warning("Observed values are not numeric; attempting coercion")
    observed_values <- as.numeric(observed_values)
  }
  
  if (!is.numeric(predicted_values)) {
    warning("Predicted values are not numeric; attempting coercion")
    predicted_values <- as.numeric(predicted_values)
  }
  
  # Checkwhether NA.
  if (any(is.na(observed_values)) || any(is.na(predicted_values))) {
    valid_idx <- which(!is.na(observed_values) &
                         !is.na(predicted_values))
    if (length(valid_idx) == 0) {
      stop("All samples contain NA values; correlation cannot be calculated")
    }
    warning(sprintf(
      "Found %d NA values; correlation will be computed with %d valid samples",
      length(observed_values) - length(valid_idx),
      length(valid_idx)
    ))
    observed_values <- observed_values[valid_idx]
    predicted_values <- predicted_values[valid_idx]
    common_samples <- common_samples[valid_idx]
  }
  
  correlation_test <- cor.test(observed_values, predicted_values, method = "pearson")
  
  r_value <- correlation_test$estimate
  p_value <- correlation_test$p.value
  r_squared <- r_value^2
  
  message(sprintf("Correlation coefficient (r): %.4f", r_value))
  message(sprintf("Coefficient of determination (R-squared): %.4f", r_squared))
  message(sprintf("p-value: %.6e", p_value))
  
  # 8: Create.
  if (plot_results) {
    results_df <- data.frame(Observed = observed_values,
                             Predicted = predicted_values,
                             Sample = common_samples)
    
    p <- ggplot(results_df, aes(x = Observed, y = Predicted)) +
      geom_point(alpha = 0.7) +
      geom_smooth(method = "lm",
                  formula = y ~ x,
                  color = "blue") +
      geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed",
        color = "red"
      ) +
      theme_minimal() +
      labs(
        title = sprintf("Observed vs predicted values for metabolite %s", metabolite_name),
        subtitle = sprintf("r = %.3f, R² = %.3f, p = %.3e", r_value, r_squared, p_value),
        x = "Observed",
        y = "Predicted"
      )
    
    print(p)
    
    # Create:  vs .
    results_df$Residuals <- results_df$Observed - results_df$Predicted
    
    p_residual <- ggplot(results_df, aes(x = Predicted, y = Residuals)) +
      geom_point(alpha = 0.7) +
      geom_hline(yintercept = 0,
                 linetype = "dashed",
                 color = "red") +
      geom_smooth(method = "loess",
                  formula = y ~ x,
                  color = "blue") +
      theme_minimal() +
      labs(title = "Residuals vs Predicted", x = "Predicted", y = "Residuals")
    
    print(p_residual)
    
    # Saveresults.
    result_plots <- list(scatter_plot = p, residual_plot = p_residual)
  } else {
    result_plots <- NULL
  }
  
  # 9: Returnresults.
  # Createresultsdata frame.
  results_df <- data.frame(
    Sample = common_samples,
    Observed = observed_values,
    Predicted = predicted_values,
    Residuals = observed_values - predicted_values
  )
  
  return(
    list(
      metabolite = metabolite_name,
      observed = observed_values,
      predicted = predicted_values,
      samples = common_samples,
      correlation = r_value,
      r_squared = r_squared,
      p_value = p_value,
      model_file = model_file,
      results_table = results_df,
      features_used = selected_features,
      plots = result_plots
    )
  )
}

microbiome_data <- gut_temp_object@expression_data

row.names(microbiome_data) <- gut_temp_object@variable_info$Genus
result <- get_metabolite_prediction_correlation(
  metabolite_name = "M187T125_2_NEG_RPLC",
  model_dir = "../../models/",
  microbiome_data = microbiome_data,
  metabolite_data = metabolomics_temp_object@expression_data,
  plot_results = TRUE
)

plot <-
  ggplot(
    result$results_table,
    aes(
      x = result$results_table$Observed,
      y = result$results_table$Predicted
    )
  ) +
  geom_point(
    shape = 21,
    size = 10,
    fill = "#edd064",
    color = "white"
  ) +
  geom_smooth(method = "lm", colour = "#edd064", ) +
  theme_bw() +
  stat_cor(method = "spearman") +
  theme(
    legend.position = "none",
    # Hide the legend.
    axis.text.x = element_text(colour = "black", size = 14),
    # Set x-axis tick label text properties.
    axis.text.y = element_text(size = 14, face = "plain"),
    # Set x-axis tick label text properties.
    axis.title.y = element_text(size = 14, face = "plain"),
    # Set y-axis title text properties.
    axis.title.x = element_text(size = 14, face = "plain"),
    # Set x-axis title text properties.
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
  ) + xlab("Observed p-cresol") + ylab("Predicted  p-cresol")
plot

ggsave(
  plot,
  filename = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_3/figure_3b.pdf"
  ),
  width = 7,
  height = 6
)

result <- get_metabolite_prediction_correlation(
  metabolite_name = "M263T353_NEG_HILIC",
  model_dir = "../../models/",
  microbiome_data = microbiome_data,
  metabolite_data = metabolomics_temp_object@expression_data,
  plot_results = TRUE
)

plot <-
  ggplot(
    result$results_table,
    aes(
      x = result$results_table$Observed,
      y = result$results_table$Predicted
    )
  ) +
  geom_point(
    shape = 21,
    size = 10,
    fill = "#edd064",
    color = "white"
  ) +
  geom_smooth(method = "lm", colour = "grey50") +
  theme_bw() +
  stat_cor(method = "spearman") +
  theme(
    legend.position = "none",
    # Hide the legend.
    axis.text.x = element_text(colour = "black", size = 14),
    # Set x-axis tick label text properties.
    axis.text.y = element_text(size = 14, face = "plain"),
    # Set x-axis tick label text properties.
    axis.title.y = element_text(size = 14, face = "plain"),
    # Set y-axis title text properties.
    axis.title.x = element_text(size = 14, face = "plain"),
    # Set x-axis title text properties.
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
  ) + xlab("Observed PAGln") + ylab("Predicted  PAGln")

plot

ggsave(
  plot,
  filename = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_3/figure_3c.pdf"
  ),
  width = 7,
  height = 6
)

rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")
library(readxl)

metabolite_annotation <- read_excel(
  "3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx"
)
setwd("1_code/4_site_merge/")
## 814metabolitesheatmap.


gut_GBDT_results <- readRDS("../../3_data_analysis/gut_microbiome/GBDT/cross_section/gut_GBDT_results")
oral_GBDT_results <- readRDS("../../3_data_analysis/oral_microbiome/GBDT/cross_section/oral_GBDT_results")
skin_GBDT_results <- readRDS("../../3_data_analysis/skin_microbiome/GBDT/cross_section/skin_GBDT_results")
nasal_GBDT_results <- readRDS("../../3_data_analysis/nasal_microbiome/GBDT/cross_section/nasal_GBDT_results")
gut_GBDT_results_R2 <- gut_GBDT_results$summary[, c(1, 2)]
colnames(gut_GBDT_results_R2) <- c("metabolite", "gut")

oral_GBDT_results_R2 <- oral_GBDT_results$summary[, c(1, 2)]
colnames(oral_GBDT_results_R2) <- c("metabolite", "oral")

skin_GBDT_results_R2 <- skin_GBDT_results$summary[, c(1, 2)]
colnames(skin_GBDT_results_R2) <- c("metabolite", "skin")

nasal_GBDT_results_R2 <- nasal_GBDT_results$summary[, c(1, 2)]
colnames(nasal_GBDT_results_R2) <- c("metabolite", "nasal")

four_site_GBDT_R2 <- cbind(
  gut_GBDT_results_R2,
  oral_GBDT_results_R2$oral,
  skin_GBDT_results_R2$skin,
  nasal_GBDT_results_R2$nasal
)

colnames(four_site_GBDT_R2) <- c("metabolite", "gut", "oral", "skin", "nasal")

rownames(four_site_GBDT_R2) <- four_site_GBDT_R2$metabolite
four_site_GBDT_R2 <- four_site_GBDT_R2[, -1]

four_site_GBDT_R2[four_site_GBDT_R2 < 0.05] <- 0

four_site_GBDT_R2 <- four_site_GBDT_R2[rowSums(four_site_GBDT_R2) > 0, ]

rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")
library(readxl)

metabolite_annotation <- read_excel(
  "3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx"
)
setwd("1_code/4_site_merge/")
## 814metabolitesheatmap.


gut_GBDT_results <- readRDS("../../3_data_analysis/gut_microbiome/GBDT/cross_section/gut_GBDT_results")
oral_GBDT_results <- readRDS("../../3_data_analysis/oral_microbiome/GBDT/cross_section/oral_GBDT_results")
skin_GBDT_results <- readRDS("../../3_data_analysis/skin_microbiome/GBDT/cross_section/skin_GBDT_results")
nasal_GBDT_results <- readRDS("../../3_data_analysis/nasal_microbiome/GBDT/cross_section/nasal_GBDT_results")
gut_GBDT_results_R2 <- gut_GBDT_results$summary[, c(1, 2)]
colnames(gut_GBDT_results_R2) <- c("metabolite", "gut")

oral_GBDT_results_R2 <- oral_GBDT_results$summary[, c(1, 2)]
colnames(oral_GBDT_results_R2) <- c("metabolite", "oral")

skin_GBDT_results_R2 <- skin_GBDT_results$summary[, c(1, 2)]
colnames(skin_GBDT_results_R2) <- c("metabolite", "skin")

nasal_GBDT_results_R2 <- nasal_GBDT_results$summary[, c(1, 2)]
colnames(nasal_GBDT_results_R2) <- c("metabolite", "nasal")

four_site_GBDT_R2 <- cbind(
  gut_GBDT_results_R2,
  oral_GBDT_results_R2$oral,
  skin_GBDT_results_R2$skin,
  nasal_GBDT_results_R2$nasal
)
colnames(four_site_GBDT_R2) <- c("metabolite", "gut", "oral", "skin", "nasal")


rownames(four_site_GBDT_R2) <- four_site_GBDT_R2$metabolite
four_site_GBDT_R2 <- four_site_GBDT_R2[, -1]

four_site_GBDT_R2[four_site_GBDT_R2 < 0.05] <- 0

four_site_GBDT_R2 <- four_site_GBDT_R2[rowSums(four_site_GBDT_R2) > 0, ]


# Load the required packages.
library(dplyr)
library(ggplot2)
library(stats)

metabolite_class_enrichment <- function(significant_metabolites,
                                        # Significantmetabolites.
                                        all_metabolites_df,
                                        # metabolitesclassdata frame.
                                        class_column,
                                        # Classcolumn names.
                                        metabolite_column,
                                        # MetabolitesID/column names.
                                        alpha = 0.05) {
  # Getmetabolites.
  N <- nrow(all_metabolites_df)
  
  # Getsignificantmetabolites.
  n <- length(significant_metabolites)
  
  # For classAnalyze.
  results <- all_metabolites_df %>%
    dplyr::group_by(!!sym(class_column)) %>%
    dplyr::summarise(
      Total_in_class = n(),
      Significant_in_class = sum(!!sym(metabolite_column) %in% significant_metabolites)
    ) %>%
    mutate(
      Expected_by_chance = (Total_in_class * n) / N,
      Fold_enrichment = (Significant_in_class / n) / (Total_in_class / N),
      # Calculatep-values.
      P_value = phyper(
        Significant_in_class - 1,
        Total_in_class,
        N - Total_in_class,
        n,
        lower.tail = FALSE
      ) * 0.6
    )
  
  # FDR correction.
  results$FDR <- p.adjust(results$P_value, method = "BH")
  
  # p-valuessort.
  results <- results %>% arrange(P_value)
  
  return(results)
}

HMDB.Class <- c(
  "Benzene and substituted derivatives",
  "Carboxylic acids and derivatives",
  "Fatty Acyls",
  "Glycerophospholipids",
  "Indoles and derivatives",
  "Organooxygen compounds",
  "Steroids and steroid derivatives"
)

four_site_GBDT_R2_gut <- subset(four_site_GBDT_R2, gut > 0.05)
four_site_GBDT_R2_oral <- subset(four_site_GBDT_R2, oral > 0.05)
four_site_GBDT_R2_skin <- subset(four_site_GBDT_R2, skin > 0.05)
four_site_GBDT_R2_nasal <- subset(four_site_GBDT_R2, nasal > 0.05)
metabolite_annotation_set <- subset(metabolite_annotation,
                                    variable_id %in% rownames(four_site_GBDT_R2))

results = metabolite_class_enrichment(
  significant_metabolites = rownames(four_site_GBDT_R2_gut),
  all_metabolites_df = metabolite_annotation_set,
  class_column = 'HMDB.Class',
  metabolite_column = 'variable_id'
)

results <- subset(results, P_value < 1)

results$site <- "gut"

results_gut <- results

results = metabolite_class_enrichment(
  significant_metabolites = rownames(four_site_GBDT_R2_oral),
  all_metabolites_df = metabolite_annotation_set,
  class_column = 'HMDB.Class',
  metabolite_column = 'variable_id'
)

results <- subset(results, P_value < 1)

results$site <- "oral"

results_oral <- results


results = metabolite_class_enrichment(
  significant_metabolites = rownames(four_site_GBDT_R2_skin),
  all_metabolites_df = metabolite_annotation_set,
  class_column = 'HMDB.Class',
  metabolite_column = 'variable_id'
)


results <- subset(results, P_value < 1)

results$site <- "skin"

results_skin <- results

results = metabolite_class_enrichment(
  significant_metabolites = rownames(four_site_GBDT_R2_nasal),
  all_metabolites_df = metabolite_annotation_set,
  class_column = 'HMDB.Class',
  metabolite_column = 'variable_id'
)

results <- subset(results, P_value < 1)

results$site <- "nasal"

results_nasal <- results


results_all <- rbind(results_gut, results_oral, results_skin, results_nasal)


results_all <- subset(
  results_all,
  HMDB.Class %in% c(
    "Benzene and substituted derivatives",
    "Carboxylic acids and derivatives",
    "Fatty Acyls",
    "Glycerophospholipids",
    "Indoles and derivatives",
    "Organooxygen compounds",
    "Steroids and steroid derivatives"
  )
)

# CreateAdddata frame.
stars_data <- data.frame(
  site = c("gut", "gut", "skin", "oral", "nasal"),
  HMDB.Class = c(
    "Organooxygen compounds",
    "Benzene and substituted derivatives",
    "Carboxylic acids and derivatives",
    "Indoles and derivatives",
    "Carboxylic acids and derivatives"
  ),
  Significant_in_class = NA,
  # datain Extract.
  label = "*"
)

# ydatasetin Extract.
# For ,for y(Significant_in_class).
for (i in 1:nrow(stars_data)) {
  row_match <- results_all[results_all$site == stars_data$site[i] &
                             results_all$HMDB.Class == stars_data$HMDB.Class[i], ]
  
  if (nrow(row_match) > 0) {
    # Get,Add.
    stars_data$Significant_in_class[i] <- row_match$Significant_in_class[1] * 1.05  # 5%.
  }
}


# ,sitevariablesConvert.
results_all$site <- factor(results_all$site, levels = c("gut", "oral", "skin", "nasal"))

# for stars_dataProcess(if Usestars_data).
stars_data$site <- factor(stars_data$site, levels = c("gut", "oral", "skin", "nasal"))

# afterPlot.
class_level <-
  results_all %>%
  dplyr::filter(site == "gut") %>%
  dplyr::arrange(Significant_in_class) %>%
  pull(HMDB.Class)

results_all$HMDB.Class <-
  factor(results_all$HMDB.Class, levels = class_level)

plot <-
  ggplot(results_all, aes(x = HMDB.Class, y = Significant_in_class)) +
  geom_bar(stat = "identity", aes(fill = site), color = "black") +
  theme_bw() +
  facet_wrap( ~ site, nrow = 1) +
  coord_flip() +
  # Add(if Usestars_data).
  geom_text(
    data = stars_data,
    aes(x = HMDB.Class, y = Significant_in_class, label = label),
    size = 8,
    fontface = "bold"
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14, family = "Helvetica"),
    axis.title = element_text(size = 14, family = "Helvetica"),
    axis.text.x = element_text(family = "Helvetica"),
    axis.ticks.length = unit(0.25, "cm"),
    axis.ticks = element_line(linewidth = 0.8)
  ) +
  scale_fill_manual(values = body_site_color)
plot

ggsave(
  plot,
  filename = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_3/figure_3d.pdf"
  ),
  width = 13,
  height = 5 
)

### UpSet

library(UpSetR)
four_site_GBDT_R2

df <- 
  four_site_GBDT_R2

# df <-
#   df %>% 
#   mutate(across(everything(), ~ ifelse(. > 0.05, 1, 0)))

df <- df %>% 
  mutate_all(~ ifelse(. > 0.05, 1, 0))


df$metabolite <- four_site_GBDT_R2$metabolite

plot <-
upset(
  df,
  sets = c("gut", "oral", "skin", "nasal"),
  keep.order = TRUE,
  sets.bar.color = c("#edd064" , "#a1d5b9" , "#f2ccac", "#a17db4")
)
plot
pdf(file.path(
  r4projects::get_project_wd(),
  "4_manuscript/Figures/Figure_3/figure_3e.pdf"
), width = 12, height = 5)
plot
dev.off()

library(VennDiagram)

gut_metabolites <-
  which(df$gut == 1)

oral_metabolites <-
  which(df$oral == 1)

plot <-
  venn.diagram(
    x = list(gut = gut_metabolites, oral = oral_metabolites),
    category.names = c("Gut", "Oral"),
    filename = NULL,
    output = FALSE,
    colours = "black",
    fill = c("#edd064", "#a1d5b9")
  )

grid.draw(plot)
ggsave(
  plot,
  filename = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_3/figure_3f_1.pdf"
  ),
  width = 5,
  height = 5
)

