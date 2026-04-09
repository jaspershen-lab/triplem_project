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


# Translated comment.
# Translated comment.
if (!requireNamespace("mediation", quietly = TRUE))
  install.packages("mediation")
if (!requireNamespace("tidyverse", quietly = TRUE))
  install.packages("tidyverse")

library(mediation)
library(tidyverse)
library(readr)

# Translated comment.
common_samples <- Reduce(intersect, list(
  colnames(gut_data),
  colnames(oral_data),
  colnames(metabolome_data)
))
gut_data <- gut_data[, common_samples]
oral_data <- oral_data[, common_samples]
metabolome_data <- metabolome_data[, common_samples]

# Translated comment.
gut_data_t <- t(gut_data)
oral_data_t <- t(oral_data)
metabolome_data_t <- t(metabolome_data)

# Translated comment.
calculate_correlations <- function(data1, data2) {
  result <- matrix(NA, nrow = ncol(data1), ncol = ncol(data2))
  rownames(result) <- colnames(data1)
  colnames(result) <- colnames(data2)
  p_values <- result
  
  for (i in 1:ncol(data1)) {
    for (j in 1:ncol(data2)) {
      cor_test <- cor.test(data1[, i], data2[, j], method = "spearman")
      result[i, j] <- cor_test$estimate
      p_values[i, j] <- cor_test$p.value
    }
  }
  
  # Translated comment.
  fdr_values <- matrix(p.adjust(p_values, method = "BH"), nrow = nrow(p_values))
  rownames(fdr_values) <- rownames(p_values)
  colnames(fdr_values) <- colnames(p_values)
  
  return(list(cor = result, p = p_values, fdr = fdr_values))
}

# Translated comment.
gut_metabolome_cor <- calculate_correlations(gut_data_t, metabolome_data_t)
gut_metabolome_associations <- which(gut_metabolome_cor$p < 0.05, arr.ind = TRUE)

oral_metabolome_cor <- calculate_correlations(oral_data_t, metabolome_data_t)
oral_metabolome_associations <- which(oral_metabolome_cor$p < 0.05, arr.ind = TRUE)

# Translated comment.
oral_gut_cor <- calculate_correlations(oral_data_t, gut_data_t)
oral_gut_associations <- which(oral_gut_cor$p < 0.05, arr.ind = TRUE)

# Translated comment.
bidirectional_mediation_results <- data.frame(
  direction = character(),
  oral_feature = character(),
  gut_feature = character(),
  metabolite = character(),
  ACME = numeric(),
  ACME_p = numeric(),
  prop_mediated = numeric(),
  stringsAsFactors = FALSE
)

# Translated comment.

# Translated comment.
for (i in 1:nrow(oral_gut_associations)) {
  oral_idx <- oral_gut_associations[i, "row"]
  gut_idx <- oral_gut_associations[i, "col"]
  
  oral_feature <- rownames(oral_gut_cor$cor)[oral_idx]
  gut_feature <- colnames(oral_gut_cor$cor)[gut_idx]
  
  # Translated comment.
  gut_metabolite_associations <- which(gut_metabolome_cor$fdr[gut_idx, ] < 0.05, arr.ind = TRUE)
  
  if (length(gut_metabolite_associations) > 0) {
    for (j in 1:length(gut_metabolite_associations)) {
      met_idx <- gut_metabolite_associations[j]
      metabolite <- colnames(gut_metabolome_cor$cor)[met_idx]
      
      # Translated comment.
      med_data <- data.frame(oral = oral_data_t[, oral_feature],
                             gut = gut_data_t[, gut_feature],
                             metabolite = metabolome_data_t[, metabolite])
      
      # Translated comment.
      med_model <- lm(gut ~ oral, data = med_data)
      out_model <- lm(metabolite ~ oral + gut + oral:gut, data = med_data)
      
      # Translated comment.
      med_result <- mediate(
        med_model,
        out_model,
        treat = "oral",
        mediator = "gut",
        boot = TRUE,
        sims = 100
      )
      
      # Translated comment.
      result_row <- data.frame(
        direction = "oral->gut->metabolite",
        oral_feature = oral_feature,
        gut_feature = gut_feature,
        metabolite = metabolite,
        ACME = med_result$d1,
        ACME_p = med_result$d1.p,
        prop_mediated = med_result$n1,
        stringsAsFactors = FALSE
      )
      
      bidirectional_mediation_results <- rbind(bidirectional_mediation_results, result_row)
    }
  }
}

# Translated comment.
for (i in 1:nrow(oral_gut_associations)) {
  oral_idx <- oral_gut_associations[i, "row"]
  gut_idx <- oral_gut_associations[i, "col"]
  
  oral_feature <- rownames(oral_gut_cor$cor)[oral_idx]
  gut_feature <- colnames(oral_gut_cor$cor)[gut_idx]
  
  # Translated comment.
  oral_metabolite_associations <- which(oral_metabolome_cor$fdr[oral_idx, ] < 0.05, arr.ind = TRUE)
  
  if (length(oral_metabolite_associations) > 0) {
    for (j in 1:length(oral_metabolite_associations)) {
      met_idx <- oral_metabolite_associations[j]
      metabolite <- colnames(oral_metabolome_cor$cor)[met_idx]
      
      # Translated comment.
      med_data <- data.frame(gut = gut_data_t[, gut_feature],
                             oral = oral_data_t[, oral_feature],
                             metabolite = metabolome_data_t[, metabolite])
      
      # Translated comment.
      med_model <- lm(oral ~ gut, data = med_data)
      out_model <- lm(metabolite ~ gut + oral + gut:oral, data = med_data)
      
      # Translated comment.
      med_result <- mediate(
        med_model,
        out_model,
        treat = "gut",
        mediator = "oral",
        boot = TRUE,
        sims = 100
      )
      
      # Translated comment.
      result_row <- data.frame(
        direction = "gut->oral->metabolite",
        oral_feature = oral_feature,
        gut_feature = gut_feature,
        metabolite = metabolite,
        ACME = med_result$d1,
        ACME_p = med_result$d1.p,
        prop_mediated = med_result$n1,
        stringsAsFactors = FALSE
      )
      
      bidirectional_mediation_results <- rbind(bidirectional_mediation_results, result_row)
    }
  }
}


# Translated comment.
if (nrow(bidirectional_mediation_results) > 0) {
  bidirectional_mediation_results$ACME_fdr <-
    p.adjust(bidirectional_mediation_results$ACME_p, method = "BH")
}

# Translated comment.

bidirectional_mediation_results_sig <- subset(bidirectional_mediation_results, ACME_p <
                                                0.1)


oral_tax <- oral_temp_object@variable_info
gut_tax <- gut_temp_object@variable_info
bidirectional_mediation_results_sig <- merge(bidirectional_mediation_results_sig,
                                             oral_tax[, c("Genus", "variable_id")],
                                             by.x = "oral_feature",
                                             by.y = "variable_id")

bidirectional_mediation_results_sig <- merge(bidirectional_mediation_results_sig,
                                             gut_tax[, c("Genus", "variable_id")],
                                             by.x = "gut_feature",
                                             by.y = "variable_id")

metabolite_annotation <- read_excel(
  "3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx"
)


bidirectional_mediation_results_sig <- merge(
  bidirectional_mediation_results_sig,
  metabolite_annotation[, c("HMDB.Name", "variable_id")],
  by.x = "metabolite",
  by.y = "variable_id"
)



# Translated comment.
direction_counts <-
  table(bidirectional_mediation_results_sig$direction)


plot <-
  ggplot(data = as.data.frame(direction_counts), aes(x = Var1, y = Freq, fill =
                                                       Var1)) +
  geom_bar(stat = "identity") +
  labs(x = "Direction", y = "Counts") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 14
    ),
    axis.text.y = element_text(size = 12)
  ) + scale_fill_manual(values = c("#Edd064", "#a1d5b9"))

plot

ggsave(
  plot,
  filename = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_5/figure_5b.pdf"
  ),
  width = 6,
  height = 6
)

bidirectional_mediation_results_sig_oral_gut_metabolite <- subset(
  bidirectional_mediation_results_sig,
  bidirectional_mediation_results_sig$direction == "oral->gut->metabolite"
)

variable_info <-
  extract_variable_info(object_cross_section)

bidirectional_mediation_results_sig_oral_gut_metabolite$HMDB.Name

which(bidirectional_mediation_results_sig_oral_gut_metabolite$HMDB.Name == "NA")

match(
  bidirectional_mediation_results_sig_oral_gut_metabolite$metabolite[30],
  variable_info$variable_id
)

bidirectional_mediation_results_sig_oral_gut_metabolite$HMDB.Name[30] <-
  variable_info$Compound.name[612]

sankey_data <- bidirectional_mediation_results_sig_oral_gut_metabolite[, 9:11]

colnames(sankey_data) <- c("Oral", "Gut", "Metabolite")

library(ggsankey)

df2 <- sankey_data %>%
  make_long(Oral, Gut, Metabolite)

plot <-
  ggplot(
    df2,
    aes(
      x = x,
      next_x = next_x,
      node = node,
      next_node = next_node,
      label = node,
      fill = factor(node)
    ),
    show.legend = FALSE
  ) +
  geom_sankey(flow.alpha = .6, node.color = "black") +
  geom_sankey_label(size = 3,
                    color = "white",
                    fill = "gray40") +
  theme_bw()

ggsave(
  plot,
  filename = file.path(r4projects::get_project_wd(), "4_manuscript/Figures/Figure_5/figure_5a.pdf"),
  width = 15,
  height = 7
)


# Translated comment.

library(networkD3)
library(dplyr)

# Translated comment.
# Translated comment.

# Translated comment.
# Translated comment.
oral_nodes <- unique(sankey_data$Oral)
gut_nodes <- unique(sankey_data$Gut)
metabolite_nodes <- unique(sankey_data$Metabolite)

# Translated comment.
nodes_df <- data.frame(
  name = c(oral_nodes, gut_nodes, metabolite_nodes),
  group = c(
    rep("Oral", length(oral_nodes)),
    rep("Gut", length(gut_nodes)),
    rep("Metabolite", length(metabolite_nodes))
  ),
  stringsAsFactors = FALSE
)

# Translated comment.
nodes_df$ID <- 0:(nrow(nodes_df) - 1)

# Translated comment.
# Translated comment.
links_oral_gut <- sankey_data %>%
  select(Oral, Gut) %>%
  group_by(Oral, Gut) %>%
  dplyr::summarise(value = n(), .groups = 'drop') %>%
  dplyr::rename("source" = "Oral", "target" = "Gut")

# Translated comment.
links_gut_metabolite <- sankey_data %>%
  select(Gut, Metabolite) %>%
  group_by(Gut, Metabolite) %>%
  dplyr::summarise(value = n(), .groups = 'drop') %>%
  dplyr::rename("source" = "Gut", "target" = "Metabolite")

# Translated comment.
# Translated comment.
links_oral_gut$group <- "Oral"
links_gut_metabolite$group <- "Gut"
links_df <- rbind(links_oral_gut, links_gut_metabolite)

# Translated comment.
links_df$source <- match(links_df$source, nodes_df$name) - 1
links_df$target <- match(links_df$target, nodes_df$name) - 1

# Translated comment.
colourScale <- JS(
  paste0(
    'd3.scaleOrdinal()
  .domain(["Oral", "Gut", "Metabolite"])
  .range(["#a1d5b9", "#Edd064", "#B6C7EA"])'
  )
)

# Translated comment.
sankey_plot <- sankeyNetwork(
  Links = links_df,
  Nodes = nodes_df,
  Source = "source",
  Target = "target",
  Value = "value",
  NodeID = "name",
  NodeGroup = "group",
  # Translated comment.
  LinkGroup = "group",
  # Translated comment.
  colourScale = colourScale,
  # Translated comment.
  fontSize = 12,
  nodeWidth = 30,
  nodePadding = 10,
  margin = list(left = 50, right = 50),
  sinksRight = TRUE,
  units = "个数"
)


sankey_plot
ggsave(
  sankey_plot,
  filename = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_5/figure_5a.pdf"
  ),
  width = 8,
  height = 5
)

htmlwidgets::saveWidget(
  sankey_plot,
  file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_5/figure_5a.html"
  ),
  selfcontained = TRUE
)

# Convert to PDF
pagedown::chrome_print(
  file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_5/figure_5a.html"
  ),
  output = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_5/figure_5a.pdf"
  )
)


# Translated comment.
colnames(bidirectional_mediation_results_sig_oral_gut_metabolite)[9:11] <-
  c("Oral", "Gut", "Metabolite")


calculate_correlations <- function(gut_data,
                                   oral_data,
                                   metabolome_data,
                                   mediation_results) {
  # Translated comment.
  correlation_results <- list()
  
  # Translated comment.
  for (i in 1:nrow(mediation_results)) {
    # Translated comment.
    gut_feature <- mediation_results$gut_feature[i]
    oral_feature <- mediation_results$oral_feature[i]
    metabolite_feature <- mediation_results$metabolite[i]
    
    # Translated comment.
    gut_values <- gut_data[gut_feature, ]
    oral_values <- oral_data[oral_feature, ]
    metabolite_values <- metabolome_data[metabolite_feature, ]
    
    # Translated comment.
    gut_values <- as.numeric(gut_values)
    oral_values <- as.numeric(oral_values)
    metabolite_values <- as.numeric(metabolite_values)
    
    # Translated comment.
    combined_data <- data.frame(gut = gut_values,
                                oral = oral_values,
                                metabolite = metabolite_values)
    
    # Translated comment.
    combined_data <- na.omit(combined_data)
    
    # Translated comment.
    if (nrow(combined_data) >= 3) {
      # Translated comment.
      cor_gut_oral <- cor.test(combined_data$gut, combined_data$oral, method = "spearman")
      cor_gut_metabolite <- cor.test(combined_data$gut, combined_data$metabolite, method = "spearman")
      cor_oral_metabolite <- cor.test(combined_data$oral, combined_data$metabolite, method = "spearman")
      
      # Translated comment.
      result <- data.frame(
        Mediation_Row = i,
        Gut_Feature = gut_feature,
        Oral_Feature = oral_feature,
        Metabolite_Feature = metabolite_feature,
        Gut_Oral_Cor = cor_gut_oral$estimate,
        Gut_Oral_Pvalue = cor_gut_oral$p.value,
        Gut_Metabolite_Cor = cor_gut_metabolite$estimate,
        Gut_Metabolite_Pvalue = cor_gut_metabolite$p.value,
        Oral_Metabolite_Cor = cor_oral_metabolite$estimate,
        Oral_Metabolite_Pvalue = cor_oral_metabolite$p.value,
        Sample_Size = nrow(combined_data)
      )
      
      # Translated comment.
      correlation_results[[i]] <- result
    } else {
      # Translated comment.
      result <- data.frame(
        Mediation_Row = i,
        Gut_Feature = gut_feature,
        Oral_Feature = oral_feature,
        Metabolite_Feature = metabolite_feature,
        Gut_Oral_Cor = NA,
        Gut_Oral_Pvalue = NA,
        Gut_Metabolite_Cor = NA,
        Gut_Metabolite_Pvalue = NA,
        Oral_Metabolite_Cor = NA,
        Oral_Metabolite_Pvalue = NA,
        Sample_Size = nrow(combined_data)
      )
      
      # Translated comment.
      correlation_results[[i]] <- result
    }
  }
  
  # Translated comment.
  final_results <- do.call(rbind, correlation_results)
  
  return(final_results)
}

# Translated comment.
correlation_results <- calculate_correlations(
  gut_data,
  oral_data,
  metabolome_data,
  bidirectional_mediation_results_sig_oral_gut_metabolite
)

bidirectional_mediation_results_sig_oral_gut_metabolite <- cbind(bidirectional_mediation_results_sig_oral_gut_metabolite,
                                                                 correlation_results[, 5:10])
