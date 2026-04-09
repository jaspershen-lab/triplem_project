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
common_samples <- Reduce(intersect, list(
  colnames(gut_data),
  colnames(oral_data),
  colnames(metabolome_data)
))

# Function: calculate Spearman correlations and apply FDR correction.
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
  
  # FDR correction.
  fdr_values <- matrix(p.adjust(p_values, method = "BH"), nrow = nrow(p_values))
  rownames(fdr_values) <- rownames(p_values)
  colnames(fdr_values) <- colnames(p_values)
  
  return(list(cor = result, p = p_values, fdr = fdr_values))
}
# AddIRISAnalyze.
# sample informationin IRIS.
rownames(gut_sample_info)<-gut_sample_info$sample_id
rownames(oral_sample_info)<-oral_sample_info$sample_id
# GetIRIS(gut_sample_infoin IRIS).
# : if IRIScolumn names"IRIS",column names.
ir_samples <- row.names(gut_sample_info)[gut_sample_info$IRIS == "IR"]
is_samples <- row.names(gut_sample_info)[gut_sample_info$IRIS == "IS"]

# samplesdatasetin .
ir_common_samples <- intersect(ir_samples, common_samples)
is_common_samples <- intersect(is_samples, common_samples)

# IRExtractdata.
gut_data_ir <- gut_data[, ir_common_samples]
oral_data_ir <- oral_data[, ir_common_samples]
metabolome_data_ir <- metabolome_data[, ir_common_samples]

# ISExtractdata.
gut_data_is <- gut_data[, is_common_samples]
oral_data_is <- oral_data[, is_common_samples]
metabolome_data_is <- metabolome_data[, is_common_samples]

# TransposeIRdataAnalyze.
gut_data_ir_t <- t(gut_data_ir)
oral_data_ir_t <- t(oral_data_ir)
metabolome_data_ir_t <- t(metabolome_data_ir)

# TransposeISdataAnalyze.
gut_data_is_t <- t(gut_data_is)
oral_data_is_t <- t(oral_data_is)
metabolome_data_is_t <- t(metabolome_data_is)

# IRcorrelation analysis.
gut_metabolome_cor_ir <- calculate_correlations(gut_data_ir_t, metabolome_data_ir_t)
oral_metabolome_cor_ir <- calculate_correlations(oral_data_ir_t, metabolome_data_ir_t)
oral_gut_cor_ir <- calculate_correlations(oral_data_ir_t, gut_data_ir_t)

# IScorrelation analysis.
gut_metabolome_cor_is <- calculate_correlations(gut_data_is_t, metabolome_data_is_t)
oral_metabolome_cor_is <- calculate_correlations(oral_data_is_t, metabolome_data_is_t)
oral_gut_cor_is <- calculate_correlations(oral_data_is_t, gut_data_is_t)

# IR.
gut_metabolome_associations_ir <- which(gut_metabolome_cor_ir$p < 0.05, arr.ind = TRUE)
oral_metabolome_associations_ir <- which(oral_metabolome_cor_ir$p < 0.05, arr.ind = TRUE)
oral_gut_associations_ir <- which(oral_gut_cor_ir$p < 0.05, arr.ind = TRUE)

# IS.
gut_metabolome_associations_is <- which(gut_metabolome_cor_is$p < 0.05, arr.ind = TRUE)
oral_metabolome_associations_is <- which(oral_metabolome_cor_is$p < 0.05, arr.ind = TRUE)
oral_gut_associations_is <- which(oral_gut_cor_is$p < 0.05, arr.ind = TRUE)

# CreateIRresultsdata frame - Setgroup.
bidirectional_mediation_results_ir <- data.frame(
  direction = character(),
  oral_feature = character(),
  gut_feature = character(),
  metabolite = character(),
  ACME = numeric(),
  ACME_p = numeric(),
  prop_mediated = numeric(),
  stringsAsFactors = FALSE
)

# CreateISresultsdata frame - Setgroup.
bidirectional_mediation_results_is <- data.frame(
  direction = character(),
  oral_feature = character(),
  gut_feature = character(),
  metabolite = character(),
  ACME = numeric(),
  ACME_p = numeric(),
  prop_mediated = numeric(),
  stringsAsFactors = FALSE
)


# AnalyzeafterAdd.
if(nrow(bidirectional_mediation_results_ir) > 0) {
  bidirectional_mediation_results_ir$group <- "IR"
}

if(nrow(bidirectional_mediation_results_is) > 0) {
  bidirectional_mediation_results_is$group <- "IS"
}
# IR: oralmicrobiome -> gutmicrobiome -> metabolites.
for (i in 1:nrow(oral_gut_associations_ir)) {
  oral_idx <- oral_gut_associations_ir[i, "row"]
  gut_idx <- oral_gut_associations_ir[i, "col"]
  
  oral_feature <- rownames(oral_gut_cor_ir$cor)[oral_idx]
  gut_feature <- colnames(oral_gut_cor_ir$cor)[gut_idx]
  
  # Checkgutmicrobiomewhether and metabolites.
  gut_metabolite_associations <- which(gut_metabolome_cor_ir$p[gut_idx, ] <  0.0001, arr.ind = TRUE)
  
  if (length(gut_metabolite_associations) > 0) {
    for (j in 1:length(gut_metabolite_associations)) {
      met_idx <- gut_metabolite_associations[j]
      metabolite <- colnames(gut_metabolome_cor_ir$cor)[met_idx]
      
      # data framein Analyze.
      med_data <- data.frame(oral = oral_data_ir_t[, oral_feature],
                             gut = gut_data_ir_t[, gut_feature],
                             metabolite = metabolome_data_ir_t[, metabolite])
      
      # Oralmicrobiome -> gutmicrobiome -> metabolites.
      med_model <- lm(gut ~ oral, data = med_data)
      out_model <- lm(metabolite ~ oral + gut + oral:gut, data = med_data)
      
      # in Analyze.
      library(mediation)
      med_result <- mediate(
        med_model,
        out_model,
        treat = "oral",
        mediator = "gut",
        boot = TRUE,
        sims = 10
      )
      
      # Saveresults.
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
      
      bidirectional_mediation_results_ir <- rbind(bidirectional_mediation_results_ir, result_row)
    }
  }
}

# IS: oralmicrobiome -> gutmicrobiome -> metabolites.
for (i in 1:nrow(oral_gut_associations_is)) {
  oral_idx <- oral_gut_associations_is[i, "row"]
  gut_idx <- oral_gut_associations_is[i, "col"]
  
  oral_feature <- rownames(oral_gut_cor_is$cor)[oral_idx]
  gut_feature <- colnames(oral_gut_cor_is$cor)[gut_idx]
  
  # Checkgutmicrobiomewhether and metabolites.
  gut_metabolite_associations <- which(gut_metabolome_cor_is$p[gut_idx, ] <  0.0001, arr.ind = TRUE)
  
  if (length(gut_metabolite_associations) > 0) {
    for (j in 1:length(gut_metabolite_associations)) {
      met_idx <- gut_metabolite_associations[j]
      metabolite <- colnames(gut_metabolome_cor_is$cor)[met_idx]
      
      # data framein Analyze.
      med_data <- data.frame(oral = oral_data_is_t[, oral_feature],
                             gut = gut_data_is_t[, gut_feature],
                             metabolite = metabolome_data_is_t[, metabolite])
      
      # Oralmicrobiome -> gutmicrobiome -> metabolites.
      med_model <- lm(gut ~ oral, data = med_data)
      out_model <- lm(metabolite ~ oral + gut + oral:gut, data = med_data)
      
      # in Analyze.
      med_result <- mediate(
        med_model,
        out_model,
        treat = "oral",
        mediator = "gut",
        boot = TRUE,
        sims = 10
      )
      
      # Saveresults.
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
      
      bidirectional_mediation_results_is <- rbind(bidirectional_mediation_results_is, result_row)
    }
  }
}

# IR: gutmicrobiome -> oralmicrobiome -> metabolites.
for (i in 1:nrow(oral_gut_associations_ir)) {
  oral_idx <- oral_gut_associations_ir[i, "row"]
  gut_idx <- oral_gut_associations_ir[i, "col"]
  
  oral_feature <- rownames(oral_gut_cor_ir$cor)[oral_idx]
  gut_feature <- colnames(oral_gut_cor_ir$cor)[gut_idx]
  
  # Checkoralmicrobiomewhether and metabolites.
  oral_metabolite_associations <- which(oral_metabolome_cor_ir$p[oral_idx, ] <  0.0001, arr.ind = TRUE)
  
  if (length(oral_metabolite_associations) > 0) {
    for (j in 1:length(oral_metabolite_associations)) {
      met_idx <- oral_metabolite_associations[j]
      metabolite <- colnames(oral_metabolome_cor_ir$cor)[met_idx]
      
      # data framein Analyze.
      med_data <- data.frame(gut = gut_data_ir_t[, gut_feature],
                             oral = oral_data_ir_t[, oral_feature],
                             metabolite = metabolome_data_ir_t[, metabolite])
      
      # Gutmicrobiome -> oralmicrobiome -> metabolites.
      med_model <- lm(oral ~ gut, data = med_data)
      out_model <- lm(metabolite ~ gut + oral + gut:oral, data = med_data)
      
      # in Analyze.
      med_result <- mediate(
        med_model,
        out_model,
        treat = "gut",
        mediator = "oral",
        boot = TRUE,
        sims = 10
      )
      
      # Saveresults.
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
      
      bidirectional_mediation_results_ir <- rbind(bidirectional_mediation_results_ir, result_row)
    }
  }
}

# IS: gutmicrobiome -> oralmicrobiome -> metabolites.
for (i in 1:nrow(oral_gut_associations_is)) {
  oral_idx <- oral_gut_associations_is[i, "row"]
  gut_idx <- oral_gut_associations_is[i, "col"]
  
  oral_feature <- rownames(oral_gut_cor_is$cor)[oral_idx]
  gut_feature <- colnames(oral_gut_cor_is$cor)[gut_idx]
  
  # Checkoralmicrobiomewhether and metabolites.
  oral_metabolite_associations <- which(oral_metabolome_cor_is$p[oral_idx, ] <  0.0001, arr.ind = TRUE)
  
  if (length(oral_metabolite_associations) > 0) {
    for (j in 1:length(oral_metabolite_associations)) {
      met_idx <- oral_metabolite_associations[j]
      metabolite <- colnames(oral_metabolome_cor_is$cor)[met_idx]
      
      # data framein Analyze.
      med_data <- data.frame(gut = gut_data_is_t[, gut_feature],
                             oral = oral_data_is_t[, oral_feature],
                             metabolite = metabolome_data_is_t[, metabolite])
      
      # Gutmicrobiome -> oralmicrobiome -> metabolites.
      med_model <- lm(oral ~ gut, data = med_data)
      out_model <- lm(metabolite ~ gut + oral + gut:oral, data = med_data)
      
      # in Analyze.
      med_result <- mediate(
        med_model,
        out_model,
        treat = "gut",
        mediator = "oral",
        boot = TRUE,
        sims = 10
      )
      
      # Saveresults.
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
      
      bidirectional_mediation_results_is <- rbind(bidirectional_mediation_results_is, result_row)
    }
  }
}

# FDR - IR.
if (nrow(bidirectional_mediation_results_ir) > 0) {
  bidirectional_mediation_results_ir$ACME_fdr <-
    p.adjust(bidirectional_mediation_results_ir$ACME_p, method = "BH")
}

# FDR - IS.
if (nrow(bidirectional_mediation_results_is) > 0) {
  bidirectional_mediation_results_is$ACME_fdr <-
    p.adjust(bidirectional_mediation_results_is$ACME_p, method = "BH")
}

# Filtersignificantresults - IR.
bidirectional_mediation_results_ir_sig <- subset(bidirectional_mediation_results_ir, ACME_p < 0.05)

# Filtersignificantresults - IS.
bidirectional_mediation_results_is_sig <- subset(bidirectional_mediation_results_is, ACME_p < 0.05)

# Addannotations - IR.
oral_tax <- oral_temp_object@variable_info
gut_tax <- gut_temp_object@variable_info

bidirectional_mediation_results_ir_sig <- merge(bidirectional_mediation_results_ir_sig,
                                                oral_tax[, c("Genus", "variable_id")],
                                                by.x = "oral_feature",
                                                by.y = "variable_id")

bidirectional_mediation_results_ir_sig <- merge(bidirectional_mediation_results_ir_sig,
                                                gut_tax[, c("Genus", "variable_id")],
                                                by.x = "gut_feature",
                                                by.y = "variable_id")

metabolite_annotation <- read_excel(
  "3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx"
)

bidirectional_mediation_results_ir_sig <- merge(
  bidirectional_mediation_results_ir_sig,
  metabolite_annotation[, c("HMDB.Name", "variable_id")],
  by.x = "metabolite",
  by.y = "variable_id"
)

# Addannotations - IS.
bidirectional_mediation_results_is_sig <- merge(bidirectional_mediation_results_is_sig,
                                                oral_tax[, c("Genus", "variable_id")],
                                                by.x = "oral_feature",
                                                by.y = "variable_id")

bidirectional_mediation_results_is_sig <- merge(bidirectional_mediation_results_is_sig,
                                                gut_tax[, c("Genus", "variable_id")],
                                                by.x = "gut_feature",
                                                by.y = "variable_id")

bidirectional_mediation_results_is_sig <- merge(
  bidirectional_mediation_results_is_sig,
  metabolite_annotation[, c("HMDB.Name", "variable_id")],
  by.x = "metabolite",
  by.y = "variable_id"
)

# Summarize.
ir_direction_counts <- table(bidirectional_mediation_results_ir_sig$direction)
is_direction_counts <- table(bidirectional_mediation_results_is_sig$direction)

# Createdata framefor .
ir_counts_df <- data.frame(
  Direction = names(ir_direction_counts),
  Count = as.numeric(ir_direction_counts),
  Group = "IR"
)

is_counts_df <- data.frame(
  Direction = names(is_direction_counts),
  Count = as.numeric(is_direction_counts),
  Group = "IS"
)

# Mergedata.
comparison_df <- rbind(ir_counts_df, is_counts_df)

# Createbar chart.
plot <-
comparison_df %>% 
  dplyr::filter(Direction == "oral->gut->metabolite") %>% 
  dplyr::mutate(Group = factor(Group, levels = c("IS", "IR"))) %>%
  ggplot(aes(x = Group, y = Count, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge", show.legend = FALSE) +
  geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.5) +
  theme_bw() +
  scale_fill_manual(values = c("IS" = "#1f77b4", "IR" = "#ff7f0e")) +
  labs(x = "", y = "Count", title = "Number of significant mediation paths") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave(plot, 
       filename = "4_manuscript/Figures/Figure_6/number of significant mediation paths.pdf",
       width = 4, height = 5)

# Extractoral->gut->metabolites.
ir_og_path <- subset(bidirectional_mediation_results_ir_sig, direction == "oral->gut->metabolite")
is_og_path <- subset(bidirectional_mediation_results_is_sig, direction == "oral->gut->metabolite")

# Createoral->gut->metabolitesSankey diagram - IR.
if(nrow(ir_og_path) > 0) {
  ir_sankey_data <- ir_og_path[, c("Genus.x", "Genus.y", "HMDB.Name")]
  colnames(ir_sankey_data) <- c("Oral", "Gut", "Metabolite")
  
  library(ggsankey)
  
  ir_df2 <- ir_sankey_data %>%
    make_long(Oral, Gut, Metabolite)
  
  ir_plot <- ggplot(
    ir_df2,
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
    theme_bw() +
    theme(legend.position = "none")+
    labs(title = "IR Group: Oral -> Gut -> Metabolite Pathway")
}

# Createoral->gut->metabolitesSankey diagram - IS.
if(nrow(is_og_path) > 0) {
  is_sankey_data <- is_og_path[, c("Genus.x", "Genus.y", "HMDB.Name")]
  colnames(is_sankey_data) <- c("Oral", "Gut", "Metabolite")
  
  library(ggsankey)
  
  is_df2 <- is_sankey_data %>%
    make_long(Oral, Gut, Metabolite)
  
  is_plot <- ggplot(
    is_df2,
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
    theme_bw() +
    theme(legend.position = "none")+
    labs(title = "IS Group: Oral -> Gut -> Metabolite Pathway")

}
is_plot
ir_plot

