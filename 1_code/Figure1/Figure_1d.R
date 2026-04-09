rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")
setwd("1_code/4_site_merge/")
library(tidymass)
library(tidyverse)
library(readxl)
###load("data)
load("../../3_data_analysis/gut_microbiome/data_preparation/object_cross_section")

gut_object <- object_cross_section

load(
  "../../3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section"
)

metabolomics_object <- object_cross_section


load("../../3_data_analysis/oral_microbiome/data_preparation/object_cross_section")

oral_object <- object_cross_section

load("../../3_data_analysis/skin_microbiome/data_preparation/object_cross_section")

skin_object <- object_cross_section

load("../../3_data_analysis/nasal_microbiome/data_preparation/object_cross_section")
nasal_object <- object_cross_section
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










# Translated comment.

# Translated comment.
# Translated comment.
gut_genus <- gut_object@expression_data
rownames(gut_genus) <- gut_object@variable_info$Genus

oral_genus <- oral_object@expression_data
rownames(oral_genus) <- oral_object@variable_info$Genus

skin_genus <- skin_object@expression_data
rownames(skin_genus) <- skin_object@variable_info$Genus

nasal_genus <- nasal_object@expression_data
rownames(nasal_genus) <- nasal_object@variable_info$Genus

# Translated comment.
library(vegan)      # translated comment
library(ggplot2)    # translated comment
library(readr)      # translated comment
library(dplyr)      # translated comment
library(tidyr)      # translated comment


# Translated comment.
gut_samples <- colnames(gut_genus)
oral_samples <- colnames(oral_genus)
skin_samples <- colnames(skin_genus)
nasal_samples <- colnames(nasal_genus)

# Translated comment.
# Translated comment.
gut_t <- t(gut_genus)
oral_t <- t(oral_genus)
skin_t <- t(skin_genus)
nasal_t <- t(nasal_genus)

# Translated comment.
# Translated comment.
rownames(gut_t) <- paste0(rownames(gut_t), "_gut")
rownames(oral_t) <- paste0(rownames(oral_t), "_oral")
rownames(skin_t) <- paste0(rownames(skin_t), "_skin")
rownames(nasal_t) <- paste0(rownames(nasal_t), "_nasal")

# Translated comment.
gut_labels <- data.frame(
  Sample = rownames(gut_t),
  Site = "Gut",
  Subject = sub("_gut$", "", rownames(gut_t))
)
oral_labels <- data.frame(
  Sample = rownames(oral_t),
  Site = "Oral",
  Subject = sub("_oral$", "", rownames(oral_t))
)
skin_labels <- data.frame(
  Sample = rownames(skin_t),
  Site = "Skin",
  Subject = sub("_skin$", "", rownames(skin_t))
)
nasal_labels <- data.frame(
  Sample = rownames(nasal_t),
  Site = "Nasal",
  Subject = sub("_nasal$", "", rownames(nasal_t))
)

# Translated comment.
# Translated comment.
all_species <- unique(c(
  colnames(gut_t),
  colnames(oral_t),
  colnames(skin_t),
  colnames(nasal_t)
))

# Translated comment.
fill_missing_species <- function(df, all_species) {
  # Translated comment.
  result <- matrix(0, nrow = nrow(df), ncol = length(all_species))
  rownames(result) <- rownames(df)
  colnames(result) <- all_species
  
  # Translated comment.
  common_species <- intersect(colnames(df), all_species)
  for (sp in common_species) {
    result[, sp] <- df[, sp]
  }
  
  # Translated comment.
  return(as.data.frame(result))
}

# Translated comment.
gut_complete <- fill_missing_species(gut_t, all_species)
oral_complete <- fill_missing_species(oral_t, all_species)
skin_complete <- fill_missing_species(skin_t, all_species)
nasal_complete <- fill_missing_species(nasal_t, all_species)

# Translated comment.
all_data <- rbind(gut_complete, oral_complete, skin_complete, nasal_complete)

# Translated comment.
sample_metadata <- rbind(gut_labels, oral_labels, skin_labels, nasal_labels)
rownames(sample_metadata) <- sample_metadata$Sample

# Translated comment.
sample_metadata <- sample_metadata[rownames(all_data), ]

# Translated comment.
library(microbiome)
library(phyloseq)
library(dplyr)
library(ggplot2)
library(gridExtra) # translated comment



# Translated comment.
otu_table <- otu_table(as.matrix(t(all_data)), taxa_are_rows = TRUE)

# Translated comment.
sample_data <- sample_data(sample_metadata)

# Translated comment.
physeq <- phyloseq(otu_table, sample_data)

# Translated comment.
site_list <- c("Gut", "Oral", "Skin", "Nasal")
core_taxa_results <- list()

# Translated comment.

# Translated comment.
for (site in site_list) {
  # Translated comment.
  site_samples <- subset_samples(physeq, Site == site)
  
  # Translated comment.
  # Translated comment.
  # Translated comment.
  core_detection <- 0.001  # translated comment
  core_prevalence <- 0.5   # translated comment
  
  # Translated comment.
  core_taxa <- core_members(site_samples, detection = core_detection, prevalence = core_prevalence)
  
  # Translated comment.
  if (length(core_taxa) > 10) {
    # Translated comment.
    taxa_sums <- taxa_sums(site_samples)
    taxa_sums <- taxa_sums[names(taxa_sums) %in% core_taxa]
    taxa_sums <- sort(taxa_sums, decreasing = TRUE)
    core_taxa <- names(taxa_sums)[1:10]
  }
  
  # Translated comment.
  core_taxa_results[[site]] <- core_taxa
  
  # Translated comment.
  cat("\n", site, "核心物种 (", length(core_taxa), "):\n", sep = "")
  print(core_taxa)
}

# Translated comment.
all_core_taxa <- unique(unlist(core_taxa_results))
cat("\n总共识别出", length(all_core_taxa), "个核心物种\n")

# Translated comment.
library(RColorBrewer)
# Translated comment.
plot_site_core <-
  function(physeq, site, core_taxa) {
    # Translated comment.
    site_samples <- subset_samples(physeq, Site == site)
    
    # Translated comment.
    site_samples_rel <- site_samples
    
    # Translated comment.
    site_samples_rel <- prune_taxa(all_core_taxa, site_samples_rel)
    
    
    prevalences <- seq(.05, 1, .02)
    
    detections <- round(10^seq(log10(0.0001), log10(.2), length = 30), 3)
    
    # Also define gray color palette
    gray <- rev(brewer.pal(5, "RdBu"))
    
    # Translated comment.
    # Translated comment.
    p <- plot_core(
      site_samples_rel,
      plot.type = "heatmap",
      colours = gray,
      prevalences = prevalences,
      taxa.order = rev(all_core_taxa),
      detections = detections
    ) +
      labs(x = "Detection threshold (Relative abundance)") +
      
      #Adjusts axis text size and legend bar height
      theme(
        axis.text.y = element_text(size = 14, face = "italic"),
        axis.text.x.bottom = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10)
      )
    
    return(p)
  }

# Translated comment.
plot_list <- list()
for (site in site_list) {
  plot_list[[site]] <- plot_site_core(physeq, site, core_taxa_results[[site]])
}

# Translated comment.
# combined_plot <- grid.arrange(
#   plot_list[["Gut"]],
#   plot_list[["Oral"]],
#   plot_list[["Skin"]],
#   plot_list[["Nasal"]],
#   nrow = 1
# )

p_gut <- plot_list[["Gut"]] +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )
p_oral <- plot_list[["Oral"]] + theme(
  legend.position = "none",
  axis.text.y  = element_blank(),
  axis.tricks.y = element_blank(),
  panel.grid = element_blank(),
  panel.background = element_blank(),
  plot.background = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
  panel.border = element_rect(colour = "black", fill = NA, size = 1)
)
p_skin <- plot_list[["Skin"]] + theme(
  legend.position = "none",
  axis.text.y  = element_blank(),
  axis.tricks.y = element_blank(),
  panel.grid = element_blank(),
  panel.background = element_blank(),
  plot.background = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
  panel.border = element_rect(colour = "black", fill = NA, size = 1)
)
p_nasal <- plot_list[["Nasal"]] + theme(
  legend.position = "none",
  axis.text.y  = element_blank(),
  axis.tricks.y = element_blank(),
  panel.grid = element_blank(),
  panel.background = element_blank(),
  plot.background = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
  panel.border = element_rect(colour = "black", fill = NA, size = 1)
)


# combined_plot <- grid.arrange(
#   p_gut,
#   p_oral,
#   p_skin,
#   p_nasal,
#   nrow = 1
# )

library(patchwork)

combined_plot <- p_gut + p_oral + p_skin + p_nasal + plot_layout(ncol = 4)
combined_plot

ggsave(
  "../../4_manuscript/Figures/Figure_1/figure_1d.pdf",
  plot = combined_plot,
  width = 13,
  height = 6
)
