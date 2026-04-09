rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")
setwd("1_code/4_site_merge/")
library(tidymass)
library(tidyverse)
library(readxl)
###load("data)
load("../../3_data_analysis/gut_microbiome/data_preparation/object_cross_section")

gut_object<-object_cross_section

load("../../3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section")

metabolomics_object<-object_cross_section


load("../../3_data_analysis/oral_microbiome/data_preparation/object_cross_section")

oral_object<-object_cross_section

load("../../3_data_analysis/skin_microbiome/data_preparation/object_cross_section")

skin_object<-object_cross_section

load("../../3_data_analysis/nasal_microbiome/data_preparation/object_cross_section")
nasal_object<-object_cross_section
metabolite_annotation<-read_excel("../../1_code/revise_code/metabolite_annotation.xlsx")
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








# Combine microbiome data from four body sites for a PCoA plot.

# Read microbiome data from four body regions.
# Assume the files are in the current working directory.
gut_genus<-gut_object@expression_data
rownames(gut_genus)<-gut_object@variable_info$Genus

oral_genus<-oral_object@expression_data
rownames(oral_genus)<-oral_object@variable_info$Genus

skin_genus<-skin_object@expression_data
rownames(skin_genus)<-skin_object@variable_info$Genus

nasal_genus<-nasal_object@expression_data
rownames(nasal_genus)<-nasal_object@variable_info$Genus

# Load required R packages.
library(vegan)      # CalculateNMDS.
library(ggplot2)    # Plot.
library(readr)      # ReadCSV.
library(dplyr)      # dataProcess.
library(tidyr)      # data.


# Add a source label for each dataset.
gut_samples <- colnames(gut_genus)
oral_samples <- colnames(oral_genus)
skin_samples <- colnames(skin_genus)
nasal_samples <- colnames(nasal_genus)

# Combine all data.
# Transpose the matrix so rows are samples and columns are taxa.
gut_t <- t(gut_genus)
oral_t <- t(oral_genus)
skin_t <- t(skin_genus)
nasal_t <- t(nasal_genus)

# Rename samples to avoid duplicates while keeping the original information.
# Assume identical sample names represent different body sites from the same person.
rownames(gut_t) <- paste0(rownames(gut_t), "_gut")
rownames(oral_t) <- paste0(rownames(oral_t), "_oral")
rownames(skin_t) <- paste0(rownames(skin_t), "_skin")
rownames(nasal_t) <- paste0(rownames(nasal_t), "_nasal")

# Create sample type labels.
gut_labels <- data.frame(Sample = rownames(gut_t), Site = "Gut", 
                         Subject = sub("_gut$", "", rownames(gut_t)))
oral_labels <- data.frame(Sample = rownames(oral_t), Site = "Oral", 
                          Subject = sub("_oral$", "", rownames(oral_t)))
skin_labels <- data.frame(Sample = rownames(skin_t), Site = "Skin", 
                          Subject = sub("_skin$", "", rownames(skin_t)))
nasal_labels <- data.frame(Sample = rownames(nasal_t), Site = "Nasal", 
                           Subject = sub("_nasal$", "", rownames(nasal_t)))

# Merge all taxa.
# First ensure all tables share the same taxa columns.
all_species <- unique(c(colnames(gut_t), colnames(oral_t), colnames(skin_t), colnames(nasal_t)))

# Update the helper that fills missing taxa to avoid index errors.
fill_missing_species <- function(df, all_species) {
  # Create a new data frame containing all possible taxa.
  result <- matrix(0, nrow = nrow(df), ncol = length(all_species))
  rownames(result) <- rownames(df)
  colnames(result) <- all_species
  
  # Fill in the existing data.
  common_species <- intersect(colnames(df), all_species)
  for (sp in common_species) {
    result[, sp] <- df[, sp]
  }
  
  # Convert to a data frame and return it.
  return(as.data.frame(result))
}

# Apply the updated function.
gut_complete <- fill_missing_species(gut_t, all_species)
oral_complete <- fill_missing_species(oral_t, all_species)
skin_complete <- fill_missing_species(skin_t, all_species)
nasal_complete <- fill_missing_species(nasal_t, all_species)

# Merge all sample data.
all_data <- rbind(gut_complete, oral_complete, skin_complete, nasal_complete)

# Merge sample labels.
sample_metadata <- rbind(gut_labels, oral_labels, skin_labels, nasal_labels)
rownames(sample_metadata) <- sample_metadata$Sample

# Ensure sample order matches.
sample_metadata <- sample_metadata[rownames(all_data), ]

# Calculate Bray-Curtis distances.
bray_dist <- vegdist(all_data, method = "bray")

# Run NMDS analysis.
set.seed(123)  # Set a random seed so results are reproducible.
nmds_result <- metaMDS(bray_dist, k = 2, trymax = 100, autotransform = FALSE)

# Check whether NMDS converged and report the stress value.
cat("NMDS Stress:", nmds_result$stress, "\n")
if(nmds_result$stress > 0.2) {
  warning("NMDS stress > 0.2, indicating poor ordination quality")
} else if(nmds_result$stress > 0.1) {
  cat("NMDS stress is between 0.1 and 0.2, indicating moderate ordination quality\n")
} else {
  cat("NMDS stress < 0.1, indicating good ordination quality\n")
}

# Extract NMDS coordinates.
nmds_df <- as.data.frame(nmds_result$points)
colnames(nmds_df) <- c("NMDS1", "NMDS2")

# Add sample information to the NMDS data.
nmds_df$Sample <- rownames(nmds_df)
nmds_df <- merge(nmds_df, sample_metadata, by = "Sample")

# Plot the NMDS results.
nmds_plot <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Site, shape = Site)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = Site), level = 0.95, linetype = 2) +
  scale_color_manual(values = c("Gut" = "#FF5733", "Oral" = "#33A8FF", 
                                "Skin" = "#33FF57", "Nasal" = "#D433FF")) +
  labs(title = "NMDS of Microbiome Samples Based on Bray-Curtis Distance",
       subtitle = paste("Stress =", round(nmds_result$stress, 3)),
       x = "NMDS1",
       y = "NMDS2") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

# Add a source label for each dataset.
gut_samples <- colnames(gut_genus)
oral_samples <- colnames(oral_genus)
skin_samples <- colnames(skin_genus)
nasal_samples <- colnames(nasal_genus)

# Combine all data.
# Transpose the matrix so rows are samples and columns are taxa.
gut_t <- t(gut_genus)
oral_t <- t(oral_genus)
skin_t <- t(skin_genus)
nasal_t <- t(nasal_genus)

# Rename samples to avoid duplicates while keeping the original information.
# Assume identical sample names represent different body sites from the same person.
rownames(gut_t) <- paste0(rownames(gut_t), "_gut")
rownames(oral_t) <- paste0(rownames(oral_t), "_oral")
rownames(skin_t) <- paste0(rownames(skin_t), "_skin")
rownames(nasal_t) <- paste0(rownames(nasal_t), "_nasal")

# Create sample type labels.
gut_labels <- data.frame(Sample = rownames(gut_t), Site = "Gut", 
                         Subject = sub("_gut$", "", rownames(gut_t)))
oral_labels <- data.frame(Sample = rownames(oral_t), Site = "Oral", 
                          Subject = sub("_oral$", "", rownames(oral_t)))
skin_labels <- data.frame(Sample = rownames(skin_t), Site = "Skin", 
                          Subject = sub("_skin$", "", rownames(skin_t)))
nasal_labels <- data.frame(Sample = rownames(nasal_t), Site = "Nasal", 
                           Subject = sub("_nasal$", "", rownames(nasal_t)))

# Merge all taxa.
# First ensure all tables share the same taxa columns.
all_species <- unique(c(colnames(gut_t), colnames(oral_t), colnames(skin_t), colnames(nasal_t)))

# Update the helper that fills missing taxa to avoid index errors.
fill_missing_species <- function(df, all_species) {
  # Create a new data frame containing all possible taxa.
  result <- matrix(0, nrow = nrow(df), ncol = length(all_species))
  rownames(result) <- rownames(df)
  colnames(result) <- all_species
  
  # Fill in the existing data.
  common_species <- intersect(colnames(df), all_species)
  for (sp in common_species) {
    result[, sp] <- df[, sp]
  }
  
  # Convert to a data frame and return it.
  return(as.data.frame(result))
}

# Apply the updated function.
gut_complete <- fill_missing_species(gut_t, all_species)
oral_complete <- fill_missing_species(oral_t, all_species)
skin_complete <- fill_missing_species(skin_t, all_species)
nasal_complete <- fill_missing_species(nasal_t, all_species)

# Merge all sample data.
all_data <- rbind(gut_complete, oral_complete, skin_complete, nasal_complete)

# Merge sample labels.
sample_metadata <- rbind(gut_labels, oral_labels, skin_labels, nasal_labels)
rownames(sample_metadata) <- sample_metadata$Sample

# Ensure sample order matches.
sample_metadata <- sample_metadata[rownames(all_data), ]

# Calculate Bray-Curtis distances.
bray_dist <- vegdist(all_data, method = "bray")

# Run NMDS analysis.
set.seed(123)  # Set a random seed so results are reproducible.
nmds_result <- metaMDS(bray_dist, k = 2, trymax = 100, autotransform = FALSE)

# Check whether NMDS converged and report the stress value.
cat("NMDS Stress:", nmds_result$stress, "\n")
if(nmds_result$stress > 0.2) {
  warning("NMDS stress > 0.2, indicating poor ordination quality")
} else if(nmds_result$stress > 0.1) {
  cat("NMDS stress is between 0.1 and 0.2, indicating moderate ordination quality\n")
} else {
  cat("NMDS stress < 0.1, indicating good ordination quality\n")
}

# Extract NMDS coordinates.
nmds_df <- as.data.frame(nmds_result$points)
colnames(nmds_df) <- c("NMDS1", "NMDS2")

# Add sample information to the NMDS data.
nmds_df$Sample <- rownames(nmds_df)
nmds_df <- merge(nmds_df, sample_metadata, by = "Sample")

# Plot the NMDS results.
nmds_plot<-ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, fill = Site,color=Site)) +
  geom_point(size = 4, alpha = 0.8,colour = "white",shape = 21) +
  stat_ellipse(aes(group = Site), level = 0.95, linetype = 2) +
  scale_fill_manual(values = c("Gut" = "#edd064", "Oral" = "#a1d5b9", 
                               "Skin" = "#f2ccac", "Nasal" = "#a17db4")) +
  scale_color_manual(values = c("Gut" = "#edd064", "Oral" = "#a1d5b9", 
                                "Skin" = "#f2ccac", "Nasal" = "#a17db4"))+
  labs(title = "NMDS of Microbiome Samples",
       subtitle = paste("Stress =", round(nmds_result$stress, 3)),
       x = "NMDS1",
       y = "NMDS2") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 12),
    panel.grid.minor = element_blank()
  )
nmds_plot
ggsave("../../4_manuscript/Figures/Figure_1/figure_1c.pdf", plot=nmds_plot, width=7, height=6, units="in")

## Summarize detected metabolite categories.


# Load required packages.
library(ggplot2)
library(dplyr)

class_counts <- as.data.frame(table(metabolite_annotation$HMDB.Class,metabolite_annotation$HMDB.Source.Microbial))
colnames(class_counts) <- c("HMDB.Class","HMDB.Source.Microbial", "Count")


class_counts<-subset(class_counts,!(HMDB.Class=="NA"))
class_counts<-subset(class_counts,!(HMDB.Source.Microbial=="NA"))
class_counts<-subset(class_counts,Count>=3)

total_by_class <- aggregate(Count ~ HMDB.Class, data = class_counts, sum)
# Sort by total count in ascending order.
class_order <- total_by_class$HMDB.Class[order(total_by_class$Count)]
# Convert HMDB.Class to an ordered factor.
class_counts$HMDB.Class <- factor(class_counts$HMDB.Class, levels = class_order)

# Then generate the plot.
p<-ggplot(class_counts, aes(x = HMDB.Class, y = Count, fill=HMDB.Source.Microbial)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#3d95d2", "#f16147"))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  coord_flip() +
  labs(x = "Metabolite class", y = "Metabolite number") +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12)
  )
p
ggsave("../../4_manuscript/Figures/Figure_1/figure_1b.pdf", plot=p, width=7, height=7, units="in")
### Source Venn diagram.

metabolite_annotation_Source<-metabolite_annotation[,29:31]

metabolite_annotation_Source<-subset(metabolite_annotation_Source,!(HMDB.Source.Endogenous=="NA"))
metabolite_annotation_Source<-subset(metabolite_annotation_Source,!(HMDB.Source.Food=="NA"))
metabolite_annotation_Source<-subset(metabolite_annotation_Source,!(HMDB.Source.Microbial=="NA"))


library(ggVennDiagram)
library(ggplot2)

# Assuming your data is already loaded as metabolite_annotation
# If not, use the following code to read it (replace with your file path)
# metabolite_annotation <- read.csv("path_to_your_file.csv", stringsAsFactors = FALSE)

# Create lists of metabolites for each source
endogenous_metabolites <- rownames(metabolite_annotation)[metabolite_annotation$HMDB.Source.Endogenous == TRUE]
food_metabolites <- rownames(metabolite_annotation)[metabolite_annotation$HMDB.Source.Food == TRUE]
microbial_metabolites <- rownames(metabolite_annotation)[metabolite_annotation$HMDB.Source.Microbial == TRUE]

# Create a list for the Venn diagram
venn_list <- list(
  Endogenous = endogenous_metabolites,
  Food = food_metabolites,
  Microbial = microbial_metabolites
)

# Create the Venn diagram
venn_plot <- ggVennDiagram(venn_list, 
                           label = "count",
                           category.names = c("Endogenous", "Food", "Microbial"))

# Customize the plot
venn_plot <- venn_plot + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none") +
  labs(title = "Distribution of Metabolite Sources",
       caption = "Data source: metabolite_annotation")

venn_plot

# Display the plot
ggsave("../../4_manuscript/Figures/Figure_1/figure_s1_venn.pdf", plot=venn_plot, width=5, height=5, units="in")
