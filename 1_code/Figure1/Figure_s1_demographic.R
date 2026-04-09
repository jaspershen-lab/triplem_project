rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")
source("1_code/mantel_Procrustes_code.R")
library(tidyverse)
library(tidymass)
library(readxl)



load("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section")
metabolomics_object<-object_cross_section


demographic_data<-data.frame(metabolomics_object@sample_info)



# Translated comment.
demographic_data_Continuous<-demographic_data[,c("adjusted_age","BMI")]
rownames(demographic_data_Continuous)<-demographic_data$sample_id

# Translated comment.
demographic_data_Discrete<-demographic_data[,c("IRIS","Gender","Ethnicity")]


rownames(demographic_data_Discrete)<-demographic_data$sample_id



library(ComplexHeatmap)
library(circlize)
library(grid)


library(ComplexHeatmap)
library(circlize)
library(grid)

# Translated comment.
alpha_value = 0.7

# Translated comment.
ha = columnAnnotation(
  # Translated comment.
  Age = anno_barplot(demographic_data$adjusted_age, 
                     gp = gpar(fill = scales::alpha("#E69F00", alpha_value)),
                     border = TRUE,
                     width = unit(2, "cm")),  # translated comment
  BMI = anno_barplot(demographic_data$BMI,
                     gp = gpar(fill = scales::alpha("#56B4E9", alpha_value)),
                     border = TRUE,
                     width = unit(2, "cm")),  # translated comment
  
  border = TRUE,
  # Translated comment.
  IRIS = demographic_data$IRIS,
  Gender = demographic_data$Gender,
  Ethnicity = demographic_data$Ethnicity,
  
  # Translated comment.
  col = list(
    IRIS = mapply(function(x) scales::alpha(x, alpha_value), iris_color, USE.NAMES = TRUE),
    Gender = mapply(function(x) scales::alpha(x, alpha_value), sex_color, USE.NAMES = TRUE),
    Ethnicity = mapply(function(x) scales::alpha(x, alpha_value), ethnicity_color, USE.NAMES = TRUE)
  ),
  
  # Translated comment.
  annotation_name_gp = gpar(fontsize = 10),
  annotation_name_side = "left",
  simple_anno_size = unit(0.5, "cm"),
  
  # Translated comment.
  show_legend = TRUE,
  annotation_legend_param = list(
    IRIS = list(title = "IRIS"),
    Gender = list(title = "Gender"),
    Ethnicity = list(title = "Ethnicity")
  ),
  
  # Translated comment.
  gap = unit(c(2, 2, 1, 1, 1), "mm")  # translated comment
)

# Translated comment.
mat = matrix(0, nrow = 1, ncol = nrow(demographic_data))

# Translated comment.
ht = Heatmap(mat,
             top_annotation = ha,
             show_row_names = FALSE,
             show_column_names = FALSE,
             show_heatmap_legend = FALSE,
             cluster_rows = FALSE,    # translated comment
             cluster_columns = FALSE, # translated comment
             height = unit(0.2, "cm"))   # translated comment

# Translated comment.
draw(ht, annotation_legend_side = "bottom")



### Translated comment.



load("3_data_analysis/gut_microbiome/data_preparation/object_cross_section")

gut_object<-object_cross_section


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







shannon_div <- diversity(t(gut_object@expression_data), index = "shannon")

# Translated comment.
results_gut <- data.frame(
  Sample = names(shannon_div),
  Shannon = shannon_div
)


load("3_data_analysis/oral_microbiome/data_preparation/object_cross_section")

oral_object<-object_cross_section


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







shannon_div <- diversity(t(oral_object@expression_data), index = "shannon")

# Translated comment.
results_oral <- data.frame(
  Sample = names(shannon_div),
  Shannon = shannon_div
)


load("3_data_analysis/skin_microbiome/data_preparation/object_cross_section")

skin_object<-object_cross_section


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







shannon_div <- diversity(t(skin_object@expression_data), index = "shannon")

# Translated comment.
results_skin <- data.frame(
  Sample = names(shannon_div),
  Shannon = shannon_div
)


load("3_data_analysis/nasal_microbiome/data_preparation/object_cross_section")

nasal_object<-object_cross_section


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







shannon_div <- diversity(t(nasal_object@expression_data), index = "shannon")

# Translated comment.
results_nasal <- data.frame(
  Sample = names(shannon_div),
  Shannon = shannon_div
)


# Translated comment.
Sample_ID<-demographic_data[,1:2]
colnames(Sample_ID)[1]<-"Sample"

alpha_diversity<-Sample_ID%>%
  full_join(results_gut, by = "Sample") %>%
  full_join(results_oral, by = "Sample") %>%
  full_join(results_skin, by = "Sample") %>%
  full_join(results_nasal, by = "Sample")

rownames(alpha_diversity)<-alpha_diversity$Sample

alpha_diversity<-alpha_diversity[,-1:-2]

colnames(alpha_diversity)<-c("gut","oral","skin","nasal")


# Translated comment.
alpha_value = 0.7

# Translated comment.
microbiome_colors <- body_site_color


# Translated comment.
ha = columnAnnotation(
  # Translated comment.
  Age = anno_barplot(demographic_data$adjusted_age, 
                     gp = gpar(fill = scales::alpha("#E69F00", alpha_value)),
                     border = TRUE,
                     width = unit(2, "cm")),
  BMI = anno_barplot(demographic_data$BMI,
                     gp = gpar(fill = scales::alpha("#56B4E9", alpha_value)),
                     border = TRUE,
                     width = unit(2, "cm")),
  
  # Translated comment.
  `Gut` = anno_points(
    alpha_diversity$gut,
    gp = gpar(col = scales::alpha(microbiome_colors["gut"], alpha_value), 
              fill = scales::alpha(microbiome_colors["gut"], alpha_value)),
    pch = 16,  # translated comment
    size = unit(2, "mm"),  # translated comment
    border = TRUE,
    width = unit(2, "cm")
  ),
  
  `Oral` = anno_points(
    alpha_diversity$oral,
    gp = gpar(col = scales::alpha(microbiome_colors["oral"], alpha_value), 
              fill = scales::alpha(microbiome_colors["oral"], alpha_value)),
    pch = 16,
    size = unit(2, "mm"),
    border = TRUE,
    width = unit(2, "cm")
  ),
  
  `Skin` = anno_points(
    alpha_diversity$skin,
    gp = gpar(col = scales::alpha(microbiome_colors["skin"], alpha_value), 
              fill = scales::alpha(microbiome_colors["skin"], alpha_value)),
    pch = 16,
    size = unit(2, "mm"),
    border = TRUE,
    width = unit(2, "cm")
  ),
  
  `Nasal` = anno_points(
    alpha_diversity$nasal,
    gp = gpar(col = scales::alpha(microbiome_colors["nasal"], alpha_value), 
              fill = scales::alpha(microbiome_colors["nasal"], alpha_value)),
    pch = 16,
    size = unit(2, "mm"),
    border = TRUE,
    width = unit(2, "cm")
  ),
  
  border = TRUE,
  # Translated comment.
  IRIS = demographic_data$IRIS,
  Gender = demographic_data$Gender,
  Ethnicity = demographic_data$Ethnicity,
  
  # Translated comment.
  col = list(
    IRIS = mapply(function(x) scales::alpha(x, alpha_value), iris_color, USE.NAMES = TRUE),
    Gender = mapply(function(x) scales::alpha(x, alpha_value), sex_color, USE.NAMES = TRUE),
    Ethnicity = mapply(function(x) scales::alpha(x, alpha_value), ethnicity_color, USE.NAMES = TRUE)
  ),
  
  # Translated comment.
  annotation_name_gp = gpar(fontsize = 10),
  annotation_name_side = "left",
  simple_anno_size = unit(0.5, "cm"),
  
  # Translated comment.
  show_legend = TRUE,
  annotation_legend_param = list(
    IRIS = list(title = "IRIS"),
    Gender = list(title = "Gender"),
    Ethnicity = list(title = "Ethnicity")
  ),
  
  # Translated comment.
  gap = unit(c(2, 2, 2, 2, 2, 2, 1, 1, 1), "mm")  # translated comment
)

# Translated comment.
mat = matrix(0, nrow = 1, ncol = nrow(demographic_data))

# Translated comment.
ht = Heatmap(mat,
             top_annotation = ha,
             show_row_names = FALSE,
             show_column_names = FALSE,
             show_heatmap_legend = FALSE,
             cluster_rows = FALSE,    # translated comment
             cluster_columns = FALSE, # translated comment
             height = unit(0.2, "cm"))   # translated comment



# Translated comment.
pdf("4_manuscript/Figures/Figure_1/figure_s1_demographic.pdf", width = 12, height = 6)
draw(ht, annotation_legend_side = "bottom")
dev.off()