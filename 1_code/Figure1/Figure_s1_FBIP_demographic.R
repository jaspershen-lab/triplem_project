rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")
library(tidyverse)
library(tidymass)
library(readxl)



demographic_data<-read.table("2_data/FBIP-main/metadata/FBIP_metadata.txt")

# Filter0.
demographic_data<-subset(demographic_data,Time=="W0")


library(ComplexHeatmap)
library(circlize)
library(grid)


demographic_data$Gender<-gsub("F","Female",demographic_data$Gender)
demographic_data$Gender<-gsub("M","Male",demographic_data$Gender)
# Set0.7(30%).
alpha_value = 0.7

# Createannotationsfor .
ha = columnAnnotation(
  # variablesUsecolors,Add.
  Age = anno_barplot(demographic_data$Age, 
                     gp = gpar(fill = scales::alpha("#E69F00", alpha_value)),
                     border = TRUE,
                     width = unit(2, "cm")),  # Use.
  BMI = anno_barplot(demographic_data$BMI,
                     gp = gpar(fill = scales::alpha("#56B4E9", alpha_value)),
                     border = TRUE,
                     width = unit(2, "cm")),  # Use.
  TC = anno_barplot(demographic_data$TC,
                     gp = gpar(fill = scales::alpha("#F21451", alpha_value)),
                     border = TRUE,
                     width = unit(2, "cm")),  
  Insulin = anno_barplot(demographic_data$Insulin,
                    gp = gpar(fill = scales::alpha("#F45112", alpha_value)),
                    border = TRUE,
                    width = unit(2, "cm")),  
  
  border = TRUE,
  # variables.
  Gender = demographic_data$Gender,

  # Setvariablescolors,Add.
  col = list(
    Gender = mapply(function(x) scales::alpha(x, alpha_value), sex_color, USE.NAMES = TRUE)
  ),
  
  # Setannotations.
  annotation_name_gp = gpar(fontsize = 10),
  annotation_name_side = "left",
  simple_anno_size = unit(0.5, "cm"),
  
  # Addlegend.
  show_legend = TRUE,
  annotation_legend_param = list(
    Gender = list(title = "Gender")
  ),
  
  # Set.
  gap = unit(c(2, 2, 1, 1, 1), "mm")  # betweenAdd.
)

# CreatematrixPlotheatmap.
mat = matrix(0, nrow = 1, ncol = nrow(demographic_data))

# CreatePlotheatmap,.
ht = Heatmap(mat,
             top_annotation = ha,
             show_row_names = FALSE,
             show_column_names = FALSE,
             show_heatmap_legend = FALSE,
             cluster_rows = FALSE,
             cluster_columns = FALSE,
             height = unit(0.2, "cm"))   # Set.

# PlotSetlegend.



pdf("4_manuscript/Figures/Figure_1/figure_s1_FBIP_demographic.pdf", width = 10, height = 6)
ht_drawn <-draw(ht, annotation_legend_side = "bottom")  # Plotfor .
dev.off()