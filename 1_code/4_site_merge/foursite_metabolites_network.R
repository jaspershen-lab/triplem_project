rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")
library(readxl)

metabolite_annotation<-read_excel("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx")


setwd("1_code/4_site_merge/")

gut_GBDT_results<-readRDS("../../3_data_analysis/gut_microbiome/GBDT/cross_section/gut_GBDT_results")

gut_GBDT_results_R2<-gut_GBDT_results$summary[,c(1,2)]
colnames(gut_GBDT_results_R2)<-c("metabolite","gut")

oral_GBDT_results<-readRDS("../../3_data_analysis/oral_microbiome/GBDT/cross_section/oral_GBDT_results")
oral_GBDT_results_R2<-oral_GBDT_results$summary[,c(1,2)]
colnames(oral_GBDT_results_R2)<-c("metabolite","oral")


skin_GBDT_results<-readRDS("../../3_data_analysis/skin_microbiome/GBDT/cross_section/skin_GBDT_results")
skin_GBDT_results_R2<-skin_GBDT_results$summary[,c(1,2)]
colnames(skin_GBDT_results_R2)<-c("metabolite","skin")

nasal_GBDT_results<-readRDS("../../3_data_analysis/nasal_microbiome/GBDT/cross_section/nasal_GBDT_results")
nasal_GBDT_results_R2<-nasal_GBDT_results$summary[,c(1,2)]
colnames(nasal_GBDT_results_R2)<-c("metabolite","nasal")



four_site_R2<-cbind(gut_GBDT_results_R2,oral_GBDT_results_R2[,2],skin_GBDT_results_R2[,2],nasal_GBDT_results_R2[,2])

rownames(four_site_R2)<-four_site_R2$metabolite
four_site_R2<-four_site_R2[,-1]
colnames(four_site_R2)<-c("gut","oral","skin","nasal")



four_site_R2<-ifelse(four_site_R2>0.05,"1","0")





library(tidygraph)
library(ggraph)
library(dplyr)
library(tidyr)

# Translated comment.
metabolite_to_network <- function(metabolite_matrix, annotation_df) {
  # Translated comment.
  if (!is.matrix(metabolite_matrix) && !is.data.frame(metabolite_matrix)) {
    stop("输入必须是矩阵或数据框")
  }
  
  # Translated comment.
  metabolite_df <- as.data.frame(metabolite_matrix)
  
  # Translated comment.
  metabolite_df$metabolite <- rownames(metabolite_matrix)
  
  # Translated comment.
  edges_df <- metabolite_df %>%
    pivot_longer(
      cols = c("gut", "oral", "skin", "nasal"),
      names_to = "body_site",
      values_to = "value"
    ) %>%
    mutate(value = as.numeric(value)) %>%  # translated comment
    filter(value > 0)  # translated comment
  
  # Translated comment.
  metabolite_nodes <- data.frame(
    name = unique(edges_df$metabolite),
    type = "metabolite"
  ) %>%
    left_join(
      annotation_df %>% 
        select(name = metabolite, HMDB.Class, HMDB.Name),
      by = "name"
    )
  
  body_site_nodes <- data.frame(
    name = c("gut", "oral", "skin", "nasal"),
    type = "body_site",
    HMDB.Class = NA,
    HMDB.Name = NA
  )
  
  # Translated comment.
  nodes <- bind_rows(metabolite_nodes, body_site_nodes)
  
  # Translated comment.
  edges <- edges_df %>%
    select(from = metabolite, to = body_site, weight = value)
  
  # Translated comment.
  graph <- tbl_graph(
    nodes = nodes,
    edges = edges,
    directed = FALSE
  )
  
  return(graph)
}

# Translated comment.
plot_metabolite_network <- function(graph) {
  # Translated comment.
  node_data <- as_tibble(graph)
  
  # Translated comment.
  top_classes <- node_data %>%
    filter(!is.na(HMDB.Class)) %>%
    count(HMDB.Class) %>%
    arrange(desc(n)) %>%
    slice_head(n = 8) %>%
    pull(HMDB.Class)
  
  # Translated comment.
  V(graph)$HMDB.Class_grouped <- ifelse(
    V(graph)$type == "body_site",
    "body_site",
    ifelse(V(graph)$HMDB.Class %in% top_classes,
           V(graph)$HMDB.Class,
           "Others")
  )
  
  # Translated comment.
  class_colors <- setNames(
    colorRampPalette(c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999"))(8),
    top_classes
  )
  # Translated comment.
  class_colors <- c(class_colors,
                    "Others" = "#CCCCCC",
                    "body_site" = "red")
  
  ggraph(graph, layout = "fr") +
    # Translated comment.
    geom_edge_link(aes(width = weight, alpha = weight),
                   color = "grey50") +
    scale_edge_width_continuous(range = c(0.2, 2)) +
    scale_edge_alpha_continuous(range = c(0.2, 0.8)) +
    # Translated comment.
    geom_node_point(aes(color = HMDB.Class_grouped, 
                        size = type)) +
    scale_color_manual(values = class_colors,
                       name = "Metabolite Class") +
    scale_size_manual(values = c(
      metabolite = 3,
      body_site = 8
    )) +
    # Translated comment.
    geom_node_text(aes(label = ifelse(type == "body_site", name,
                                      ifelse(HMDB.Class == "Indoles and derivatives", 
                                             HMDB.Name, ""))),
                   repel = TRUE,
                   size = 3) +
    # Translated comment.
    theme_graph() +
    theme(
      legend.position = "bottom",
      plot.margin = margin(1, 1, 1, 1, "cm")
    ) +
    labs(
      edge_width = "关联强度",
      edge_alpha = "关联强度",
      color = "节点类型",
      size = "节点类型"
    )
}

# Translated comment.
# Translated comment.
example_data <- matrix(
  c(1, 0, 1, 0,
    0, 1, 1, 1,
    1, 0, 0, 1),
  nrow = 3,
  byrow = TRUE,
  dimnames = list(
    c("Metabolite1", "Metabolite2", "Metabolite3"),
    c("gut", "oral", "skin", "nasal")
  )
)

# Translated comment.
graph <- metabolite_to_network(four_site_R2,metabolite_annotation)

# Translated comment.
plot_metabolite_network(graph)