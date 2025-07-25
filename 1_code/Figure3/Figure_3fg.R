rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")
library(ggbeeswarm)
setwd("1_code/4_site_merge/")
##### 合并 gut 和 oral 的GBDT数据
gut_GBDT_results <- readRDS("../../3_data_analysis/gut_microbiome/GBDT/cross_section/gut_GBDT_results")
oral_GBDT_results <- readRDS("../../3_data_analysis/oral_microbiome/GBDT/cross_section/oral_GBDT_results")
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
  metabolite_annotation[, c("variable_id", "HMDB.Name", "HMDB.Source.Microbial")],
  by.x = "metabolite",
  by.y = "variable_id"
)


gut_oral_results_summary <- gut_oral_results_summary %>%
  mutate(
    group = case_when(
      gut_R2 > 0.05 & oral_R2 < 0.05 ~ "gut",
      oral_R2 > 0.05 & gut_R2 < 0.05 ~ "oral",
      gut_R2 > 0.05 & oral_R2 > 0.05 ~ "co-influence",
      TRUE ~ "none"  # 这个处理其他情况，比如两个都是0的情况
    )
  )
# Create a subset for labeled points
labeled_data <- subset(gut_oral_results_summary, gut_R2 > 0.1 &
                         oral_R2 > 0.1)

library(ggrepel)

plot <-
  ggplot(
    gut_oral_results_summary,
    aes(
      x = gut_R2,
      y = oral_R2,
      fill = group,
      shape = HMDB.Source.Microbial
    )
  ) +
  geom_point(size = 6, color = "white") +
  scale_shape_manual(values = c(21, 22, 24)) +
  scale_fill_manual(values = c("#f5ccca", "#edd064", "grey70", "#a1d5b9")) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 4)),
         shape = guide_legend(override.aes = list(fill = "grey50", size = 4))) +
  theme_bw() +
  ylim(c(0, 0.3)) + geom_text_repel(
    data = labeled_data,
    aes(label = HMDB.Name),
    box.padding = 0.5,
    # Adjust label distance from points
    max.overlaps = 10,
    # Allowable label overlaps
    force = 10,
    # Increase label dispersion force
    segment.color = "grey50"
  ) + theme(
    legend.position = "right",
    axis.text.x = element_text(colour = "black", size = 14),
    axis.text.y = element_text(size = 14, face = "plain"),
    axis.title.y = element_text(size = 14, face = "plain"),
    axis.title.x = element_text(size = 14, face = "plain"),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
  )

plot

ggsave(
  plot,
  filename = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_3/figure_3f.pdf"
  ),
  width = 8,
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

gut_GBDT_results <- readRDS("../../3_data_analysis/gut_microbiome/GBDT/cross_section/gut_GBDT_results")
gut_GBDT_results_R2 <- gut_GBDT_results$summary[, c(1, 2)]
colnames(gut_GBDT_results_R2) <- c("metabolite", "gut")


gut_GBDT_results_R2$HMDB.ID <- metabolite_annotation$HMDB.ID


#

gut_GBDT_results_R2 <- subset(gut_GBDT_results_R2, gut > 0.05)

library(metpath)
library(tidyverse)
data("hmdb_pathway", package = "metpath")
data("query_id_hmdb", package = "metpath")

#get the class of pathways
pathway_class =
  hmdb_pathway@pathway_class

remain_idx = which(unlist(pathway_class) == "Metabolic;primary_pathway")

hmdb_pathway =
  hmdb_pathway[remain_idx]

compound_list <- hmdb_pathway@compound_list

# 需要筛选的 HMDB.ID 列表
selected_hmdb_ids <- metabolite_annotation$HMDB.ID

# 筛选 compound_list 中每个 data.frame 的元素
filtered_compound_list <- lapply(compound_list, function(df) {
  df[df$HMDB.ID %in% selected_hmdb_ids, ]
})

# 查看筛选结果
filtered_compound_list

hmdb_pathway@compound_list <- filtered_compound_list


result =
  enrich_hmdb(
    query_id = gut_GBDT_results_R2$HMDB.ID,
    query_type = "compound",
    id_type = "HMDB",
    pathway_database = hmdb_pathway,
    only_primary_pathway = TRUE,
    p_cutoff = 0.05,
    p_adjust_method = "none",
    threads = 3
  )

enrich_results <- result@result

enrich_results <- subset(enrich_results, mapped_number >= 5 &
                           mapped_percentage >= 50)

enrich_results_gut <- enrich_results


#######


oral_GBDT_results <- readRDS("../../3_data_analysis/oral_microbiome/GBDT/cross_section/oral_GBDT_results")
oral_GBDT_results_R2 <- oral_GBDT_results$summary[, c(1, 2)]
colnames(oral_GBDT_results_R2) <- c("metabolite", "oral")
oral_GBDT_results_R2$HMDB.ID <- metabolite_annotation$HMDB.ID


#

oral_GBDT_results_R2 <- subset(oral_GBDT_results_R2, oral > 0.05)

library(metpath)
library(tidyverse)
data("hmdb_pathway", package = "metpath")
data("query_id_hmdb", package = "metpath")

#get the class of pathways
pathway_class =
  hmdb_pathway@pathway_class

remain_idx = which(unlist(pathway_class) == "Metabolic;primary_pathway")

hmdb_pathway =
  hmdb_pathway[remain_idx]

compound_list <- hmdb_pathway@compound_list

# 需要筛选的 HMDB.ID 列表
selected_hmdb_ids <- metabolite_annotation$HMDB.ID

# 筛选 compound_list 中每个 data.frame 的元素
filtered_compound_list <- lapply(compound_list, function(df) {
  df[df$HMDB.ID %in% selected_hmdb_ids, ]
})

# 查看筛选结果
filtered_compound_list

hmdb_pathway@compound_list <- filtered_compound_list


result =
  enrich_hmdb(
    query_id = oral_GBDT_results_R2$HMDB.ID,
    query_type = "compound",
    id_type = "HMDB",
    pathway_database = hmdb_pathway,
    only_primary_pathway = TRUE,
    p_cutoff = 0.05,
    p_adjust_method = "none",
    threads = 3
  )


enrich_results <- result@result

enrich_results <- subset(enrich_results, mapped_number >= 3 &
                           mapped_percentage >= 30)

enrich_results_oral <- enrich_results


#######


skin_GBDT_results <- readRDS("../../3_data_analysis/skin_microbiome/GBDT/cross_section/skin_GBDT_results")
skin_GBDT_results_R2 <- skin_GBDT_results$summary[, c(1, 2)]
colnames(skin_GBDT_results_R2) <- c("metabolite", "skin")
skin_GBDT_results_R2$HMDB.ID <- metabolite_annotation$HMDB.ID


#

skin_GBDT_results_R2 <- subset(skin_GBDT_results_R2, skin > 0.05)

library(metpath)
library(tidyverse)
data("hmdb_pathway", package = "metpath")
data("query_id_hmdb", package = "metpath")

#get the class of pathways
pathway_class =
  hmdb_pathway@pathway_class

remain_idx = which(unlist(pathway_class) == "Metabolic;primary_pathway")

hmdb_pathway =
  hmdb_pathway[remain_idx]

compound_list <- hmdb_pathway@compound_list

# 需要筛选的 HMDB.ID 列表
selected_hmdb_ids <- metabolite_annotation$HMDB.ID

# 筛选 compound_list 中每个 data.frame 的元素
filtered_compound_list <- lapply(compound_list, function(df) {
  df[df$HMDB.ID %in% selected_hmdb_ids, ]
})

# 查看筛选结果
filtered_compound_list

hmdb_pathway@compound_list <- filtered_compound_list

result =
  enrich_hmdb(
    query_id = skin_GBDT_results_R2$HMDB.ID,
    query_type = "compound",
    id_type = "HMDB",
    pathway_database = hmdb_pathway,
    only_primary_pathway = TRUE,
    p_cutoff = 0.05,
    p_adjust_method = "none",
    threads = 3
  )

enrich_results <- result@result

enrich_results <- subset(enrich_results, mapped_number >= 3 &
                           mapped_percentage >= 30)

enrich_results_skin <- enrich_results



#######


nasal_GBDT_results <- readRDS("../../3_data_analysis/nasal_microbiome/GBDT/cross_section/nasal_GBDT_results")
nasal_GBDT_results_R2 <- nasal_GBDT_results$summary[, c(1, 2)]
colnames(nasal_GBDT_results_R2) <- c("metabolite", "nasal")
nasal_GBDT_results_R2$HMDB.ID <- metabolite_annotation$HMDB.ID


#

nasal_GBDT_results_R2 <- subset(nasal_GBDT_results_R2, nasal > 0.05)

library(metpath)
library(tidyverse)
data("hmdb_pathway", package = "metpath")
data("query_id_hmdb", package = "metpath")

#get the class of pathways
pathway_class =
  hmdb_pathway@pathway_class

remain_idx = which(unlist(pathway_class) == "Metabolic;primary_pathway")

hmdb_pathway =
  hmdb_pathway[remain_idx]

compound_list <- hmdb_pathway@compound_list

# 需要筛选的 HMDB.ID 列表
selected_hmdb_ids <- metabolite_annotation$HMDB.ID

# 筛选 compound_list 中每个 data.frame 的元素
filtered_compound_list <- lapply(compound_list, function(df) {
  df[df$HMDB.ID %in% selected_hmdb_ids, ]
})

# 查看筛选结果
filtered_compound_list

hmdb_pathway@compound_list <- filtered_compound_list

result =
  enrich_hmdb(
    query_id = nasal_GBDT_results_R2$HMDB.ID,
    query_type = "compound",
    id_type = "HMDB",
    pathway_database = hmdb_pathway,
    only_primary_pathway = TRUE,
    p_cutoff = 0.05,
    p_adjust_method = "none",
    threads = 3
  )

enrich_results <- result@result

enrich_results <- subset(enrich_results, mapped_number >= 3 &
                           mapped_percentage >= 40)

enrich_results_nasal <- enrich_results


#########合并四个位点代谢物通路富集结果生成网络图


enrich_results_gut$site <- "gut"
enrich_results_oral$site <- "oral"
enrich_results_skin$site <- "skin"
enrich_results_nasal$site <- "nasal"

enrich_results_merge <- rbind(
  enrich_results_gut,
  enrich_results_oral,
  enrich_results_skin,
  enrich_results_nasal
)


#
# enrich_results_merge <-
#   enrich_results_merge[, c("pathway_name", "site", "mapped_number")]
#
# edges <- data.frame(from = enrich_results_merge$site, to = enrich_results_merge$pathway_name)
# edges$color <- edges$from
# vertices <- data.frame(
#   name = unique(
#     c(
#       enrich_results_merge$site,
#       enrich_results_merge$pathway_name
#     )
#   ),
#   type = c(rep("site", 4), rep("pathway", 26)),
#   size = c(20, 20, 20, 20, rep(5, 26))
# )
#
# vertices$color <- 'black'
#
# vertices$show_name <- vertices$name
#
# ggraph_data <- igraph::graph_from_data_frame(d = edges,
#                                              vertices = vertices,
#                                              directed = T)
#
# ggraph_data <- igraph::graph_from_data_frame(d = edges,
#                                              vertices = vertices,
#                                              directed = T)
#
# #加载包
# library(ggplot2)
# library(RColorBrewer)
# library(tidyverse)
# library(ggsci)
# library(igraph)
# library(ggraph)
#
#
# layout_in_circles <- function(g, group = 1) {
#   layout <- lapply(split(V(g), group), function(x) {
#     layout_in_circle(induced_subgraph(g, x))
#   })
#   layout <- Map(`*`, layout, seq_along(layout))
#   x <- matrix(0, nrow = vcount(g), ncol = 2)
#   split(x, group) <- layout
#   x
# }
#
# ggraph(ggraph_data, layout = layout_in_circles(ggraph_data, group = vertices$type)) +
#   geom_edge_link(
#     aes(edge_color = color),
#     edge_alpha = 0.8,
#     arrow = arrow(
#       type = "closed",
#       angle = 15,
#       ends = "last",
#       length = unit(0.1, "inches")
#     )
#   ) +
#   geom_node_point(aes(
#     shape = type,
#     fill = type,
#     size = size,
#     color = color
#   ), alpha = .75) +
#
#   geom_node_text(aes(x = x, y = y + 0.1 , label  = show_name), size = 3) +
#   theme_void() +
#   scale_color_manual(values = c('black', pal_startrek()(2)[c(2, 1)])) +
#   scale_fill_manual(values = pal_startrek()(2)[c(2, 1)]) +
#   scale_size_continuous(range = c(2, 4)) +
#   scale_shape_manual(values = c(22, 21)) +
#   coord_equal() +
#   theme(legend.position = 'none') + scale_edge_color_manual(values = body_site_color)



# enrich_results_merge <-
#   enrich_results_merge[, c("pathway_name", "site", "mapped_number", "mapped_id")]
# 
# save(enrich_results_merge, file = "enrich_results_merge.rda")
load("enrich_results_merge.rda")
load("hmdb_ms1.rda")

edges1 <-
  data.frame(from = enrich_results_merge$site, to = enrich_results_merge$pathway_name) %>%
  dplyr::mutate(edge_class = from)

edges2 <-
  seq_len(nrow(enrich_results_merge)) %>%
  purrr::map(function(i) {
    data.frame(
      from = enrich_results_merge$pathway_name[i],
      to = stringr::str_split(enrich_results_merge$mapped_id[i], ";")[[1]]
    )
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::mutate(edge_class = "pathway")

edge_data <-
  rbind(edges1, edges2)

node_data <-
  rbind(
    data.frame(
      node = unique(enrich_results_merge$site),
      node_class = unique(enrich_results_merge$site)
    ),
    data.frame(
      node = unique(enrich_results_merge$pathway_name),
      node_class = "pathway"
    ),
    data.frame(
      node = stringr::str_split(enrich_results_merge$mapped_id, ";") %>%
        unlist() %>%
        unique(),
      node_class = "metabolite"
    )
  )

node_data <-
  node_data %>%
  dplyr::left_join(hmdb_ms1@spectra.info[, c("HMDB.ID", "Compound.name")], by = c("node" = "HMDB.ID"))  %>%
  dplyr::mutate(node_name =
                  case_when(!is.na(Compound.name) ~ Compound.name, TRUE ~ node))


graph_data <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%
  tidygraph::activate(what = "nodes") %>%
  dplyr::mutate(degree = tidygraph::centrality_degree())


g <-
  graph_data

igraph::V(g)$type <- igraph::bipartite_mapping(g)$type

coords <-
  ggraph::create_layout(g, layout = "bipartite")

coords$index = 1:nrow(coords)

coords$x <-
  coords$x + 1

coords$y[coords$node_class == "microbiome"] <- 1
coords$y[coords$node_class == "pathway"] <- 2
coords$y[coords$node_class == "metabolite"] <- 3

# temp_coords <-
#   coords %>%
#   dplyr::filter(class != "metabolite") %>%
#   dplyr::arrange(x)
#
# temp_coords$x <-
#   seq(min(coords$x) + 0.1 * max(coords$x),
#       0.9 * max(coords$x),
#       length.out = nrow(temp_coords))
#
# coords[coords$class != "metabolite", ] <-
#   temp_coords

coords <-
  coords %>%
  dplyr::arrange(index)

coords <-
  coords %>%
  dplyr::mutate(x1 = y, y1 = x) %>%
  dplyr::select(-c(x, y)) %>%
  dplyr::mutate(x = x1, y = y1)

# if (circle) {
#   coords <-
#     coords %>%
#     dplyr::mutate(
#       theta = y / (max(y) + 1) * 2 * pi,
#       r = x + 1,
#       x = r * cos(theta),
#       y = r * sin(theta)
#     )
# }

my_graph <-
  ggraph::create_layout(
    graph = g,
    layout = "manual",
    x = coords$x,
    y = coords$y
    # node.position = coords
  )

plot <-
  ggraph::ggraph(my_graph, layout = 'bipartite') +
  ggraph::geom_edge_diagonal(
    strength = 1,
    aes(color = edge_class),
    edge_width = 0.5,
    alpha = 0.8,
    show.legend = FALSE
  ) +
  ggraph::scale_edge_color_manual(values = c(site_color, pathway = ggsci::pal_aaas()(n =
                                                                                       3)[2])) +
  ggraph::geom_node_point(aes(size = degree, fill = node_class),
                          color = "black",
                          shape = 21) +
  scale_fill_manual(values = c(
    site_color,
    pathway = ggsci::pal_aaas()(n = 3)[2],
    metabolite = "black"
  )) +
  ggraph::geom_node_text(aes(x = x, y = y, label = node_name),
                         size = 4,
                         show.legend = FALSE) +
  ggraph::theme_graph() +
  scale_size_continuous(range = c(2, 10)) +
  theme(plot.background = element_blank(), panel.background = element_blank())

plot

ggsave(
  plot,
  filename = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_3/figure_3g.pdf"
  ),
  width = 12,
  height = 6
)
