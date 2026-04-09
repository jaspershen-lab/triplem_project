rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")
setwd("1_code/4_site_merge/")
library(tidymass)
library(tidyverse)
library(readxl)
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
metabolite_annotation<-read_excel("../../3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx")
metabolite_annotation<-subset(metabolite_annotation,!(HMDB.Name=="NA"))
metabolite_annotation_micro<-subset(metabolite_annotation,HMDB.Source.Microbial=="TRUE")
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



##
##adjust BMI, sex, and IRIS, ethnicity
library(tidyverse)
library(ggpubr)
library(rstatix)

skin_expression_data <-
  extract_expression_data(skin_object) %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

library(plyr)

skin_sample_info <-
  skin_object@sample_info

#######adjust BMI, sex, and IRIS, ethnicity
skin_expression_data <-
  lm_adjust(expression_data = skin_expression_data,
            sample_info = skin_sample_info,
            threads = 3)

skin_temp_object <- skin_object
skin_temp_object@expression_data <- skin_expression_data



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



##
##adjust BMI, sex, and IRIS, ethnicity
library(tidyverse)
library(ggpubr)
library(rstatix)

nasal_expression_data <-
  extract_expression_data(nasal_object) %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

library(plyr)

nasal_sample_info <-
  nasal_object@sample_info

#######adjust BMI, sex, and IRIS, ethnicity
nasal_expression_data <-
  lm_adjust(expression_data = nasal_expression_data,
            sample_info = nasal_sample_info,
            threads = 3)

nasal_temp_object <- nasal_object
nasal_temp_object@expression_data <- nasal_expression_data

#######adjust BMI, sex, and IRIS, ethnicity
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

# Load required packages.
library(Hmisc)      # Calculatecorrelation.
library(igraph)     # Remove incomplete cases.
library(ggplot2)    # Plot.
library(reshape2)   # data.
library(tidyverse)  # Process the data.

# dataunder data framein:.
# Gut_data: gutmicrobiome data.
# Oral_data: oralmicrobiome data.
# Skin_data: skinmicrobiome data.
# Nasal_data: nasalmicrobiome data.
# Metabolites: metabolite data.

# : CalculatecorrelationFiltersignificantresults.
calculate_correlations <- function(microbiome_data,taxdata, metabolite_data,metabolite_anno, site, threshold = 0.3, p_value = 0.05) {
  # Ensure sample order matches.
  common_samples <- intersect(rownames(microbiome_data), rownames(metabolite_data))
  
  # data.
  microbiome_subset <- microbiome_data[common_samples,]
  
  colnames(microbiome_subset)<-  paste(site,taxdata$Genus,sep = "_")
  
  metabolite_subset <- metabolite_data[common_samples,metabolite_anno$variable_id]
  
  
  colnames(metabolite_subset)<-  metabolite_anno$HMDB.Name
  # Calculatecorrelationmatrixp-values.
  cors <- matrix(NA, nrow = ncol(microbiome_subset), ncol = ncol(metabolite_subset))
  pvals <- matrix(NA, nrow = ncol(microbiome_subset), ncol = ncol(metabolite_subset))
  
  # Calculatecorrelation.
  for(i in 1:ncol(microbiome_subset)) {
    for(j in 1:ncol(metabolite_subset)) {
      test_result <- cor.test(microbiome_subset[,i], metabolite_subset[,j], 
                              method = "spearman", exact = FALSE)
      cors[i,j] <- test_result$estimate
      pvals[i,j] <- test_result$p.value
    }
  }
  
  # Setcolumn names.
  rownames(cors) <- colnames(microbiome_subset)
  colnames(cors) <- colnames(metabolite_subset)
  rownames(pvals) <- colnames(microbiome_subset)
  colnames(pvals) <- colnames(metabolite_subset)
  
  # Convertdata.
  cors_df <- melt(cors)
  pvals_df <- melt(pvals)
  
  # Mergecorrelation coefficientp-values.
  result_df <- data.frame(
    Microbiome = cors_df$Var1,
    Metabolite = cors_df$Var2,
    Correlation = cors_df$value,
    P_value = pvals_df$value,
    Site = site
  )
  
  # Filtersignificantresults.
  significant_cors <- result_df %>%
    filter(abs(Correlation) >= threshold & P_value <= p_value)
  
  return(significant_cors)
}

# For body siteCalculatecorrelation.
gut_cors <- calculate_correlations(t(gut_temp_object@expression_data),gut_temp_object@variable_info,t(metabolomics_temp_object@expression_data),metabolite_annotation, "Gut")
oral_cors <- calculate_correlations(t(oral_temp_object@expression_data),oral_temp_object@variable_info, t(metabolomics_temp_object@expression_data),metabolite_annotation ,"Oral")
skin_cors <- calculate_correlations(t(skin_temp_object@expression_data),skin_temp_object@variable_info ,t(metabolomics_temp_object@expression_data),metabolite_annotation, "Skin")
nasal_cors <- calculate_correlations(t(nasal_temp_object@expression_data),nasal_temp_object@variable_info, t(metabolomics_temp_object@expression_data),metabolite_annotation, "Nasal")

# Mergecorrelationresults.
all_cors <- rbind(gut_cors, oral_cors, skin_cors, nasal_cors)



correlation_data<-all_cors


# UseggraphCreatenetwork plot.
edges <- correlation_data %>%
  select(Microbiome, Metabolite, Correlation, Site)

# Create.
nodes <- data.frame(
  name = unique(c(correlation_data$Microbiome, correlation_data$Metabolite)),
  type = ifelse(unique(c(correlation_data$Microbiome, correlation_data$Metabolite)) %in% correlation_data$Microbiome, "Microbiome", "Metabolite"),
  site = ifelse(unique(c(correlation_data$Microbiome, correlation_data$Metabolite)) %in% correlation_data$Microbiome, 
                correlation_data$Site[match(unique(c(correlation_data$Microbiome, correlation_data$Metabolite)), correlation_data$Microbiome)], 
                "Metabolite")
)

# Createigraphfor .
graph <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
# Calculatedegree.
# Calculatedegree.
node_degrees <- degree(graph)
nodes$degree <- node_degrees[match(nodes$name, names(node_degrees))]

# Filtertop 20.
top_nodes <- nodes %>%
  group_by(site) %>%
  top_n(10, degree) %>%
  ungroup()

# Filtertop.
filtered_edges <- edges %>%
  filter(Microbiome %in% top_nodes$name | Metabolite %in% top_nodes$name)

# Create(Filterafterin ).
filtered_nodes <- nodes %>%
  filter(name %in% unique(c(filtered_edges$Microbiome, filtered_edges$Metabolite)))


# Getsitetop 10.
top_10_nodes <- filtered_nodes %>%
  group_by(site) %>%
  top_n(10, degree) %>%
  ungroup() %>%
  pull(name)

# Addwhether labels.
filtered_nodes$show_label <- filtered_nodes$name %in% top_10_nodes

# CreateFilterafterigraphfor .
filtered_graph <- graph_from_data_frame(d = filtered_edges, vertices = filtered_nodes, directed = FALSE)

# Calculatedegree.
filtered_degrees <- degree(filtered_graph)
V(filtered_graph)$degree <- filtered_degrees
V(filtered_graph)$show_label <- filtered_nodes$show_label
# Createggraph.

library(ggraph)
p<-ggraph(filtered_graph, layout = "stress") +
  # Add.
  geom_edge_link(aes(edge_alpha = abs(Correlation),
                     edge_width = abs(Correlation),
                     color = Correlation > 0),
                 show.legend = TRUE) +
  # Add.
  geom_node_point(aes(fill = site, 
                      size = degree,
                      shape = type),color="white") +
  # Addlabels.
  geom_node_text(aes(label = ifelse(show_label, name, "")), 
                 repel = TRUE, 
                 size = 2,
                 max.overlaps = 20)+
  # Set.
  scale_edge_color_manual(values = c("TRUE" = "#FF9999", "FALSE" = "#9999FF"),
                          name = "Correlation",
                          labels = c("Negative", "Positive")) +
  scale_fill_manual(values = c("Gut" = "#edd064",
                               "Oral" = "#a1d5b9",
                               "Skin" = "#f2ccac",
                               "Nasal" = "#a17db4",
                               "Metabolite" = "#FF9999")) +
  # Set.
  scale_size_continuous(range = c(2, 10), name = "Degree")+
  # Set.
  scale_shape_manual(values = c("Microbiome" = 21, "Metabolite" = 22)) +
  # Set.
  scale_edge_alpha(range = c(0.2, 0.8)) +
  scale_edge_width(range = c(0.3, 2)) +
  # Theme settings.
  theme_graph() +
  theme(legend.position = "none") +
  guides(
    color = guide_legend(title = "Type"),
    size = guide_legend(title = "Node Type"),
    shape = guide_legend(title = "Node Type"),
    edge_alpha = guide_legend(title = "Correlation Strength"),
    edge_width = guide_legend(title = "Correlation Strength")
  ) +
  labs(title = "")

setwd(r4projects::get_project_wd())
ggsave(p,
       filename = "4_manuscript/Figures/Figure_2/figure_s2_network.pdf",
       width = 10,
       height = 10)





##########


# Import.
library(igraph)
library(dplyr)

# Analyze.
# filtered_graphfiltered_nodes.

# 1. Extractimportance.
node_importance <- data.frame(
  name = V(filtered_graph)$name,
  site = V(filtered_graph)$site,
  type = V(filtered_graph)$type,
  degree = V(filtered_graph)$degree,
  betweenness = betweenness(filtered_graph),
  closeness = closeness(filtered_graph),
  eigen_centrality = eigen_centrality(filtered_graph)$vector
)

# 2. Summarizemicrobesmetabolites.
# Microbes(body sitesort).
important_microbes <- node_importance %>%
  filter(type == "Microbiome") %>%
  group_by(site) %>%
  arrange(desc(degree)) %>%
  slice_head(n = 10) %>%
  ungroup()

# Metabolites(sort).
important_metabolites <- node_importance %>%
  filter(type == "Metabolite") %>%
  arrange(desc(degree)) %>%
  slice_head(n = 20)



# 4. Summarizebody sitefor metabolitesmicrobes.
# and metabolites.

# ,Calculatemicrobesand metabolitescorrelation.
microbe_influence <- all_cors %>%
  dplyr::group_by(Microbiome, Site) %>%
  dplyr::summarize(
    avg_correlation = mean(abs(Correlation)),
    total_connections = n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Site, desc(avg_correlation))

# Getbody sitemicrobes(5).
top_influential_microbes <- microbe_influence %>%
  dplyr::group_by(Site) %>%
  dplyr::slice_head(n = 5) %>%
  dplyr::ungroup()

cat("\nTop 5 microbes with the strongest effect on metabolites for each body site:\n")
print(top_influential_microbes, n = nrow(top_influential_microbes))

# 5. body sitemetabolites.
metabolites_by_site <- all_cors %>%
  dplyr::group_by(Metabolite) %>%
  dplyr::summarize(
    sites = list(unique(Site)),
    site_count = length(unique(Site)),
    total_microbes = length(unique(Microbiome)),
    avg_correlation = mean(abs(Correlation)),
    .groups = "drop"
  ) %>%
  dplyr:: arrange(desc(site_count), desc(avg_correlation))

# Getbody sitemetabolites(2or morebody site).
shared_metabolites <- metabolites_by_site %>%
  dplyr::filter(site_count > 1) %>%
  dplyr::slice_head(n = 15)

cat("\nTop 15 metabolites influenced by microbes from multiple body sites:\n")
print(shared_metabolites, n = nrow(shared_metabolites))

# 6. Create.
library(ggplot2)
library(forcats)

# body sitemicrobes.
p1 <- important_microbes %>%
  group_by(site) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(name = fct_reorder(name, degree)) %>%
  ggplot(aes(x = degree, y = name, fill = site)) +
  geom_col() +
  scale_fill_manual(values = c("Gut" = "#edd064",
                               "Oral" = "#a1d5b9",
                               "Skin" = "#f2ccac",
                               "Nasal" = "#a17db4")) +
  facet_wrap(~site, scales = "free_y",nrow = 1) +
  labs(title = "",
       x = "Degree",
       y = "") +
  theme_bw()

# Saveresults.
ggsave("4_manuscript/Figures/Figure_2/Figure_s2_important_microbes.pdf", 
       p1, width = 15, height = 4)

# metabolites.
p2 <- important_metabolites %>%
  slice_head(n = 10) %>%
  mutate(name = fct_reorder(name, degree)) %>%
  ggplot(aes(x = degree, y = name, fill = "Metabolite")) +
  geom_col(fill = "#FF9999") +
  labs(title = "",
       x = "Degree",
       y = "") +
  theme_bw()

ggsave("4_manuscript/Figures/Figure_2/Figure_s2_important_Metabolite.pdf", 
       p2, width = 4, height = 4)