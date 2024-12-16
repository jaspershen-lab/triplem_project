rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")
library(readxl)

metabolite_annotation<-read_excel("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx")
setwd("1_code/4_site_merge/")

gut_GBDT_results<-readRDS("../../3_data_analysis/gut_microbiome/GBDT/cross_section/gut_GBDT_results")
gut_GBDT_results_R2<-gut_GBDT_results$summary[,c(1,2)]
colnames(gut_GBDT_results_R2)<-c("metabolite","gut")


gut_GBDT_results_R2$HMDB.ID<-metabolite_annotation$HMDB.ID


# 

gut_GBDT_results_R2<-subset(gut_GBDT_results_R2,gut>0.05)

library(metpath)
library(tidyverse)
data("hmdb_pathway", package = "metpath")
data("query_id_hmdb", package = "metpath")

#get the class of pathways
pathway_class = 
  metpath::pathway_class(hmdb_pathway) 

remain_idx = which(unlist(pathway_class) == "Metabolic;primary_pathway")


hmdb_pathway = 
  hmdb_pathway[remain_idx]

compound_list<-hmdb_pathway@compound_list





# 需要筛选的 HMDB.ID 列表
selected_hmdb_ids <- metabolite_annotation$HMDB.ID

# 筛选 compound_list 中每个 data.frame 的元素
filtered_compound_list <- lapply(compound_list, function(df) {
  df[df$HMDB.ID %in% selected_hmdb_ids, ]
})

# 查看筛选结果
filtered_compound_list

hmdb_pathway@compound_list<-filtered_compound_list


result = 
  enrich_hmdb(query_id = gut_GBDT_results_R2$HMDB.ID, 
              query_type = "compound", 
              id_type = "HMDB",
              pathway_database = hmdb_pathway,
              only_primary_pathway = TRUE,
              p_cutoff = 0.05, 
              p_adjust_method = "none", 
              threads = 3)


enrich_results<-result@result

enrich_results<-subset(enrich_results,mapped_number>=5&mapped_percentage>=50)

enrich_results_gut<-enrich_results




#######


oral_GBDT_results<-readRDS("../../3_data_analysis/oral_microbiome/GBDT/cross_section/oral_GBDT_results")
oral_GBDT_results_R2<-oral_GBDT_results$summary[,c(1,2)]
colnames(oral_GBDT_results_R2)<-c("metabolite","oral")
oral_GBDT_results_R2$HMDB.ID<-metabolite_annotation$HMDB.ID


# 

oral_GBDT_results_R2<-subset(oral_GBDT_results_R2,oral>0.05)

library(metpath)
library(tidyverse)
data("hmdb_pathway", package = "metpath")
data("query_id_hmdb", package = "metpath")

#get the class of pathways
pathway_class = 
  metpath::pathway_class(hmdb_pathway) 

remain_idx = which(unlist(pathway_class) == "Metabolic;primary_pathway")


hmdb_pathway = 
  hmdb_pathway[remain_idx]

compound_list<-hmdb_pathway@compound_list





# 需要筛选的 HMDB.ID 列表
selected_hmdb_ids <- metabolite_annotation$HMDB.ID

# 筛选 compound_list 中每个 data.frame 的元素
filtered_compound_list <- lapply(compound_list, function(df) {
  df[df$HMDB.ID %in% selected_hmdb_ids, ]
})

# 查看筛选结果
filtered_compound_list

hmdb_pathway@compound_list<-filtered_compound_list


result = 
  enrich_hmdb(query_id = oral_GBDT_results_R2$HMDB.ID, 
              query_type = "compound", 
              id_type = "HMDB",
              pathway_database = hmdb_pathway,
              only_primary_pathway = TRUE,
              p_cutoff = 0.05, 
              p_adjust_method = "none", 
              threads = 3)


enrich_results<-result@result

enrich_results<-subset(enrich_results,mapped_number>=3&mapped_percentage>=30)

enrich_results_oral<-enrich_results


#######


skin_GBDT_results<-readRDS("../../3_data_analysis/skin_microbiome/GBDT/cross_section/skin_GBDT_results")
skin_GBDT_results_R2<-skin_GBDT_results$summary[,c(1,2)]
colnames(skin_GBDT_results_R2)<-c("metabolite","skin")
skin_GBDT_results_R2$HMDB.ID<-metabolite_annotation$HMDB.ID


# 

skin_GBDT_results_R2<-subset(skin_GBDT_results_R2,skin>0.05)

library(metpath)
library(tidyverse)
data("hmdb_pathway", package = "metpath")
data("query_id_hmdb", package = "metpath")

#get the class of pathways
pathway_class = 
  metpath::pathway_class(hmdb_pathway) 

remain_idx = which(unlist(pathway_class) == "Metabolic;primary_pathway")


hmdb_pathway = 
  hmdb_pathway[remain_idx]

compound_list<-hmdb_pathway@compound_list





# 需要筛选的 HMDB.ID 列表
selected_hmdb_ids <- metabolite_annotation$HMDB.ID

# 筛选 compound_list 中每个 data.frame 的元素
filtered_compound_list <- lapply(compound_list, function(df) {
  df[df$HMDB.ID %in% selected_hmdb_ids, ]
})

# 查看筛选结果
filtered_compound_list

hmdb_pathway@compound_list<-filtered_compound_list


result = 
  enrich_hmdb(query_id = skin_GBDT_results_R2$HMDB.ID, 
              query_type = "compound", 
              id_type = "HMDB",
              pathway_database = hmdb_pathway,
              only_primary_pathway = TRUE,
              p_cutoff = 0.05, 
              p_adjust_method = "none", 
              threads = 3)


enrich_results<-result@result

enrich_results<-subset(enrich_results,mapped_number>=3&mapped_percentage>=30)

enrich_results_skin<-enrich_results






#######


nasal_GBDT_results<-readRDS("../../3_data_analysis/nasal_microbiome/GBDT/cross_section/nasal_GBDT_results")
nasal_GBDT_results_R2<-nasal_GBDT_results$summary[,c(1,2)]
colnames(nasal_GBDT_results_R2)<-c("metabolite","nasal")
nasal_GBDT_results_R2$HMDB.ID<-metabolite_annotation$HMDB.ID


# 

nasal_GBDT_results_R2<-subset(nasal_GBDT_results_R2,nasal>0.05)

library(metpath)
library(tidyverse)
data("hmdb_pathway", package = "metpath")
data("query_id_hmdb", package = "metpath")

#get the class of pathways
pathway_class = 
  metpath::pathway_class(hmdb_pathway) 

remain_idx = which(unlist(pathway_class) == "Metabolic;primary_pathway")


hmdb_pathway = 
  hmdb_pathway[remain_idx]

compound_list<-hmdb_pathway@compound_list





# 需要筛选的 HMDB.ID 列表
selected_hmdb_ids <- metabolite_annotation$HMDB.ID

# 筛选 compound_list 中每个 data.frame 的元素
filtered_compound_list <- lapply(compound_list, function(df) {
  df[df$HMDB.ID %in% selected_hmdb_ids, ]
})

# 查看筛选结果
filtered_compound_list

hmdb_pathway@compound_list<-filtered_compound_list


result = 
  enrich_hmdb(query_id = nasal_GBDT_results_R2$HMDB.ID, 
              query_type = "compound", 
              id_type = "HMDB",
              pathway_database = hmdb_pathway,
              only_primary_pathway = TRUE,
              p_cutoff = 0.05, 
              p_adjust_method = "none", 
              threads = 3)


enrich_results<-result@result

enrich_results<-subset(enrich_results,mapped_number>=3&mapped_percentage>=40)

enrich_results_nasal<-enrich_results







#########合并四个位点代谢物通路富集结果生成网络图


enrich_results_gut$site<-"gut"
enrich_results_oral$site<-"oral"
enrich_results_skin$site<-"skin"
enrich_results_nasal$site<-"nasal"



enrich_results_merge<-rbind(enrich_results_gut,enrich_results_oral,enrich_results_skin,enrich_results_nasal)

enrich_results_merge<-enrich_results_merge[,c("pathway_name","site","mapped_number")]




edges<-data.frame(from=enrich_results_merge$site,to=enrich_results_merge$pathway_name)
edges$color<-edges$from
vertices<-data.frame(name=unique(c(enrich_results_merge$site,enrich_results_merge$pathway_name)),type=c(rep("site",4),rep("pathway",26)),size=c(20,20,20,20,rep(5,26)))



vertices$color <- 'black'

vertices$show_name <- vertices$name

ggraph_data <- igraph::graph_from_data_frame(d = edges,vertices = vertices,directed = T )


ggraph_data <- igraph::graph_from_data_frame(d = edges,
                                             vertices = vertices,
                                             directed = T
)



#加载包
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(ggsci)
library(igraph)
library(ggraph)

ggraph(graph = ggraph_data,
       layout = 'lgl',
       circular = F)+
  geom_edge_link(aes(edge_color = color),
                 edge_alpha = 0.8,
                 arrow = arrow(type = "closed",
                               angle = 15,
                               ends = "last",
                               length = unit(0.1, "inches")))+
  geom_node_point(aes(shape = type,
                      fill = type,
                      size = size,
                      color = color),
                  alpha = .75)+
  
  geom_node_text(aes(x = x,
                     y = y-0.5 ,
                     label  = show_name),
                 size = 3)+
  theme_void()+
  scale_color_manual(values = c('black',pal_startrek()(2)[c(2,1)]))+
  scale_fill_manual(values = pal_startrek()(2)[c(2,1)])+
  scale_size_continuous(range = c(2,4))+
  scale_shape_manual(values = c(22,21))+
  coord_equal()+
  theme(legend.position = 'none')+scale_edge_color_manual(values = body_site_color)
