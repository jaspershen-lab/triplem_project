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

gut_GBDT_results_R2<-subset(gut_GBDT_results_R2,gut>0.01)

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

enrich_results<-subset(enrich_results,mapped_num>=5&mapped_percentage>=70)

saveRDS(enrich_results,"gut_enrich_results")
