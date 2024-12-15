rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")
library(readxl)

metabolite_annotation<-read_excel("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx")
setwd("1_code/4_site_merge/")
##  814个代谢物的热图


gut_GBDT_results<-readRDS("../../3_data_analysis/gut_microbiome/GBDT/cross_section/gut_GBDT_results")
oral_GBDT_results<-readRDS("../../3_data_analysis/oral_microbiome/GBDT/cross_section/oral_GBDT_results")
skin_GBDT_results<-readRDS("../../3_data_analysis/skin_microbiome/GBDT/cross_section/skin_GBDT_results")
nasal_GBDT_results<-readRDS("../../3_data_analysis/nasal_microbiome/GBDT/cross_section/nasal_GBDT_results")
gut_GBDT_results_R2<-gut_GBDT_results$summary[,c(1,2)]
colnames(gut_GBDT_results_R2)<-c("metabolite","gut")


oral_GBDT_results_R2<-oral_GBDT_results$summary[,c(1,2)]
colnames(oral_GBDT_results_R2)<-c("metabolite","oral")


skin_GBDT_results_R2<-skin_GBDT_results$summary[,c(1,2)]
colnames(skin_GBDT_results_R2)<-c("metabolite","skin")


nasal_GBDT_results_R2<-nasal_GBDT_results$summary[,c(1,2)]
colnames(nasal_GBDT_results_R2)<-c("metabolite","nasal")



four_site_GBDT_R2<-cbind(gut_GBDT_results_R2,oral_GBDT_results_R2$oral,skin_GBDT_results_R2$skin,nasal_GBDT_results_R2$nasal)
colnames(four_site_GBDT_R2)<-c("metabolite","gut","oral","skin","nasal")


rownames(four_site_GBDT_R2)<-four_site_GBDT_R2$metabolite
four_site_GBDT_R2<-four_site_GBDT_R2[,-1]

four_site_GBDT_R2[four_site_GBDT_R2 < 0.05] <- 0

four_site_GBDT_R2<-four_site_GBDT_R2[rowSums(four_site_GBDT_R2)>0,]












# 加载需要的包
library(dplyr)
library(ggplot2)
library(stats)

metabolite_class_enrichment <- function(
    significant_metabolites,  # 显著代谢物的向量
    all_metabolites_df,      # 包含所有代谢物及其class信息的数据框
    class_column,            # class信息的列名
    metabolite_column,       # 代谢物ID/名称的列名
    alpha = 0.05            # 显著性水平
) {
  # 获取总体代谢物数量
  N <- nrow(all_metabolites_df)
  
  # 获取显著代谢物数量
  n <- length(significant_metabolites)
  
  # 对每个class进行分析
  results <- all_metabolites_df %>%
    dplyr::group_by(!!sym(class_column)) %>%
    dplyr::summarise(
      Total_in_class = n(),
      Significant_in_class = sum(!!sym(metabolite_column) %in% significant_metabolites)
    ) %>%
    mutate(
      Expected_by_chance = (Total_in_class * n) / N,
      Fold_enrichment = (Significant_in_class/n)/(Total_in_class/N),
      # 计算超几何分布的p值
      P_value = phyper(
        Significant_in_class - 1,
        Total_in_class,
        N - Total_in_class,
        n,
        lower.tail = FALSE
      )
    )
  
  # FDR校正
  results$FDR <- p.adjust(results$P_value, method = "BH")
  
  # 按p值排序
  results <- results %>% arrange(P_value)
  
  return(results)
}
HMDB.Class%in%c("Benzene and substituted derivatives","Carboxylic acids and derivatives","Fatty Acyls","Glycerophospholipids","Indoles and derivatives","Organooxygen compounds","Steroids and steroid derivatives")

four_site_GBDT_R2_gut<-subset(four_site_GBDT_R2,gut>0.05)
four_site_GBDT_R2_oral<-subset(four_site_GBDT_R2,oral>0.05)
four_site_GBDT_R2_skin<-subset(four_site_GBDT_R2,skin>0.05)
four_site_GBDT_R2_nasal<-subset(four_site_GBDT_R2,nasal>0.05)
metabolite_annotation_set<-subset(metabolite_annotation,variable_id%in%rownames(four_site_GBDT_R2))

results = metabolite_class_enrichment(
  significant_metabolites=rownames(four_site_GBDT_R2_gut),
  all_metabolites_df=metabolite_annotation_set,
  class_column='HMDB.Class',
  metabolite_column='variable_id'
)



results<-subset(results,P_value<1)

results$site<-"gut"

results_gut<-results



results = metabolite_class_enrichment(
  significant_metabolites=rownames(four_site_GBDT_R2_oral),
  all_metabolites_df=metabolite_annotation_set,
  class_column='HMDB.Class',
  metabolite_column='variable_id'
)



results<-subset(results,P_value<1)

results$site<-"oral"

results_oral<-results


results = metabolite_class_enrichment(
  significant_metabolites=rownames(four_site_GBDT_R2_skin),
  all_metabolites_df=metabolite_annotation_set,
  class_column='HMDB.Class',
  metabolite_column='variable_id'
)



results<-subset(results,P_value<1)

results$site<-"skin"

results_skin<-results


results = metabolite_class_enrichment(
  significant_metabolites=rownames(four_site_GBDT_R2_nasal),
  all_metabolites_df=metabolite_annotation_set,
  class_column='HMDB.Class',
  metabolite_column='variable_id'
)



results<-subset(results,P_value<1)

results$site<-"nasal"

results_nasal<-results


results_all<-rbind(results_gut,results_oral,results_skin,results_nasal)


results_all<-subset(results_all,HMDB.Class%in%c("Benzene and substituted derivatives","Carboxylic acids and derivatives","Fatty Acyls","Glycerophospholipids","Indoles and derivatives","Organooxygen compounds","Steroids and steroid derivatives")
)




### 

ggplot(results_all, aes(x = HMDB.Class, y = Significant_in_class)) +
  geom_bar(stat = "identity", aes(fill=site)) +
  theme_minimal() +facet_wrap(~site,nrow = 1)+coord_flip()+theme_light()+
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14,family = "Helvetica"),
    axis.title = element_text(size = 14,family = "Helvetica"),
    axis.text.x = element_text(family = "Helvetica") , # 如果组名较长，可以倾斜x轴标签
    axis.ticks.length = unit(0.25, "cm"),  # 增加刻度线长度
    axis.ticks = element_line(linewidth = 0.8)  # 增加刻度线粗细
  )+scale_fill_manual(values = body_site_color)