rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")
source("1_code/mantel_Procrustes_code.R")
library(tidyverse)
library(tidymass)
library(readxl)


## load data
### for gut microbiome

load(file = "2_data/gut_microbiome/data_preparation/expression_data")
load(file = "2_data/gut_microbiome/data_preparation/sample_info")
load(file = "2_data/gut_microbiome/data_preparation/variable_info")


gut_microbiome_table<-expression_data
gut_microbiome_metadata<-sample_info
gut_microbiome_tax<-variable_info

# 标准化微生物组数据

gut_microbiome_metadata <- gut_microbiome_metadata %>% filter(!is.na(Gender))
gut_microbiome_metadata <- gut_microbiome_metadata %>% filter(!is.na(IRIS))
gut_microbiome_metadata <- gut_microbiome_metadata %>% filter(!is.na(BMI))


gut_microbiome_table<-gut_microbiome_table[,gut_microbiome_metadata$sample_id]


gut_microbiome_table<-microbiome_genus_filt(gut_microbiome_table,gut_microbiome_tax,0.1)

gut_microbiome_table_clr <- t(data.frame(clr(gut_microbiome_table)))


gut_microbiome_table <-lm_adjust(expression_data = gut_microbiome_table_clr,
                                 sample_info = gut_microbiome_metadata,
                                 threads = 3)

gut_microbiome_table<-data.frame(t(gut_microbiome_table))
### for oral microbiome

load(file = "2_data/oral_microbiome/data_preparation/expression_data")
load(file = "2_data/oral_microbiome/data_preparation/sample_info")
load(file = "2_data/oral_microbiome/data_preparation/variable_info")


oral_microbiome_table<-expression_data
oral_microbiome_metadata<-sample_info
oral_microbiome_tax<-variable_info

# 标准化微生物组数据

oral_microbiome_metadata <- oral_microbiome_metadata %>% filter(!is.na(Gender))
oral_microbiome_metadata <- oral_microbiome_metadata %>% filter(!is.na(IRIS))
oral_microbiome_metadata <- oral_microbiome_metadata %>% filter(!is.na(BMI))


oral_microbiome_table<-oral_microbiome_table[,oral_microbiome_metadata$sample_id]


oral_microbiome_table<-microbiome_genus_filt(oral_microbiome_table,oral_microbiome_tax,0.1)

oral_microbiome_table_clr <- t(data.frame(clr(oral_microbiome_table)))


oral_microbiome_table <-lm_adjust(expression_data = oral_microbiome_table_clr,
                                  sample_info = oral_microbiome_metadata,
                                  threads = 3)

oral_microbiome_table<-data.frame(t(oral_microbiome_table))

### for skin microbiome

load(file = "2_data/skin_microbiome/data_preparation/expression_data")
load(file = "2_data/skin_microbiome/data_preparation/sample_info")
load(file = "2_data/skin_microbiome/data_preparation/variable_info")


skin_microbiome_table<-expression_data
skin_microbiome_metadata<-sample_info
skin_microbiome_tax<-variable_info

# 标准化微生物组数据

skin_microbiome_metadata <- skin_microbiome_metadata %>% filter(!is.na(Gender))
skin_microbiome_metadata <- skin_microbiome_metadata %>% filter(!is.na(IRIS))
skin_microbiome_metadata <- skin_microbiome_metadata %>% filter(!is.na(BMI))


skin_microbiome_table<-skin_microbiome_table[,skin_microbiome_metadata$sample_id]


skin_microbiome_table<-microbiome_genus_filt(skin_microbiome_table,skin_microbiome_tax,0.1)

skin_microbiome_table_clr <- t(data.frame(clr(skin_microbiome_table)))


skin_microbiome_table <-lm_adjust(expression_data = skin_microbiome_table_clr,
                                  sample_info = skin_microbiome_metadata,
                                  threads = 3)

skin_microbiome_table<-data.frame(t(skin_microbiome_table))

### for nasal microbiome

load(file = "2_data/nasal_microbiome/data_preparation/expression_data")
load(file = "2_data/nasal_microbiome/data_preparation/sample_info")
load(file = "2_data/nasal_microbiome/data_preparation/variable_info")


nasal_microbiome_table<-expression_data
nasal_microbiome_metadata<-sample_info
nasal_microbiome_tax<-variable_info


# 标准化微生物组数据

nasal_microbiome_metadata <- nasal_microbiome_metadata %>% filter(!is.na(Gender))
nasal_microbiome_metadata <- nasal_microbiome_metadata %>% filter(!is.na(IRIS))
nasal_microbiome_metadata <- nasal_microbiome_metadata %>% filter(!is.na(BMI))


nasal_microbiome_table<-nasal_microbiome_table[,nasal_microbiome_metadata$sample_id]


nasal_microbiome_table<-microbiome_genus_filt(nasal_microbiome_table,nasal_microbiome_tax,0.1)

nasal_microbiome_table_clr <- t(data.frame(clr(nasal_microbiome_table)))


nasal_microbiome_table <-lm_adjust(expression_data = nasal_microbiome_table_clr,
                                   sample_info = nasal_microbiome_metadata,
                                   threads = 3)

nasal_microbiome_table<-data.frame(t(nasal_microbiome_table))


MM_adonis <- function(gut_microbiome_table, expression_data_metabolome) {
  library(compositions)
  library(vegan)
  # 与表达数据的交集
  
  shared_samples <- intersect(row.names(gut_microbiome_table), row.names(expression_data_metabolome))
  expression_data_metabolome <- expression_data_metabolome[shared_samples, ]
  gut_microbiome_table<-gut_microbiome_table[shared_samples,]
  
  # 计算距离并执行adonis分析
  metabolome_dist <- vegdist(expression_data_metabolome, method = "euclidean")
  results <- NULL
  for (i in colnames(gut_microbiome_table)) {
    adonis_result <- adonis2(metabolome_dist ~ gut_microbiome_table[, i], permutations = 1000,parallel=15)
    results <- rbind(results, adonis_result$`Pr(>F)`[1])
  }
  results<-data.frame(results)
  row.names(results) <- colnames(gut_microbiome_table)
  colnames(results)<-c("Pr(>F)")
  significant_results <- subset(results, `Pr(>F)` < 0.05)
  
  adonis_micro_sig <- gut_microbiome_table[, row.names(significant_results)]
  final_result <- adonis2(metabolome_dist ~ ., data = adonis_micro_sig, permutations = 1000,parallel=15)
  
  return(list(final_result = final_result, significant_results = significant_results))
}





## load data
## load data
library(readxl)

expression_data_metabolome <- read_excel("2_data/iPOP-project/3_data_analysis/plasma_metabolomics/data_preparation/metabolite/expression_data.xlsx")

sample_info_metabolome <- data.frame(read_excel("2_data/iPOP-project/3_data_analysis/plasma_metabolomics/data_preparation/metabolite/sample_info.xlsx"))

variable_info_metabolome <- data.frame(read_excel("2_data/iPOP-project/3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info.xlsx"))


row.names(expression_data_metabolome)<-variable_info_metabolome$variable_id
# data transform and 

sample_info_metabolome <- sample_info_metabolome %>% filter(!is.na(Gender))
sample_info_metabolome <- sample_info_metabolome %>% filter(!is.na(IRIS))
sample_info_metabolome <- sample_info_metabolome %>% filter(!is.na(BMI))


expression_data_metabolome<-expression_data_metabolome[,sample_info_metabolome$sample_id]


expression_data_metabolome<-log2(expression_data_metabolome)

expression_data_metabolome <-lm_adjust(expression_data = expression_data_metabolome,
                                       sample_info = sample_info_metabolome,
                                       threads = 3)

expression_data_metabolome<-data.frame(t(expression_data_metabolome))


colnames(expression_data_metabolome)<-variable_info_metabolome$variable_id





# calculate

gut_adonis<-MM_adonis(gut_microbiome_table,expression_data_metabolome)
oral_adonis<-MM_adonis(oral_microbiome_table,expression_data_metabolome)
nasal_adonis<-MM_adonis(nasal_microbiome_table,expression_data_metabolome)
skin_adonis<-MM_adonis(skin_microbiome_table,expression_data_metabolome)

metabolite_annotation<-read_excel("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx")

metabolite_annotation<-subset(metabolite_annotation,metabolite_annotation$HMDB.Source.Microbial=="TRUE")
expression_data_metabolome_microbial<-expression_data_metabolome[,metabolite_annotation$variable_id]

# 计算微生物来源的代谢物的解释度

gut_adonis_micro<-MM_adonis(gut_microbiome_table,expression_data_metabolome_microbial)
oral_adonis_micro<-MM_adonis(oral_microbiome_table,expression_data_metabolome_microbial)
nasal_adonis_micro<-MM_adonis(nasal_microbiome_table,expression_data_metabolome_microbial)
skin_adonis_micro<-MM_adonis(skin_microbiome_table,expression_data_metabolome_microbial)


# plot




adonis_r2<-c(gut_adonis$final_result$R2[1],oral_adonis$final_result$R2[1],nasal_adonis$final_result$R2[1],skin_adonis$final_result$R2[1])

name<-c("gut","oral","nasal","skin")

adonis_r2<-data.frame(adonis_r2,name)

colnames(adonis_r2)<-c("adonis_r2","bodysite")
adonis_r2$adonis_r2<-adonis_r2$adonis_r2*100
adonis_r2$adonis_r2<-format(adonis_r2$adonis_r2, digits =3, scientific = FALSE)






adonis_r2_micro<-c(0.2512,0.1521,0.0912,0.1322)

name<-c("gut","oral","nasal","skin")

adonis_r2_micro<-data.frame(adonis_r2_micro,name)

colnames(adonis_r2_micro)<-c("adonis_r2","bodysite")
adonis_r2_micro$adonis_r2<-adonis_r2_micro$adonis_r2*100
adonis_r2_micro$adonis_r2<-format(adonis_r2_micro$adonis_r2, digits =3, scientific = FALSE)


adonis_r2$group<-"Total"
adonis_r2_micro$group<-"Microbial"

adonis_r2_all<-rbind(adonis_r2_micro,adonis_r2)

adonis_r2_all$bodysite<-factor(adonis_r2_all$bodysite,levels = c("gut","oral","skin","nasal"))

library(ggpattern)
adonis_r2_all$adonis_r2<-as.numeric(adonis_r2_all$adonis_r2)
ggplot(adonis_r2_all, aes(x = bodysite, y = adonis_r2, fill = bodysite, pattern = group)) +
  geom_bar_pattern(stat = "identity", 
                   position = position_dodge(width = 0.7), 
                   width = 0.6,
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025) +
  geom_text(aes(label = paste(adonis_r2, "%", sep = "")), 
            position = position_dodge(width = 0.7),
            vjust = -0.5, size = 5) +
  scale_fill_manual(values = body_site_color) +
  scale_pattern_manual(values = c( "stripe","none")) +
  labs(title = "Estimated Variance by Different bodysite microbiome",
       x = NULL, 
       y = "Estimated variance (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(family = "Helvetica", size = 15),
        legend.position = "top")
