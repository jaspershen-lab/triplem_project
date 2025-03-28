setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")

library(tidyverse)
library(tidymass)

###load("data)
load("3_data_analysis/gut_microbiome/data_preparation/object_cross_section")

gut_object<-object_cross_section

load("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section")

metabolomics_object<-object_cross_section

dir.create("3_data_analysis/gut_microbiome/spearman/cross_section/",recursive = TRUE)

setwd("3_data_analysis/gut_microbiome/spearman/cross_section/")


metabolite_annotation<-read_excel("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx")
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



# 
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



##### run pca for metabolomics and gut microbiome


metabolomics_pca_object <-
  metabolomics_temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  massstat::run_pca()

## 计算所有genus与代谢组数据PC1之间的相关性
metabolomics_pc<-t(metabolomics_temp_object@expression_data)
gut_pc<-t(gut_temp_object@expression_data)


share_index<-intersect(row.names(metabolomics_pc),row.names(gut_pc))



metabolomics_pc<-metabolomics_pc[share_index,]
gut_pc<-gut_pc[share_index,]


# 假设 'metabolomics_pc' 和 'gut_pc' 是您的数据框，并已正确设置
# 确保两个数据框的样本（行名）相同且顺序一致

# 将数据框转换为矩阵以便进行相关性计算
metabolomics_matrix <- as.matrix(metabolomics_pc)
gut_matrix <- as.matrix(gut_pc)

# 计算相关性矩阵
correlation_results <- apply(metabolomics_matrix, 2, function(metabolite) {
  apply(gut_matrix, 2, function(species) {
    cor.test(metabolite, species, method = "spearman")
  })
})

# 提取相关系数和p值
cor_values <- sapply(correlation_results, function(x) sapply(x, function(y) y$estimate))
p_values <- sapply(correlation_results, function(x) sapply(x, function(y) y$p.value))

# 将结果转换为数据框


cor_values[p_values > 0.05] <- 0

cor_values_filtered <- cor_values[rowSums(cor_values) != 0, ]

cor_values_filtered <- cor_values_filtered[, colSums(cor_values_filtered) != 0]



# 选取top100的代谢物
num_count<-apply(cor_values_filtered, 2, function(x) length(which(x != 0)))%>%sort(decreasing = TRUE)%>%data.frame()
metabolomics_top100<-rownames(num_count)[1:100]



metabolomics_pca_object <-
  metabolomics_temp_object[metabolomics_top100,] %>%
  activate_mass_dataset(what = "variable_info") %>%
  massstat::run_pca()


metabolomics_top100_PC1<-metabolomics_pca_object$x[,1:5]
metabolomics_top100_PC1<-metabolomics_top100_PC1[share_index,]

##cor(gut_pc,metabolomics_top100_PC1)%>%View()


metabolomics_matrix <- as.matrix(metabolomics_top100_PC1)
gut_matrix <- as.matrix(gut_pc)

# 计算相关性矩阵
correlation_results <- apply(metabolomics_matrix, 2, function(metabolite) {
  apply(gut_matrix, 2, function(species) {
    cor.test(metabolite, species, method = "spearman")
  })
})

# 提取相关系数和p值
cor_values <- sapply(correlation_results, function(x) sapply(x, function(y) y$estimate))
p_values <- sapply(correlation_results, function(x) sapply(x, function(y) y$p.value))

cor_values[p_values > 0.05] <- 0


# top20 genus correlation with metabolomics PC1


genus_top20_mete<-data.frame(cbind(cor_values[,1],p_values[,1]))

rownames(genus_top20_mete)<-rownames(p_values)

colnames(genus_top20_mete)<-c("rho","p_val")

genus_top20_mete <- genus_top20_mete %>% rownames_to_column("row_name")
genus_top30_mete <- genus_top20_mete %>%
  mutate(abs_rho = abs(rho)) %>% # 创建一个绝对值列
  arrange(desc(abs_rho))

genus_top30_mete<-merge(genus_top30_mete,gut_temp_object@variable_info,by.x="row_name",by.y="variable_id")

gut_cor_PC1<-genus_top30_mete

##################


###load("data)

setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")

library(tidyverse)
library(tidymass)
load("3_data_analysis/oral_microbiome/data_preparation/object_cross_section")

gut_object<-object_cross_section

load("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section")

metabolomics_object<-object_cross_section

dir.create("3_data_analysis/oral_microbiome/spearman/cross_section/",recursive = TRUE)

setwd("3_data_analysis/oral_microbiome/spearman/cross_section/")



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



# 
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



##### run pca for metabolomics and gut microbiome


metabolomics_pca_object <-
  metabolomics_temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  massstat::run_pca()

## 计算所有genus与代谢组数据PC1之间的相关性
metabolomics_pc<-t(metabolomics_temp_object@expression_data)
gut_pc<-t(gut_temp_object@expression_data)


share_index<-intersect(row.names(metabolomics_pc),row.names(gut_pc))



metabolomics_pc<-metabolomics_pc[share_index,]
gut_pc<-gut_pc[share_index,]


# 假设 'metabolomics_pc' 和 'gut_pc' 是您的数据框，并已正确设置
# 确保两个数据框的样本（行名）相同且顺序一致

# 将数据框转换为矩阵以便进行相关性计算
metabolomics_matrix <- as.matrix(metabolomics_pc)
gut_matrix <- as.matrix(gut_pc)

# 计算相关性矩阵
correlation_results <- apply(metabolomics_matrix, 2, function(metabolite) {
  apply(gut_matrix, 2, function(species) {
    cor.test(metabolite, species, method = "spearman")
  })
})

# 提取相关系数和p值
cor_values <- sapply(correlation_results, function(x) sapply(x, function(y) y$estimate))
p_values <- sapply(correlation_results, function(x) sapply(x, function(y) y$p.value))

# 将结果转换为数据框


cor_values[p_values > 0.05] <- 0

cor_values_filtered <- cor_values[rowSums(cor_values) != 0, ]

cor_values_filtered <- cor_values_filtered[, colSums(cor_values_filtered) != 0]



# 选取top100的代谢物
num_count<-apply(cor_values_filtered, 2, function(x) length(which(x != 0)))%>%sort(decreasing = TRUE)%>%data.frame()
metabolomics_top100<-rownames(num_count)[1:100]



metabolomics_pca_object <-
  metabolomics_temp_object[metabolomics_top100,] %>%
  activate_mass_dataset(what = "variable_info") %>%
  massstat::run_pca()


metabolomics_top100_PC1<-metabolomics_pca_object$x[,1:5]
metabolomics_top100_PC1<-metabolomics_top100_PC1[share_index,]

##cor(gut_pc,metabolomics_top100_PC1)%>%View()


metabolomics_matrix <- as.matrix(metabolomics_top100_PC1)
gut_matrix <- as.matrix(gut_pc)

# 计算相关性矩阵
correlation_results <- apply(metabolomics_matrix, 2, function(metabolite) {
  apply(gut_matrix, 2, function(species) {
    cor.test(metabolite, species, method = "spearman")
  })
})

# 提取相关系数和p值
cor_values <- sapply(correlation_results, function(x) sapply(x, function(y) y$estimate))
p_values <- sapply(correlation_results, function(x) sapply(x, function(y) y$p.value))

cor_values[p_values > 0.05] <- 0


# top20 genus correlation with metabolomics PC1


genus_top20_mete<-data.frame(cbind(cor_values[,1],p_values[,1]))

rownames(genus_top20_mete)<-rownames(p_values)

colnames(genus_top20_mete)<-c("rho","p_val")

genus_top20_mete <- genus_top20_mete %>% rownames_to_column("row_name")
genus_top30_mete <- genus_top20_mete %>%
  mutate(abs_rho = abs(rho)) %>% # 创建一个绝对值列
  arrange(desc(abs_rho))

genus_top30_mete<-merge(genus_top30_mete,gut_temp_object@variable_info,by.x="row_name",by.y="variable_id")


oral_cor_PC1<-genus_top30_mete



##################

setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")

library(tidyverse)
library(tidymass)
load("3_data_analysis/skin_microbiome/data_preparation/object_cross_section")

gut_object<-object_cross_section

load("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section")

metabolomics_object<-object_cross_section

dir.create("3_data_analysis/skin_microbiome/spearman/cross_section/",recursive = TRUE)

setwd("3_data_analysis/skin_microbiome/spearman/cross_section/")



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



# 
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



##### run pca for metabolomics and gut microbiome


metabolomics_pca_object <-
  metabolomics_temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  massstat::run_pca()

## 计算所有genus与代谢组数据PC1之间的相关性
metabolomics_pc<-t(metabolomics_temp_object@expression_data)
gut_pc<-t(gut_temp_object@expression_data)


share_index<-intersect(row.names(metabolomics_pc),row.names(gut_pc))



metabolomics_pc<-metabolomics_pc[share_index,]
gut_pc<-gut_pc[share_index,]


# 假设 'metabolomics_pc' 和 'gut_pc' 是您的数据框，并已正确设置
# 确保两个数据框的样本（行名）相同且顺序一致

# 将数据框转换为矩阵以便进行相关性计算
metabolomics_matrix <- as.matrix(metabolomics_pc)
gut_matrix <- as.matrix(gut_pc)

# 计算相关性矩阵
correlation_results <- apply(metabolomics_matrix, 2, function(metabolite) {
  apply(gut_matrix, 2, function(species) {
    cor.test(metabolite, species, method = "spearman")
  })
})

# 提取相关系数和p值
cor_values <- sapply(correlation_results, function(x) sapply(x, function(y) y$estimate))
p_values <- sapply(correlation_results, function(x) sapply(x, function(y) y$p.value))

# 将结果转换为数据框


cor_values[p_values > 0.05] <- 0

cor_values_filtered <- cor_values[rowSums(cor_values) != 0, ]

cor_values_filtered <- cor_values_filtered[, colSums(cor_values_filtered) != 0]



# 选取top100的代谢物
num_count<-apply(cor_values_filtered, 2, function(x) length(which(x != 0)))%>%sort(decreasing = TRUE)%>%data.frame()
metabolomics_top100<-rownames(num_count)[1:100]



metabolomics_pca_object <-
  metabolomics_temp_object[metabolomics_top100,] %>%
  activate_mass_dataset(what = "variable_info") %>%
  massstat::run_pca()


metabolomics_top100_PC1<-metabolomics_pca_object$x[,1:5]
metabolomics_top100_PC1<-metabolomics_top100_PC1[share_index,]

##cor(gut_pc,metabolomics_top100_PC1)%>%View()


metabolomics_matrix <- as.matrix(metabolomics_top100_PC1)
gut_matrix <- as.matrix(gut_pc)

# 计算相关性矩阵
correlation_results <- apply(metabolomics_matrix, 2, function(metabolite) {
  apply(gut_matrix, 2, function(species) {
    cor.test(metabolite, species, method = "spearman")
  })
})

# 提取相关系数和p值
cor_values <- sapply(correlation_results, function(x) sapply(x, function(y) y$estimate))
p_values <- sapply(correlation_results, function(x) sapply(x, function(y) y$p.value))

cor_values[p_values > 0.05] <- 0


# top20 genus correlation with metabolomics PC1


genus_top20_mete<-data.frame(cbind(cor_values[,1],p_values[,1]))

rownames(genus_top20_mete)<-rownames(p_values)

colnames(genus_top20_mete)<-c("rho","p_val")

genus_top20_mete <- genus_top20_mete %>% rownames_to_column("row_name")
genus_top30_mete <- genus_top20_mete %>%
  mutate(abs_rho = abs(rho)) %>% # 创建一个绝对值列
  arrange(desc(abs_rho))

genus_top30_mete<-merge(genus_top30_mete,gut_temp_object@variable_info,by.x="row_name",by.y="variable_id")


skin_cor_PC1<-genus_top30_mete


###########


setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")

library(tidyverse)
library(tidymass)
load("3_data_analysis/nasal_microbiome/data_preparation/object_cross_section")

gut_object<-object_cross_section

load("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section")

metabolomics_object<-object_cross_section

dir.create("3_data_analysis/nasal_microbiome/spearman/cross_section/",recursive = TRUE)

setwd("3_data_analysis/nasal_microbiome/spearman/cross_section/")



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



# 
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



##### run pca for metabolomics and gut microbiome


metabolomics_pca_object <-
  metabolomics_temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  massstat::run_pca()

## 计算所有genus与代谢组数据PC1之间的相关性
metabolomics_pc<-t(metabolomics_temp_object@expression_data)
gut_pc<-t(gut_temp_object@expression_data)


share_index<-intersect(row.names(metabolomics_pc),row.names(gut_pc))



metabolomics_pc<-metabolomics_pc[share_index,]
gut_pc<-gut_pc[share_index,]


# 假设 'metabolomics_pc' 和 'gut_pc' 是您的数据框，并已正确设置
# 确保两个数据框的样本（行名）相同且顺序一致

# 将数据框转换为矩阵以便进行相关性计算
metabolomics_matrix <- as.matrix(metabolomics_pc)
gut_matrix <- as.matrix(gut_pc)

# 计算相关性矩阵
correlation_results <- apply(metabolomics_matrix, 2, function(metabolite) {
  apply(gut_matrix, 2, function(species) {
    cor.test(metabolite, species, method = "spearman")
  })
})

# 提取相关系数和p值
cor_values <- sapply(correlation_results, function(x) sapply(x, function(y) y$estimate))
p_values <- sapply(correlation_results, function(x) sapply(x, function(y) y$p.value))

# 将结果转换为数据框


cor_values[p_values > 0.05] <- 0

cor_values_filtered <- cor_values[rowSums(cor_values) != 0, ]

cor_values_filtered <- cor_values_filtered[, colSums(cor_values_filtered) != 0]



# 选取top100的代谢物
num_count<-apply(cor_values_filtered, 2, function(x) length(which(x != 0)))%>%sort(decreasing = TRUE)%>%data.frame()
metabolomics_top100<-rownames(num_count)[1:100]



metabolomics_pca_object <-
  metabolomics_temp_object[metabolomics_top100,] %>%
  activate_mass_dataset(what = "variable_info") %>%
  massstat::run_pca()


metabolomics_top100_PC1<-metabolomics_pca_object$x[,1:5]
metabolomics_top100_PC1<-metabolomics_top100_PC1[share_index,]

##cor(gut_pc,metabolomics_top100_PC1)%>%View()


metabolomics_matrix <- as.matrix(metabolomics_top100_PC1)
gut_matrix <- as.matrix(gut_pc)

# 计算相关性矩阵
correlation_results <- apply(metabolomics_matrix, 2, function(metabolite) {
  apply(gut_matrix, 2, function(species) {
    cor.test(metabolite, species, method = "spearman")
  })
})

# 提取相关系数和p值
cor_values <- sapply(correlation_results, function(x) sapply(x, function(y) y$estimate))
p_values <- sapply(correlation_results, function(x) sapply(x, function(y) y$p.value))

cor_values[p_values > 0.05] <- 0


# top20 genus correlation with metabolomics PC1


genus_top20_mete<-data.frame(cbind(cor_values[,1],p_values[,1]))

rownames(genus_top20_mete)<-rownames(p_values)

colnames(genus_top20_mete)<-c("rho","p_val")

genus_top20_mete <- genus_top20_mete %>% rownames_to_column("row_name")
genus_top30_mete <- genus_top20_mete %>%
  mutate(abs_rho = abs(rho)) %>% # 创建一个绝对值列
  arrange(desc(abs_rho))

genus_top30_mete<-merge(genus_top30_mete,gut_temp_object@variable_info,by.x="row_name",by.y="variable_id")


nasal_cor_PC1<-genus_top30_mete





##########
all_genus_list<-list(c(gut_cor_PC1$Genus),c(oral_cor_PC1$Genus),c(skin_cor_PC1$Genus),c(nasal_cor_PC1$Genus))

all_genus<-Reduce(union, all_genus_list)



all_genus<-data.frame(all_genus)


all_genus<-merge(all_genus,gut_cor_PC1[,c("abs_rho","Genus")],by.x="all_genus",by.y="Genus",all.x=TRUE)
all_genus<-merge(all_genus,oral_cor_PC1[,c("abs_rho","Genus")],by.x="all_genus",by.y="Genus",all.x=TRUE)
all_genus<-merge(all_genus,skin_cor_PC1[,c("abs_rho","Genus")],by.x="all_genus",by.y="Genus",all.x=TRUE)
all_genus<-merge(all_genus,nasal_cor_PC1[,c("abs_rho","Genus")],by.x="all_genus",by.y="Genus",all.x=TRUE)

colnames(all_genus)<-c("Genus","gut","oral","skin","nasal")


# 筛选全部为0的物种


all_genus[is.na(all_genus)] <- 0
all_genus <- all_genus[rowSums(all_genus[,2:5] != 0) > 0, ]

all_genus_tax<-rbind(gut_cor_PC1[,5:10],oral_cor_PC1[,5:10],skin_cor_PC1[,5:10],nasal_cor_PC1[,5:10])
all_genus_tax<-unique(all_genus_tax)

all_genus<-merge(all_genus,all_genus_tax,by="Genus")

all_genus<-all_genus[-128,]