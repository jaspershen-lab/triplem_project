rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")


## read MetaCardis data
library(readxl)

FBIP_metagenomic <- read.table("2_data/FBIP-main/otu_table/merged_abundance_table_species.txt",
                               header = TRUE)

FBIP_metagenomic <- aggregate(FBIP_metagenomic[, 8:669], list(FBIP_metagenomic$Genus), sum)

FBIP_metagenomic$Group.1 <- gsub("g__", "", FBIP_metagenomic$Group.1)


row.names(FBIP_metagenomic) <- FBIP_metagenomic$Group.1

FBIP_metagenomic <- FBIP_metagenomic[, -1]

FBIP_metagenomic <- data.frame(t(FBIP_metagenomic), check.names = FALSE)



FBIP_metabolome <- read.table(
  "2_data/FBIP-main/metabolites data/plasma metabolites data.txt",
  header = TRUE,
  sep = "\t"
)

FBIP_metadata <- read.table("2_data/FBIP-main/metadata/FBIP_metadata.txt",
                            header = TRUE,
                            sep = "\t")
FBIP_metadata <- subset(FBIP_metadata, time_num %in% c("W0", "W4", "W16"))

FBIP_metagenomic <- FBIP_metagenomic[FBIP_metadata$SampleID, ]


FBIP_metabolome <- subset(FBIP_metabolome, Time %in% c("W0", "W4", "W16"))

FBIP_metabolome <- merge(FBIP_metabolome,
                         FBIP_metadata[, c("Time", "SampleID")],
                         by.x = "ST",
                         by.y = "Time")
rownames(FBIP_metabolome) <- FBIP_metabolome$SampleID
FBIP_metabolome <- FBIP_metabolome[, c(-1, -2, -3, -197)]

FBIP_metagenomic <- FBIP_metagenomic[rownames(FBIP_metabolome), ]

FBIP_metabolome_annotation <- read.table(
  "2_data/FBIP-main/metabolites data/plasma_metabolites_name_group.txt",
  header = TRUE,
  sep = "\t"
)

# 过滤没有HMDB ID的代谢物

FBIP_metabolome_annotation <- subset(FBIP_metabolome_annotation, HMDB !=
                                       c(""))

FBIP_metabolome <- FBIP_metabolome[, FBIP_metabolome_annotation$ID]

colnames(FBIP_metabolome) <- FBIP_metabolome_annotation$HMDB


# 读取iPOP数据


library(tidyverse)
library(tidymass)
library(plyr)
library(microbiomedataset)
###load("data)
load("3_data_analysis/gut_microbiome/data_preparation/object_cross_section")

gut_object <- object_cross_section

load(
  "3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section"
)

metabolomics_object <- object_cross_section

metabolite_annotation <- read_excel(
  "3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx"
)




####only remain the genus level


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
  which(non_zero_per > 0.2)

gut_object <-
  gut_object[idx, ]


gut_object <-
  gut_object %>%
  transform2relative_intensity()

gut_expression_data <-
  extract_expression_data(gut_object) %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()



gut_sample_info <-
  gut_object@sample_info


gut_temp_object <- gut_object
gut_temp_object@expression_data <- gut_expression_data


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



##筛选共同的属
com_tax <- intersect(gut_temp_object@variable_info$Genus,
                     colnames(FBIP_metagenomic))


microbiome_data_ipop <- gut_temp_object@expression_data
microbiome_tax_ipop <- gut_temp_object@variable_info

rownames(microbiome_data_ipop) <- microbiome_tax_ipop$Genus

microbiome_data_ipop <- microbiome_data_ipop[com_tax, ]

FBIP_metagenomic <- data.frame(t(FBIP_metagenomic), check.names = FALSE)

FBIP_metagenomic <- FBIP_metagenomic[com_tax, ]

FBIP_metagenomic <-
  FBIP_metagenomic %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()


### 挑选共同的代谢物

metabolome_data_ipop <- metabolomics_temp_object@expression_data
metabolome_ann_ipop <- data.frame(metabolite_annotation)

com_metabolites <- intersect(metabolome_ann_ipop$HMDB, FBIP_metabolome_annotation$HMDB)

metabolome_ann_ipop <- subset(metabolome_ann_ipop, HMDB %in% com_metabolites)

metabolome_data_ipop <- metabolome_data_ipop[metabolome_ann_ipop$variable_id, ]


FBIP_metabolome_annotation <- subset(FBIP_metabolome_annotation, HMDB %in%
                                       com_metabolites)

FBIP_metabolome <- FBIP_metabolome[, com_metabolites]

common_samples <- intersect(colnames(metabolome_data_ipop),
                            colnames(microbiome_data_ipop))
metabolome_data <- metabolome_data_ipop[, common_samples]
microbiome_data <- microbiome_data_ipop[, common_samples]

# 转置数据矩阵以便进行相关性计算
metabolome_t <- t(metabolome_data)
microbiome_t <- t(microbiome_data)

# 初始化一个数据框来存储所有相关性结果
correlation_results <- data.frame()

# 计算每个代谢物与每个细菌之间的相关性
for (i in 1:ncol(metabolome_t)) {
  metabolite_name <- colnames(metabolome_t)[i]
  metabolite_data <- metabolome_t[, i]
  
  for (j in 1:ncol(microbiome_t)) {
    microbe_name <- colnames(microbiome_t)[j]
    microbe_data <- microbiome_t[, j]
    
    # 使用complete.cases移除任何含有NA值的样本
    valid_indices <- complete.cases(metabolite_data, microbe_data)
    
    if (sum(valid_indices) > 5) {
      # 确保至少有足够的有效样本
      # 计算Spearman相关系数和p值
      cor_test <- cor.test(metabolite_data[valid_indices], microbe_data[valid_indices], method = "spearman")
      
      # 存储结果
      result <- data.frame(
        Metabolite = metabolite_name,
        Microbe = microbe_name,
        Correlation = cor_test$estimate,
        P_value = cor_test$p.value
      )
      
      correlation_results <- rbind(correlation_results, result)
    }
  }
}

# 计算校正后的p值（FDR校正）
correlation_results$FDR <- p.adjust(correlation_results$P_value, method = "BH")

# 计算相关性的绝对值，用于排序
correlation_results$Abs_Correlation <- abs(correlation_results$Correlation)

# 按相关性绝对值降序排序
correlation_results <- correlation_results[order(correlation_results$Abs_Correlation, decreasing = TRUE), ]

# 选择相关性最强的前100个
top_100_correlations <- head(correlation_results, 500)


top_100_correlations <- merge(top_100_correlations,
                              metabolite_annotation[, c("variable_id", "HMDB", "HMDB.Name")],
                              by.x = "Metabolite",
                              by.y = "variable_id")

top_100_correlations$bac_meta <- paste(top_100_correlations$Microbe,
                                       top_100_correlations$HMDB,
                                       sep = "_")


top_100_correlations <- top_100_correlations %>%
  group_by(bac_meta) %>%
  dplyr::summarise(across(where(is.numeric), mean, na.rm = TRUE))






library(tidyverse)  # 数据处理
library(Hmisc)

#######FBIP
# 读取数据
# 假设您的数据存储在CSV文件中，请根据实际情况修改文件路径
metagenomic_data <- FBIP_metagenomic
metabolome_data <- data.frame(t(FBIP_metabolome))

# 检查数据结构
cat("代谢物数据维度:", dim(metagenomic_data), "\n")
cat("细菌数据维度:", dim(metabolome_data), "\n")

# 确保样本列名一致并按相同顺序排列
# 提取公共样本
common_samples <- intersect(colnames(metagenomic_data), colnames(metabolome_data))

# 如果没有共同样本，则停止程序
if (length(common_samples) == 0) {
  stop("没有共同的样本名在两个数据集中")
}

# 使用公共样本筛选数据集
metagenomic_filtered <- metagenomic_data[, common_samples]
metabolome_filtered <- metabolome_data[, common_samples]

# 相关性分析
# 转置使样本在行，特征在列
metagenomic_t <- t(metagenomic_filtered)
metabolome_t <- t(metabolome_filtered)

# 计算相关性矩阵（Spearman相关系数）
correlation_result <- rcorr(as.matrix(metagenomic_t), as.matrix(metabolome_t), type = "spearman")

# 提取相关系数和P值
cor_coef <- correlation_result$r
cor_pval <- correlation_result$P

# 提取代谢物和细菌之间的相关系数部分
# 假设metagenomic_t包含代谢物，metabolome_t包含细菌
cor_subset <- cor_coef[1:ncol(metagenomic_t), (ncol(metagenomic_t) + 1):(ncol(metagenomic_t) +
                                                                           ncol(metabolome_t))]
pval_subset <- cor_pval[1:ncol(metagenomic_t), (ncol(metagenomic_t) + 1):(ncol(metagenomic_t) +
                                                                            ncol(metabolome_t))]

# 将相关矩阵转换为表格形式
correlation_table <- data.frame()

for (i in 1:nrow(cor_subset)) {
  for (j in 1:ncol(cor_subset)) {
    metabolite <- rownames(cor_subset)[i]
    bacteria <- colnames(cor_subset)[j]
    rho <- cor_subset[i, j]
    pvalue <- pval_subset[i, j]
    
    correlation_table <- rbind(
      correlation_table,
      data.frame(
        Metabolite = metabolite,
        Bacteria = bacteria,
        Rho = rho,
        Pvalue = pvalue
      )
    )
  }
}

# 按p值排序
correlation_table <- correlation_table[order(correlation_table$Pvalue), ]

# 添加FDR校正的p值
correlation_table$AdjustedPvalue <- p.adjust(correlation_table$Pvalue, method = "BH")

# 添加显著性标记
correlation_table$Significance <- ""
correlation_table$Significance[correlation_table$Pvalue < 0.05] <- "*"
correlation_table$Significance[correlation_table$Pvalue < 0.01] <- "**"
correlation_table$Significance[correlation_table$Pvalue < 0.001] <- "***"


correlation_table_FBIP <- correlation_table

correlation_table_FBIP$bac_meta <- paste(correlation_table_FBIP$Metabolite,
                                         correlation_table_FBIP$Bacteria,
                                         sep = "_")

correlation_table_FBIP <- subset(correlation_table_FBIP,
                                 bac_meta %in% top_100_correlations$bac_meta)



correlation_table_FBIP_ipop <- merge(correlation_table_FBIP, top_100_correlations, by =
                                       "bac_meta")
correlation_table_FBIP_ipop <- correlation_table_FBIP_ipop[which(
  sign(correlation_table_FBIP_ipop$Rho) * sign(correlation_table_FBIP_ipop$Correlation) > 0
), ]

correlation_table_FBIP_ipop <- subset(correlation_table_FBIP_ipop, abs(Rho) >
                                        0.1)


plot <-
  ggplot(
    correlation_table_FBIP_ipop,
    aes(x = correlation_table_FBIP_ipop$Rho, y = correlation_table_FBIP_ipop$Correlation)
  ) +
  geom_point(
    shape = 21,
    size = 4,
    fill = "#A1D0C7",
    color = "white"
  ) +
  geom_smooth(method = "lm", colour = "grey50") + theme_light() + stat_cor(method = "spearman") +
  theme(
    legend.position = "none",
    #不需要图例
    axis.text.x =
      element_text(colour = "black", size = 14),
    #设置x轴刻度标签的字体属性
    axis.text.y =
      element_text(size = 14, face = "plain"),
    #设置x轴刻度标签的字体属性
    axis.title.y =
      element_text(size = 14, face = "plain"),
    #设置y轴的标题的字体属性
    axis.title.x =
      element_text(size = 14, face = "plain"),
    #设置x轴的标题的字体属性
    plot.title = element_text(size =
                                15, face = "bold", hjust = 0.5)
  ) + xlab("iPOP_Rho") + ylab("FBIP_Rho")

ggsave(
  plot,
  filename = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_2/figure_2d.pdf"
  ),
  width = 7,
  height = 7
)
