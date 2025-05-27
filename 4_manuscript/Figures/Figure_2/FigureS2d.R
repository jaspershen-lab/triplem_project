rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")

library(tidyverse)
library(tidymass)

###load("data)
load("3_data_analysis/oral_microbiome/data_preparation/object_cross_section")

oral_object <- object_cross_section

load(
  "3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section"
)

metabolomics_object <- object_cross_section

dir.create("3_data_analysis/oral_microbiome/spearman/cross_section/",
           recursive = TRUE)

setwd("3_data_analysis/oral_microbiome/spearman/cross_section/")

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


##### run pca for metabolomics and oral microbiome


metabolomics_pca_object <-
  metabolomics_temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  massstat::run_pca()

## 计算所有genus与代谢组数据PC1之间的相关性
metabolomics_pc <- t(metabolomics_temp_object@expression_data)
oral_pc <- t(oral_temp_object@expression_data)


share_index <- intersect(row.names(metabolomics_pc), row.names(oral_pc))

metabolomics_pc <- metabolomics_pc[share_index, ]
oral_pc <- oral_pc[share_index, ]


# 假设 'metabolomics_pc' 和 'oral_pc' 是您的数据框，并已正确设置
# 确保两个数据框的样本（行名）相同且顺序一致

# 将数据框转换为矩阵以便进行相关性计算
metabolomics_matrix <- as.matrix(metabolomics_pc)
oral_matrix <- as.matrix(oral_pc)

# 计算相关性矩阵
correlation_results <- apply(metabolomics_matrix, 2, function(metabolite) {
  apply(oral_matrix, 2, function(species) {
    cor.test(metabolite, species, method = "spearman")
  })
})

# 提取相关系数和p值
cor_values <- sapply(correlation_results, function(x)
  sapply(x, function(y)
    y$estimate))
p_values <- sapply(correlation_results, function(x)
  sapply(x, function(y)
    y$p.value))

# 将结果转换为数据框


cor_values[p_values > 0.05] <- 0

cor_values_filtered <- cor_values[rowSums(cor_values) != 0, ]

cor_values_filtered <- cor_values_filtered[, colSums(cor_values_filtered) != 0]



# 选取top100的代谢物
num_count <- apply(cor_values_filtered, 2, function(x)
  length(which(x != 0))) %>% sort(decreasing = TRUE) %>% data.frame()
metabolomics_top100 <- rownames(num_count)[1:100]



metabolomics_pca_object <-
  metabolomics_temp_object[metabolomics_top100, ] %>%
  activate_mass_dataset(what = "variable_info") %>%
  massstat::run_pca()


metabolomics_top100_PC1 <- metabolomics_pca_object$x[, 1:5]
metabolomics_top100_PC1 <- metabolomics_top100_PC1[share_index, ]

##cor(oral_pc,metabolomics_top100_PC1)%>%View()


metabolomics_matrix <- as.matrix(metabolomics_top100_PC1)
oral_matrix <- as.matrix(oral_pc)

# 计算相关性矩阵
correlation_results <- apply(metabolomics_matrix, 2, function(metabolite) {
  apply(oral_matrix, 2, function(species) {
    cor.test(metabolite, species, method = "spearman")
  })
})

# 提取相关系数和p值
cor_values <- sapply(correlation_results, function(x)
  sapply(x, function(y)
    y$estimate))
p_values <- sapply(correlation_results, function(x)
  sapply(x, function(y)
    y$p.value))

cor_values[p_values > 0.05] <- 0


# top20 genus correlation with metabolomics PC1


genus_top20_mete <- data.frame(cbind(cor_values[, 1], p_values[, 1]))

rownames(genus_top20_mete) <- rownames(p_values)

colnames(genus_top20_mete) <- c("rho", "p_val")

genus_top20_mete <- genus_top20_mete %>% rownames_to_column("row_name")
genus_top30_mete <- genus_top20_mete %>%
  dplyr::mutate(abs_rho = abs(rho)) %>% # 创建一个绝对值列
  dplyr::arrange(desc(abs_rho)) %>%       # 根据绝对值降序排序
  dplyr::slice_head(n = 20)           # 选择前20行

variable_info <- oral_temp_object@variable_info

genus_top30_mete <- merge(genus_top30_mete,
                          variable_info,
                          by.x = "row_name",
                          by.y = "variable_id")

genus_top30_mete <- genus_top30_mete %>%
  arrange(rho)

# 定义门水平的颜色方案（如果phylum_color未定义）
if (!exists("phylum_color")) {
  phylum_color <- setNames(
    c(
      "#E41A1C",
      "#377EB8",
      "#4DAF4A",
      "#984EA3",
      "#FF7F00",
      "#FFFF33",
      "#A65628",
      "#F781BF"
    ),
    unique(genus_top30_mete$Phylum)
  )
}

p<-ggplot(genus_top30_mete, aes(
  x = reorder(Genus, rho),
  y = rho,
  group = Genus
)) +
  geom_segment(aes(y = 0, yend = rho, xend = Genus),
               color = "#a1d5b9",
               size = 1.5) +  # 棒部分
  geom_point(aes(color = Phylum), size = 3) +  # 圆球部分，颜色根据phylum
  scale_color_manual(values = phylum_color) +  # 自定义phylum颜色
  labs(x = "Bacteria", y = "Rho Value") +    # 添加坐标轴标签
  theme_light(base_size = 14) +  # 使用light主题，基础字体大小为14
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 10,
      color = "grey30",
      face = "italic"
    ),
    axis.text.y = element_text(size = 12, color = "grey30"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )




ggsave(
  p,
  filename = file.path(r4projects::get_project_wd(), "4_manuscript/Figures/Figure_2/figure_s2d1.pdf"),
  width = 5,
  height = 10
)