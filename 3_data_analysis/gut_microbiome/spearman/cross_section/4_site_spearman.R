rm(list = ls())
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


rownames(cor_values_filtered)<-rownames(gut_expression_data)



## 绘制桑基图
library(reshape2)
sangkey_data<-melt(cor_values_filtered)


sangkey_data<-subset(sangkey_data,!(value==0))

variable_info<-gut_temp_object@variable_info

sangkey_data<-merge(sangkey_data,variable_info,by.x="Var1",by.y="variable_id")

variable_info<-read_excel("/Users/zhangjingxiang/project/multi_microbiome_metabolome/2_data/iPOP-project/variable_info_metabolome.xlsx")

sangkey_data<-merge(sangkey_data,variable_info[,c("variable_id","HMDB.ID")],by.x="Var2",by.y="variable_id")

variable_info<-read_excel("/Users/zhangjingxiang/project/multi_microbiome_metabolome/2_data/iPOP-project/meta_class.xlsx")

variable_info<-subset(variable_info,!(variable_info$HMDB=="NA"))


variable_info<-unique(variable_info)
sangkey_data<-merge(sangkey_data,variable_info,by.x="HMDB.ID",by.y="HMDB")

sangkey_data$Site<-"gut"
sangkey_data_gut<-sangkey_data









sangkey_data<-rbind(sangkey_data_gut,sangkey_data_oral,sangkey_data_skin,sangkey_data_nasal)



plot_data<-sangkey_data_gut%>%make_long(Phylum,HMDB.Super.Class)


ggplot(plot_data,aes(x=x,next_x=next_x,
                     node=node,next_node=next_node,
                     fill=factor(node),label=node))+
  geom_sankey(flow.alpha=.6,
              node.color="gray30")+
  geom_sankey_label(size=3,color="white",fill="gray40")+
  scale_fill_manual(values=c(phylum_color,body_site_color))+
  theme_sankey(base_size=18)+
  labs(x=NULL)+
  theme(legend.position="none",
        plot.title=element_text(hjust=.5))






plot_data<-data.frame(table(sangkey_data$Site,sangkey_data$Genus))

plot_data<-plot_data %>%arrange(desc(Freq)) %>% slice_head(n = 50)

plot_data$site_genus <- paste(plot_data$Var1, plot_data$Var2, sep = "_")




plot_data_sorted <- plot_data[order(plot_data$Freq), ]

# 更新因子水平以反映新的排序
plot_data_sorted$site_genus <- factor(plot_data_sorted$site_genus , levels = unique(plot_data_sorted$site_genus ))


plot_data <- plot_data %>%
  arrange(site_genus, desc(Freq))

ggplot(plot_data, aes(x = site_genus, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_text(aes(label = paste(Freq)), vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("#DF8F44FF","#79AF97FF","#B24745FF","#00A1D5FF" )) +
  labs(title = "Estimated Variance by Different bodysite microbiome",
       x = NULL, 
       y = "Freq") +theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,family = "Helvetica", size = 6),
        legend.position = "none")