setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")

library(tidyverse)
library(tidymass)
library(readxl)

# 读取数据
## 肠道数据
gut_microbiome<-data.frame(read_excel("2_data/2023_Eransegal_NC/gut_table.xlsx"))
rownames(gut_microbiome)<-gut_microbiome$Sample_id
gut_microbiome<-gut_microbiome[,-1]
gut_metadata<-data.frame(read_excel("2_data/2023_Eransegal_NC/gut_metadata.xlsx"))
gut_tax<-data.frame(read_excel("2_data/2023_Eransegal_NC/gut_tax.xlsx"))
## 口腔数据
oral_microbiome<-data.frame(read_excel("2_data/2023_Eransegal_NC/oral_table.xlsx"))
rownames(oral_microbiome)<-oral_microbiome$Sample_id
oral_microbiome<-oral_microbiome[,-1]
oral_metadata<-data.frame(read_excel("2_data/2023_Eransegal_NC/oral_metadata.xlsx"))
oral_tax<-data.frame(read_excel("2_data/2023_Eransegal_NC/oral_tax.xlsx"))

## 代谢组数据
metabolome<-data.frame(read_excel("2_data/2023_Eransegal_NC/metabolome_table.xlsx"),check.names = FALSE)
rownames(metabolome)<-metabolome$Samole_id
metabolome<-metabolome[,-1]
metabolome_metadata<-data.frame(read_excel("2_data/2023_Eransegal_NC/metabolome_metadata.xlsx"))
metabolome_annotation<-data.frame(read_excel("2_data/2023_Eransegal_NC/metabolome_annotation.xlsx"))


core_metabolite<-c("p-cresol sulfate","p-cresol glucuronide*","2-hydroxyphenylacetate","4-hydroxyphenylacetate","phenylacetylglutamine","chenodeoxycholate")

# 双向中介效应分析：口腔菌群、肠道菌群与代谢物的关系
# 安装并加载必要的包
if (!requireNamespace("mediation", quietly = TRUE)) install.packages("mediation")
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")

library(mediation)
library(tidyverse)
library(readr)


gut_microbiome<-data.frame(t(gut_microbiome))
oral_microbiome<-data.frame(t(oral_microbiome))
metabolome<-data.frame(t(metabolome),check.names = FALSE)
# 确保样本名一致
common_samples <- Reduce(intersect, list(colnames(gut_microbiome), colnames(oral_microbiome), colnames(metabolome)))
gut_microbiome <- gut_microbiome[, common_samples]
oral_microbiome <- oral_microbiome[, common_samples]
metabolome <- metabolome[, common_samples]

# 转置数据以便于分析（行为样本，列为特征）
gut_data_t <- t(gut_microbiome)
oral_data_t <- t(oral_microbiome)
metabolome_data_t <- t(metabolome[core_metabolite,])

calculate_correlations <- function(data1, data2) {
  result <- matrix(NA, nrow = ncol(data1), ncol = ncol(data2))
  rownames(result) <- colnames(data1)
  colnames(result) <- colnames(data2)
  p_values <- result
  
  for (i in 1:ncol(data1)) {
    for (j in 1:ncol(data2)) {
      # 使用tryCatch处理可能的错误情况
      tryCatch({
        # 只使用两个向量中都有有效值的观测
        complete_obs <- complete.cases(data1[, i], data2[, j])
        if(sum(complete_obs) > 3) {  # 确保至少有足够的观测值进行相关性分析
          cor_test <- cor.test(data1[complete_obs, i], data2[complete_obs, j], method = "spearman")
          result[i, j] <- cor_test$estimate
          p_values[i, j] <- cor_test$p.value
        } else {
          result[i, j] <- NA
          p_values[i, j] <- NA
        }
      }, error = function(e) {
        result[i, j] <- NA
        p_values[i, j] <- NA
      })
    }
  }
  
  # FDR校正 (只对非NA值进行校正)
  p_val_vector <- as.vector(p_values)
  non_na_idx <- !is.na(p_val_vector)
  if(sum(non_na_idx) > 0) {  # 确保有非NA值需要校正
    adjusted_p <- p.adjust(p_val_vector[non_na_idx], method = "BH")
    p_val_vector[non_na_idx] <- adjusted_p
  }
  fdr_values <- matrix(p_val_vector, nrow = nrow(p_values))
  rownames(fdr_values) <- rownames(p_values)
  colnames(fdr_values) <- colnames(p_values)
  
  return(list(cor = result, p = p_values, fdr = fdr_values))
}

# 第一步：找出与代谢物相关联的微生物特征（肠道菌群和口腔菌群）
gut_metabolome_cor <- calculate_correlations(gut_data_t, metabolome_data_t)
gut_metabolome_associations <- which(gut_metabolome_cor$p< 0.05, arr.ind = TRUE)

oral_metabolome_cor <- calculate_correlations(oral_data_t, metabolome_data_t)
oral_metabolome_associations <- which(oral_metabolome_cor$p < 0.05, arr.ind = TRUE)

# 第二步：找出口腔菌群与肠道菌群之间的关联
oral_gut_cor <- calculate_correlations(oral_data_t, gut_data_t)
oral_gut_associations <- which(oral_gut_cor$p < 0.05, arr.ind = TRUE)

bidirectional_mediation_results <- data.frame(
  direction = character(),
  oral_feature = character(),
  gut_feature = character(),
  metabolite = character(),
  ACME = numeric(),
  ACME_p = numeric(),
  prop_mediated = numeric(),
  stringsAsFactors = FALSE
)
######

for (i in 1:nrow(oral_gut_associations)) {
  oral_idx <- oral_gut_associations[i, "row"]
  gut_idx <- oral_gut_associations[i, "col"]
  
  oral_feature <- rownames(oral_gut_cor$cor)[oral_idx]
  gut_feature <- colnames(oral_gut_cor$cor)[gut_idx]
  
  # 检查该肠道菌群是否与代谢物相关联
  gut_metabolite_associations <- which(gut_metabolome_cor$p[gut_idx, ] < 0.05, arr.ind = TRUE)
  
  if (length(gut_metabolite_associations) > 0) {
    for (j in 1:length(gut_metabolite_associations)) {
      met_idx <- gut_metabolite_associations[j]
      metabolite <- colnames(gut_metabolome_cor$cor)[met_idx]
      
      # 构建数据框用于中介分析
      med_data <- data.frame(
        oral = oral_data_t[, oral_feature],
        gut = gut_data_t[, gut_feature],
        metabolite = metabolome_data_t[, metabolite]
      )
      
      # 在进行分析之前移除任何含有NA的行
      med_data_complete <- na.omit(med_data)
      
      # 确保有足够的观测进行分析
      if (nrow(med_data_complete) >= 10) {  # 根据你的研究需要调整最小样本量
        # 口腔菌群 → 肠道菌群 → 代谢物
        tryCatch({
          med_model <- lm(gut ~ oral, data = med_data_complete)
          out_model <- lm(metabolite ~ oral + gut, data = med_data_complete)  # 移除交互项以简化模型
          
          # 进行中介分析 - 使用非交互式bootstrap
          med_result <- mediate(med_model, out_model, treat = "oral", mediator = "gut",
                                boot = TRUE, sims = 10)
          
          # 保存结果
          result_row <- data.frame(
            direction = "oral->gut->metabolite",
            oral_feature = oral_feature,
            gut_feature = gut_feature,
            metabolite = metabolite,
            ACME = med_result$d1,
            ACME_p = med_result$d1.p,
            prop_mediated = med_result$n1,
            stringsAsFactors = FALSE
          )
          
          bidirectional_mediation_results <- rbind(bidirectional_mediation_results, result_row)
        }, error = function(e) {
          # 记录错误但继续循环
          cat("Error in mediation analysis for:", oral_feature, gut_feature, metabolite, "\n")
          cat("Error message:", e$message, "\n")
        })
      }
    }
  }
}





######

for (i in 1:nrow(oral_gut_associations)) {
  oral_idx <- oral_gut_associations[i, "row"]
  gut_idx <- oral_gut_associations[i, "col"]
  
  oral_feature <- rownames(oral_gut_cor$cor)[oral_idx]
  gut_feature <- colnames(oral_gut_cor$cor)[gut_idx]
  
  # 检查该肠道菌群是否与代谢物相关联
  oral_metabolite_associations <- which(oral_metabolome_cor$p[oral_idx, ] < 0.05, arr.ind = TRUE)
  
  if (length(oral_metabolite_associations) > 0) {
    for (j in 1:length(oral_metabolite_associations)) {
      met_idx <- oral_metabolite_associations[j]
      metabolite <- colnames(oral_metabolome_cor$cor)[met_idx]
      
      # 构建数据框用于中介分析
      med_data <- data.frame(
        oral = oral_data_t[, oral_feature],
        gut = gut_data_t[, gut_feature],
        metabolite = metabolome_data_t[, metabolite]
      )
      
      # 在进行分析之前移除任何含有NA的行
      med_data_complete <- na.omit(med_data)
      
      # 确保有足够的观测进行分析
      if (nrow(med_data_complete) >= 10) {  # 根据你的研究需要调整最小样本量
        # 口腔菌群 → 肠道菌群 → 代谢物
        tryCatch({
          med_model <- lm(oral ~ gut, data = med_data_complete)
          out_model <- lm(metabolite ~ gut + oral , data = med_data_complete)
          
          # 进行中介分析 - 使用非交互式bootstrap或将interaction设为FALSE
          med_result <- mediate(med_model, out_model, treat = "gut", mediator = "oral",boot = TRUE, sims = 10)
          # 保存结果
          result_row <- data.frame(
            direction = "gut->oral->metabolite",
            oral_feature = oral_feature,
            gut_feature = gut_feature,
            metabolite = metabolite,
            ACME = med_result$d1,
            ACME_p = med_result$d1.p,
            prop_mediated = med_result$n1,
            stringsAsFactors = FALSE
          )
          
          bidirectional_mediation_results <- rbind(bidirectional_mediation_results, result_row)
        }, error = function(e) {
          # 记录错误但继续循环
          cat("Error in mediation analysis for:", oral_feature, gut_feature, metabolite, "\n")
          cat("Error message:", e$message, "\n")
        })
      }
    }
  }
}





# FDR校正
if (nrow(bidirectional_mediation_results) > 0) {
  bidirectional_mediation_results$ACME_fdr <- 
    p.adjust(bidirectional_mediation_results$ACME_p, method = "BH")
}

# 统计两个方向的存在显著交互效应的比例

bidirectional_mediation_results_sig<-subset(bidirectional_mediation_results,ACME_p<0.1)


oral_tax<-oral_temp_object@variable_info
gut_tax<-gut_temp_object@variable_info
# FDR校正
if (nrow(bidirectional_mediation_results) > 0) {
  bidirectional_mediation_results$ACME_fdr <- 
    p.adjust(bidirectional_mediation_results$ACME_p, method = "BH")
}





# 统计direction列的频数
direction_counts <- table(bidirectional_mediation_results_sig$direction)




ggplot(data=as.data.frame(direction_counts), aes(x=Var1, y=Freq,fill=Var1)) +
  geom_bar(stat="identity") +
  labs(x="Direction", y="Counts") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle=45, hjust=1,size=14),
        axis.text.y = element_text(size=12))+scale_fill_manual(values = c("#Edd064","#a1d5b9"))

bidirectional_mediation_results_sig<-subset(bidirectional_mediation_results,ACME_p<0.05)


bidirectional_mediation_results_sig<-merge(bidirectional_mediation_results_sig,oral_tax[,c("Genus","MGSS")],by.x="oral_feature",by.y="MGSS")

bidirectional_mediation_results_sig<-merge(bidirectional_mediation_results_sig,gut_tax[,c("Genus","MGSG")],by.x="gut_feature",by.y="MGSG")

bidirectional_mediation_results_sig<-subset(bidirectional_mediation_results_sig,!(Genus.x=="g__unknown"))

bidirectional_mediation_results_sig<-subset(bidirectional_mediation_results_sig,!(Genus.y%in%c("g__unknown","g__Firmicutes_unclassified")))


sankey_data<-bidirectional_mediation_results_sig[,c(8,9,4)]

colnames(sankey_data)<-c("Oral","Gut","Metabolite")
sankey_data$Oral<-gsub("g__","Oral_",sankey_data$Oral)
sankey_data$Gut<-gsub("g__","Gut_",sankey_data$Gut)
# 安装并加载必要的包

library(networkD3)
library(dplyr)

# 假设你的数据已经存在名为sankey_data的数据框中
# sankey_data包含三列: "Oral", "Gut", "Metabolite"

# 第一步：准备节点数据
# 收集所有唯一的节点名称
oral_nodes <- unique(sankey_data$Oral)
gut_nodes <- unique(sankey_data$Gut)
metabolite_nodes <- unique(sankey_data$Metabolite)

# 创建节点数据框
nodes_df <- data.frame(
  name = c(oral_nodes, gut_nodes, metabolite_nodes),
  group = c(rep("Oral", length(oral_nodes)), 
            rep("Gut", length(gut_nodes)), 
            rep("Metabolite", length(metabolite_nodes))),
  stringsAsFactors = FALSE
)

# 为每个节点分配一个唯一ID
nodes_df$ID <- 0:(nrow(nodes_df) - 1)

# 第二步：准备连接(links)数据
# 创建从Oral到Gut的连接
links_oral_gut <- sankey_data %>%
  dplyr::select(Oral, Gut) %>%
  dplyr::group_by(Oral, Gut) %>%
  dplyr::summarise(value = n(), .groups = 'drop') %>%
  dplyr::rename("source" = "Oral", "target" = "Gut")

# 创建从Gut到Metabolite的连接
links_gut_metabolite <- sankey_data %>%
  dplyr::select(Gut, Metabolite) %>%
  dplyr::group_by(Gut, Metabolite) %>%
  dplyr::summarise(value = n(), .groups = 'drop') %>%
  dplyr::rename("source" = "Gut", "target" = "Metabolite")

# 合并所有连接
# 为links添加组标识
links_oral_gut$group <- "Oral"
links_gut_metabolite$group <- "Gut"
links_df <- rbind(links_oral_gut, links_gut_metabolite)

# 将节点名称转换为节点ID
links_df$source <- match(links_df$source, nodes_df$name) - 1
links_df$target <- match(links_df$target, nodes_df$name) - 1

# 定义颜色函数，为每个组分配固定颜色
colourScale <- JS(paste0('d3.scaleOrdinal()
  .domain(["Oral", "Gut", "Metabolite"])
  .range(["#a1d5b9", "#Edd064", "#B6C7EA"])'))

# 第三步：创建桑基图
sankey_plot<-sankeyNetwork(
  Links = links_df, 
  Nodes = nodes_df,
  Source = "source", 
  Target = "target",
  Value = "value", 
  NodeID = "name",
  NodeGroup = "group",  # 根据组分配颜色
  LinkGroup = "group",  # 根据源节点的组分配连接颜色
  colourScale = colourScale,  # 使用自定义颜色比例
  fontSize = 12,
  nodeWidth = 30,
  nodePadding = 10,
  margin = list(left = 50, right = 50),
  sinksRight = TRUE,
  units = "个数"
)

final_sankey <- htmlwidgets::onRender(sankey_plot, nodePositionJS)

# 保存为HTML文件
html_file <- "sankey_diagram.html"
htmlwidgets::saveWidget(final_sankey, html_file)

# 将HTML文件转换为PDF
pdf_file <- "sankey_diagram.pdf"
webshot2::webshot(
  url = html_file, 
  file = pdf_file,
  delay = 1,        # 等待渲染完成的时间（秒）
  zoom = 2,         # 增加分辨率
  vwidth = 800,     # 视口宽度
  vheight = 600     # 视口高度
)

calculate_correlations <- function(gut_data, oral_data, metabolome_data, 
                                   mediation_results) {
  
  # 创建一个空的列表来存储结果
  correlation_results <- list()
  
  # 对mediation_results的每一行进行循环
  for (i in 1:nrow(mediation_results)) {
    # 获取当前行的特征名称
    gut_feature <- mediation_results$gut_feature[i]
    oral_feature <- mediation_results$oral_feature[i]
    metabolite_feature <- mediation_results$metabolite[i]
    
    # 提取相应的数据
    gut_values <- gut_data[gut_feature, ]
    oral_values <- oral_data[oral_feature, ]
    metabolite_values <- metabolome_data[metabolite_feature, ]
    
    # 确保所有数据都是数值型
    gut_values <- as.numeric(gut_values)
    oral_values <- as.numeric(oral_values)
    metabolite_values <- as.numeric(metabolite_values)
    
    # 创建一个组合数据框，只包含有完整观测的样本
    combined_data <- data.frame(
      gut = gut_values,
      oral = oral_values,
      metabolite = metabolite_values
    )
    
    # 移除含有NA的行
    combined_data <- na.omit(combined_data)
    
    # 如果有足够的数据点来计算相关性
    if (nrow(combined_data) >= 3) {
      # 计算Pearson相关系数
      cor_gut_oral <- cor.test(combined_data$gut, combined_data$oral, method = "spearman")
      cor_gut_metabolite <- cor.test(combined_data$gut, combined_data$metabolite, method = "spearman")
      cor_oral_metabolite <- cor.test(combined_data$oral, combined_data$metabolite, method = "spearman")
      
      # 存储结果
      result <- data.frame(
        Mediation_Row = i,
        Gut_Feature = gut_feature,
        Oral_Feature = oral_feature,
        Metabolite_Feature = metabolite_feature,
        Gut_Oral_Cor = cor_gut_oral$estimate,
        Gut_Oral_Pvalue = cor_gut_oral$p.value,
        Gut_Metabolite_Cor = cor_gut_metabolite$estimate,
        Gut_Metabolite_Pvalue = cor_gut_metabolite$p.value,
        Oral_Metabolite_Cor = cor_oral_metabolite$estimate,
        Oral_Metabolite_Pvalue = cor_oral_metabolite$p.value,
        Sample_Size = nrow(combined_data)
      )
      
      # 添加到结果列表
      correlation_results[[i]] <- result
    } else {
      # 如果数据点不足，添加一个包含NA的行
      result <- data.frame(
        Mediation_Row = i,
        Gut_Feature = gut_feature,
        Oral_Feature = oral_feature,
        Metabolite_Feature = metabolite_feature,
        Gut_Oral_Cor = NA,
        Gut_Oral_Pvalue = NA,
        Gut_Metabolite_Cor = NA,
        Gut_Metabolite_Pvalue = NA,
        Oral_Metabolite_Cor = NA,
        Oral_Metabolite_Pvalue = NA,
        Sample_Size = nrow(combined_data)
      )
      
      # 添加到结果列表
      correlation_results[[i]] <- result
    }
  }
  
  # 将所有结果合并成一个数据框
  final_results <- do.call(rbind, correlation_results)
  
  return(final_results)
}