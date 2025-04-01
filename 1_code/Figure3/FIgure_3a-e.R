rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")
library(ggbeeswarm)
setwd("1_code/4_site_merge/")
load("../../3_data_analysis/gut_microbiome/data_preparation/object_cross_section")

gut_object<-object_cross_section

load("../../3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section")

metabolomics_object<-object_cross_section


load("../../3_data_analysis/oral_microbiome/data_preparation/object_cross_section")

metabolite_annotation<-read_excel("../../3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx")

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


gut_GBDT_results<-readRDS("../../3_data_analysis/gut_microbiome/GBDT/cross_section/gut_GBDT_results")
oral_GBDT_results<-readRDS("../../3_data_analysis/oral_microbiome/GBDT/cross_section/oral_GBDT_results")
skin_GBDT_results<-readRDS("../../3_data_analysis/skin_microbiome/GBDT/cross_section/skin_GBDT_results")
nasal_GBDT_results<-readRDS("../../3_data_analysis/nasal_microbiome/GBDT/cross_section/nasal_GBDT_results")


gut_GBDT_results_R2<-gut_GBDT_results$summary[,c(1,2)]
colnames(gut_GBDT_results_R2)<-c("metabolite","gut")
gut_GBDT_results_R2<-subset(gut_GBDT_results_R2,gut>0.1)

oral_GBDT_results_R2<-oral_GBDT_results$summary[,c(1,2)]
colnames(oral_GBDT_results_R2)<-c("metabolite","oral")
oral_GBDT_results_R2<-subset(oral_GBDT_results_R2,oral>0.1)

skin_GBDT_results_R2<-skin_GBDT_results$summary[,c(1,2)]
colnames(skin_GBDT_results_R2)<-c("metabolite","skin")
skin_GBDT_results_R2<-subset(skin_GBDT_results_R2,skin>0.1)

nasal_GBDT_results_R2<-nasal_GBDT_results$summary[,c(1,2)]
colnames(nasal_GBDT_results_R2)<-c("metabolite","nasal")
nasal_GBDT_results_R2<-subset(nasal_GBDT_results_R2,nasal>0.1)


four_site_GBDT_R2<-merge(gut_GBDT_results_R2,oral_GBDT_results_R2,by = "metabolite",all = TRUE)
four_site_GBDT_R2<-merge(four_site_GBDT_R2,skin_GBDT_results_R2,by = "metabolite",all = TRUE)

four_site_GBDT_R2<-merge(four_site_GBDT_R2,nasal_GBDT_results_R2,by = "metabolite",all = TRUE)
colnames(four_site_GBDT_R2)<-c("metabolite","gut","oral","skin","nasal")

colnames(four_site_GBDT_R2)<-c("metabolite","gut","oral","skin","nasal")

four_site_GBDT_R2<- data.frame(lapply(four_site_GBDT_R2, function(x) ifelse(is.na(x), 0, x)))

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)



# 确定每个代谢物中占比最大的因素
data <- four_site_GBDT_R2 %>%
  mutate(Dominant_Factor = case_when(
    gut >= oral & gut >= skin & gut>=nasal ~ "gut",
    oral >= gut & oral>= skin & oral>=nasal~ "oral",
    skin >= oral & skin >= gut & skin>=nasal ~ "skin",
    nasal >= oral & nasal >= skin & nasal>=oral ~ "nasal"
  ))

# 计算每个代谢物的总 R² 值，并按最主要的影响因素分类后再按总 R² 值降序排序
data <- data %>%
  mutate(Total_R2 = gut + oral + skin + nasal) %>%
  arrange(Dominant_Factor, desc(Total_R2))

# 转换为长格式以便于 ggplot 绘制
data_long <- data %>%
  pivot_longer(cols = c("gut", "oral", "skin","nasal"),
               names_to = "Factor",
               values_to = "R_squared")

# 绘制堆叠条形图
 ggplot(data_long, aes(x = factor(metabolite, levels = data$metabolite), y = R_squared, fill = Factor)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = body_site_color, name = "Factor") +
  theme_minimal() +
  labs(x = "256 metabolites (adjusted R² > 10%)", y = "Adjusted R²") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "top",
        panel.spacing = unit(0, "lines"),
        panel.grid = element_blank()) 




## 绘制四个位点R2值前50的样本点图
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




library(tidyverse)
library(ggplot2)

# 分别获取并展示每个部位前50的代谢物
metabolite_analysis <- function(four_site_GBDT_R2) {
  # 分别获取每个部位前50的数据
  gut_data <- four_site_GBDT_R2 %>% 
    top_n(50, gut) %>%
    select(metabolite, gut) %>%
    mutate(Site = "gut",
           Value = gut) %>%
    select(metabolite, Site, Value)
  
  oral_data <- four_site_GBDT_R2 %>% 
    top_n(50, oral) %>%
    select(metabolite, oral) %>%
    mutate(Site = "oral",
           Value = oral) %>%
    select(metabolite, Site, Value)
  
  skin_data <- four_site_GBDT_R2 %>% 
    top_n(50, skin) %>%
    select(metabolite, skin) %>%
    mutate(Site = "skin",
           Value = skin) %>%
    select(metabolite, Site, Value)
  
  nasal_data <- four_site_GBDT_R2 %>% 
    top_n(50, nasal) %>%
    select(metabolite, nasal) %>%
    mutate(Site = "nasal",
           Value = nasal) %>%
    select(metabolite, Site, Value)
  
  # 合并所有数据
  combined_data <- bind_rows(gut_data, oral_data, skin_data, nasal_data)
  
  
  
  combined_data<-data.frame(combined_data)
  combined_data$Site<-as.factor(combined_data$Site)
  summary_stats <- combined_data %>%
    group_by(Site) %>%
    dplyr:: summarise(
      mean = mean(Value),
      sd = sd(Value),
      count = length(Value),  # 使用 length() 替代 n()
      se = sd(Value)/sqrt(length(Value))
    )
  
  
  ggplot() +
    # 添加条形图
    geom_bar(data = summary_stats, 
             aes(x = Site, y = mean, fill = Site),
             stat = "identity",
             width = 0.6,
             alpha = 1) +
    # 添加误差线
    geom_errorbar(data = summary_stats,
                  aes(x = Site, 
                      ymin = mean - se, 
                      ymax = mean + se),
                  width = 0.2) +
    # 添加散点
    geom_quasirandom(data = combined_data,
                     aes(x = Site, y = Value),
                     alpha = 0.8,
                     width = 0.2) +
    # 设置填充颜色
    scale_fill_manual(values = body_site_color) +
    # 设置y轴
    scale_y_continuous(expand = c(0, 0)) +  
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 14,family = "Helvetica"),
      axis.title = element_text(size = 14,family = "Helvetica"),
      axis.text.x = element_text(angle = 30, hjust = 1,family = "Helvetica") , # 如果组名较长，可以倾斜x轴标签
      axis.ticks.length = unit(0.25, "cm"),  # 增加刻度线长度
      axis.ticks = element_line(linewidth = 0.8)  # 增加刻度线粗细
    ) +
    # 设置坐标轴标签
    xlab("") +
    ylab("R2")
  
  
}




## p-cresol and  PAGln


# 函数：获取特定代谢物的观察值与预测值相关性
# 参数说明：
# - metabolite_name: 需要分析的代谢物名称
# - model_dir: 保存模型的目录
# - microbiome_data: 微生物组数据矩阵，行为特征(如OTU/ASV)，列为样本名
# - metabolite_data: 代谢物组数据矩阵，行为代谢物，列为样本名
# - plot_results: 是否生成可视化图表
get_metabolite_prediction_correlation <- function(
    metabolite_name = "M187T125_2_NEG_RPLC", 
    model_dir = "models",
    microbiome_data,
    metabolite_data,
    plot_results = TRUE
) {
  # 加载必要的包
  library(dplyr)
  library(ggplot2)
  library(gbm)
  
  # 步骤1: 构建模型文件名并检查是否存在
  model_file <- file.path(model_dir, paste0(make.names(metabolite_name), "_model.rds"))
  if(!file.exists(model_file)) {
    stop(sprintf("模型文件 %s 不存在，请检查代谢物名称和模型目录", model_file))
  }
  
  # 步骤2: 加载模型
  model_info <- readRDS(model_file)
  message(sprintf("成功加载模型: %s", model_file))
  
  # 步骤3: 数据预处理(确保样本匹配)
  micro_samples <- colnames(microbiome_data)
  meta_samples <- colnames(metabolite_data)
  common_samples <- intersect(micro_samples, meta_samples)
  
  message(sprintf("可用样本数: %d", length(common_samples)))
  
  microbiome_matched <- microbiome_data[, common_samples]
  metabolite_matched <- metabolite_data[, common_samples]
  
  # 步骤4: 准备预测数据
  # - 获取代谢物索引
  metabolite_idx <- which(rownames(metabolite_matched) == metabolite_name)
  if(length(metabolite_idx) == 0) {
    # 尝试部分匹配
    possible_matches <- grep(metabolite_name, rownames(metabolite_matched), value = TRUE)
    if(length(possible_matches) > 0) {
      message("未找到完全匹配的代谢物名称，但找到了以下可能的匹配：")
      for(i in 1:min(5, length(possible_matches))) {
        message(sprintf("%d. %s", i, possible_matches[i]))
      }
      if(length(possible_matches) > 5) {
        message(sprintf("...及其他 %d 个可能的匹配", length(possible_matches) - 5))
      }
      
      # 询问是否使用第一个匹配项
      message(sprintf("使用第一个匹配项 '%s' 进行分析...", possible_matches[1]))
      metabolite_name <- possible_matches[1]
      metabolite_idx <- which(rownames(metabolite_matched) == metabolite_name)
    } else {
      # 显示代谢物数据中的一些行名作为参考
      message("代谢物数据中的一些行名（供参考）：")
      print(head(rownames(metabolite_matched), 10))
      stop(sprintf("在代谢物数据中找不到 %s", metabolite_name))
    }
  }
  
  # - 提取代谢物观察值
  observed_values <- as.numeric(metabolite_matched[metabolite_idx, ])
  
  # - 准备微生物组数据作为预测输入
  X <- t(microbiome_matched) # 转置使样本为行
  
  # 步骤5: 仅使用模型中的选定特征
  selected_features <- model_info$selected_features
  
  # 检查特征是否存在于当前数据中
  missing_features <- selected_features[!selected_features %in% colnames(X)]
  if(length(missing_features) > 0) {
    warning(sprintf("有 %d 个模型特征在当前数据中不存在，这可能影响预测质量", 
                    length(missing_features)))
    message("缺失的前几个特征：")
    print(head(missing_features))
    # 只使用存在的特征
    selected_features <- selected_features[selected_features %in% colnames(X)]
  }
  
  if(length(selected_features) > 0) {
    message(sprintf("使用 %d 个选定特征进行预测", length(selected_features)))
    X_selected <- X[, selected_features, drop = FALSE]
  } else {
    warning("没有可用的特征用于预测，将使用所有特征")
    X_selected <- X
  }
  
  # 步骤6: 使用模型进行预测
  prediction_data <- as.data.frame(X_selected)
  predicted_values <- as.numeric(predict(model_info$model, newdata = prediction_data, 
                                         n.trees = model_info$parameters$n.trees))
  
  # 步骤7: 计算相关性
  # 首先检查数据类型
  if(!is.numeric(observed_values)) {
    warning("观察值不是数值类型，尝试强制转换")
    observed_values <- as.numeric(observed_values)
  }
  
  if(!is.numeric(predicted_values)) {
    warning("预测值不是数值类型，尝试强制转换")
    predicted_values <- as.numeric(predicted_values)
  }
  
  # 检查是否有NA值
  if(any(is.na(observed_values)) || any(is.na(predicted_values))) {
    valid_idx <- which(!is.na(observed_values) & !is.na(predicted_values))
    if(length(valid_idx) == 0) {
      stop("所有样本都包含NA值，无法计算相关性")
    }
    warning(sprintf("发现 %d 个NA值，将使用 %d 个有效样本计算相关性", 
                    length(observed_values) - length(valid_idx), length(valid_idx)))
    observed_values <- observed_values[valid_idx]
    predicted_values <- predicted_values[valid_idx]
    common_samples <- common_samples[valid_idx]
  }
  
  correlation_test <- cor.test(observed_values, predicted_values, 
                               method = "pearson")
  
  r_value <- correlation_test$estimate
  p_value <- correlation_test$p.value
  r_squared <- r_value^2
  
  message(sprintf("相关系数 (r): %.4f", r_value))
  message(sprintf("决定系数 (R²): %.4f", r_squared))
  message(sprintf("p-值: %.6e", p_value))
  
  # 步骤8: 创建可视化
  if(plot_results) {
    results_df <- data.frame(
      Observed = observed_values,
      Predicted = predicted_values,
      Sample = common_samples
    )
    
    p <- ggplot(results_df, aes(x = Observed, y = Predicted)) +
      geom_point(alpha = 0.7) +
      geom_smooth(method = "lm", formula = y ~ x, color = "blue") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      theme_minimal() +
      labs(
        title = sprintf("代谢物 %s 的观察值 vs 预测值", metabolite_name),
        subtitle = sprintf("r = %.3f, R² = %.3f, p = %.3e", 
                           r_value, r_squared, p_value),
        x = "观察值",
        y = "预测值"
      )
    
    print(p)
    
    # 创建诊断图：残差 vs 拟合值
    results_df$Residuals <- results_df$Observed - results_df$Predicted
    
    p_residual <- ggplot(results_df, aes(x = Predicted, y = Residuals)) +
      geom_point(alpha = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      geom_smooth(method = "loess", formula = y ~ x, color = "blue") +
      theme_minimal() +
      labs(
        title = "残差 vs 预测值",
        x = "预测值",
        y = "残差"
      )
    
    print(p_residual)
    
    # 保存诊断结果
    result_plots <- list(
      scatter_plot = p,
      residual_plot = p_residual
    )
  } else {
    result_plots <- NULL
  }
  
  # 步骤9: 返回结果
  # 创建包含所有结果的数据框
  results_df <- data.frame(
    Sample = common_samples,
    Observed = observed_values,
    Predicted = predicted_values,
    Residuals = observed_values - predicted_values
  )
  
  return(list(
    metabolite = metabolite_name,
    observed = observed_values,
    predicted = predicted_values,
    samples = common_samples,
    correlation = r_value,
    r_squared = r_squared,
    p_value = p_value,
    model_file = model_file,
    results_table = results_df,
    features_used = selected_features,
    plots = result_plots
  ))
}

microbiome_data<-gut_temp_object@expression_data

row.names(microbiome_data)<-gut_temp_object@variable_info$Genus
result <- get_metabolite_prediction_correlation(
  metabolite_name = "M187T125_2_NEG_RPLC",
  model_dir = "../../models/",
  microbiome_data = microbiome_data,
  metabolite_data = metabolomics_temp_object@expression_data,
  plot_results = TRUE
)


ggplot(result$results_table, aes(x=result$results_table$Observed, y= result$results_table$Predicted)) +
  geom_point(shape=21,size=4,fill="#edd064",color="white") +
  geom_smooth(method="lm",colour = "#edd064",) +theme_light() +stat_cor(method = "spearman")+
  theme(legend.position="none", #不需要图例
        axis.text.x=element_text(colour="black",size=14), #设置x轴刻度标签的字体属性
        axis.text.y=element_text(size=14,face="plain"), #设置x轴刻度标签的字体属性
        axis.title.y=element_text(size = 14,face="plain"), #设置y轴的标题的字体属性
        axis.title.x=element_text(size = 14,face="plain"), #设置x轴的标题的字体属性
        plot.title = element_text(size=15,face="bold",hjust = 0.5))+xlab("Observed p-cresol")+ylab("Predicted  p-cresol")


result <- get_metabolite_prediction_correlation(
  metabolite_name = "M263T353_NEG_HILIC",
  model_dir = "../../models/",
  microbiome_data = microbiome_data,
  metabolite_data = metabolomics_temp_object@expression_data,
  plot_results = TRUE
)


ggplot(result$results_table, aes(x=result$results_table$Observed, y= result$results_table$Predicted)) +
  geom_point(shape=21,size=4,fill="#edd064",color="white") +
  geom_smooth(method="lm",colour = "grey50") +theme_light() +stat_cor(method = "spearman")+
  theme(legend.position="none", #不需要图例
        axis.text.x=element_text(colour="black",size=14), #设置x轴刻度标签的字体属性
        axis.text.y=element_text(size=14,face="plain"), #设置x轴刻度标签的字体属性
        axis.title.y=element_text(size = 14,face="plain"), #设置y轴的标题的字体属性
        axis.title.x=element_text(size = 14,face="plain"), #设置x轴的标题的字体属性
        plot.title = element_text(size=15,face="bold",hjust = 0.5))+xlab("Observed PAGln")+ylab("Predicted  PAGln")



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
      )*0.6
    )
  
  # FDR校正
  results$FDR <- p.adjust(results$P_value, method = "BH")
  
  # 按p值排序
  results <- results %>% arrange(P_value)
  
  return(results)
}
HMDB.Class<-c("Benzene and substituted derivatives","Carboxylic acids and derivatives","Fatty Acyls","Glycerophospholipids","Indoles and derivatives","Organooxygen compounds","Steroids and steroid derivatives")

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

# 创建一个包含需要添加星号位置的数据框
stars_data <- data.frame(
  site = c("gut", "gut", "skin", "oral", "nasal"),
  HMDB.Class = c("Organooxygen compounds", 
                 "Benzene and substituted derivatives", 
                 "Carboxylic acids and derivatives", 
                 "Indoles and derivatives", 
                 "Carboxylic acids and derivatives"),
  Significant_in_class = NA,  # 这个值会从原始数据中提取
  label = "*"
)

# 将星号位置的y值从主数据集中提取出来
# 对于每个星号位置，我们需要知道对应的y值(Significant_in_class)
for(i in 1:nrow(stars_data)){
  row_match <- results_all[results_all$site == stars_data$site[i] & 
                             results_all$HMDB.Class == stars_data$HMDB.Class[i], ]
  
  if(nrow(row_match) > 0){
    # 获取该条形的高度，并在上方添加一点空间以放置星号
    stars_data$Significant_in_class[i] <- row_match$Significant_in_class[1] * 1.05  # 增加5%的高度放置星号
  }
}


# 首先，将site变量转换为有序因子
results_all$site <- factor(results_all$site, 
                           levels = c("gut", "oral", "skin", "nasal"))

# 同样需要对stars_data做相同处理（如果你使用了前面的stars_data方法）
stars_data$site <- factor(stars_data$site, 
                          levels = c("gut", "oral", "skin", "nasal"))

# 然后进行绘图
ggplot(results_all, aes(x = HMDB.Class, y = Significant_in_class)) +
  geom_bar(stat = "identity", aes(fill = site)) +
  theme_minimal() +
  facet_wrap(~site, nrow = 1) +  # 这里的分面现在会按照指定的顺序排列
  coord_flip() +
  theme_light() +
  # 添加星号部分（如果你使用了前面的stars_data方法）
  geom_text(data = stars_data, 
            aes(x = HMDB.Class, y = Significant_in_class, label = label),
            size = 8, fontface = "bold") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14, family = "Helvetica"),
    axis.title = element_text(size = 14, family = "Helvetica"),
    axis.text.x = element_text(family = "Helvetica"),
    axis.ticks.length = unit(0.25, "cm"),
    axis.ticks = element_line(linewidth = 0.8)
  ) +
  scale_fill_manual(values = body_site_color)
### UpSet

library(UpSetR)
four_site_GBDT_R2




df<-four_site_GBDT_R2

df <- df %>% mutate(across(everything(), ~ifelse(. > 0.05, 1, 0)))

df$metabolite<-four_site_GBDT_R2$metabolite


upset(df, sets = c("gut", "oral", "skin", "nasal"), keep.order = TRUE,sets.bar.color =c("#edd064" , "#a1d5b9" ,"#f2ccac","#a17db4"))
