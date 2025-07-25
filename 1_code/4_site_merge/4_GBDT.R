rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")
library(ggbeeswarm)
setwd("1_code/4_site_merge/")



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
g1 <- ggplot(data_long, aes(x = factor(metabolite, levels = data$metabolite), y = R_squared, fill = Factor)) +
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

# 使用函数
 result <- metabolite_analysis(four_site_GBDT_R2)
 print(result)




  
  
  
  
  
### UpSet
  
  library(UpSetR)
  four_site_GBDT_R2

  metabolite_annotation_mainclass<-subset(metabolite_annotation,HMDB.Class%in%c("Benzene and substituted derivatives","Carboxylic acids and derivatives","Fatty Acyls","Glycerophospholipids","Indoles and derivatives","Organooxygen compounds","Steroids and steroid derivatives"))
  
  
  df<-four_site_GBDT_R2
  
  df <- df %>% mutate(across(everything(), ~ifelse(. > 0.05, 1, 0)))
  
  df$metabolite<-four_site_GBDT_R2$metabolite
  
  
  upset(df, sets = c("gut", "oral", "skin", "nasal"), keep.order = TRUE,sets.bar.color =c("#edd064" , "#a1d5b9" ,"#f2ccac","#a17db4"))
  
  
  
 
  
  
  
  
  
  
  
  
