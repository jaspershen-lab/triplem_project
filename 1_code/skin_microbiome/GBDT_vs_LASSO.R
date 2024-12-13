# 对比GBDT和LASSO的效果
rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")

library(tidyverse)
library(tidymass)
library(plyr)
library(microbiomedataset)


# 导入数据

GBDT_results<-readRDS("3_data_analysis/skin_microbiome/GBDT/cross_section/skin_GBDT_results")

LASSO_results<-readRDS("3_data_analysis/skin_microbiome/Lasso/cross_section/skin_lasso_results")


GBDT_results<-GBDT_results$summary
LASSO_results<-LASSO_results$results_df



GBDT_results<-subset(GBDT_results,r2_mean>0.01)
LASSO_results<-subset(LASSO_results,test_r2>0.01)


rownames(GBDT_results)<-GBDT_results$metabolite
rownames(LASSO_results)<-LASSO_results$metabolite

common_metabolite<-intersect(GBDT_results$metabolite,LASSO_results$metabolite)


GBDT_results<-GBDT_results[common_metabolite,]
LASSO_results<-LASSO_results[common_metabolite,]

two_model<-cbind(GBDT_results[,1:2],LASSO_results[,2])

colnames(two_model)<-c("metabolite","GBDT_R2","lasso_R2")


two_model$delt_R2<-two_model$GBDT_R2-two_model$lasso_R2

# 绘图


p1<-ggplot(two_model, aes(x=lasso_R2, y=GBDT_R2)) +
  geom_point(shape=21, size=4, fill="#f2ccac", color="white") +
  geom_abline(intercept=0, slope=1, color="grey50", linetype="dashed", size=1) +  # 添加对角线
  scale_x_continuous(limits=c(0, 0.2)) +
  scale_y_continuous(limits=c(0, 0.2)) +
  theme_light()+
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14,family = "Helvetica"),
    axis.title = element_text(size = 14,family = "Helvetica"),
    axis.text.x = element_text(family = "Helvetica") , # 如果组名较长，可以倾斜x轴标签
    axis.ticks.length = unit(0.25, "cm"),  # 增加刻度线长度
    axis.ticks = element_line(linewidth = 0.8)  # 增加刻度线粗细
  ) +xlab("lasso R²")+ylab("GBDT R²")




p2<-ggplot(two_model, aes(x=delt_R2)) +
  geom_histogram(binwidth=0.01, fill="#f2ccac", color="white") +  # 绘制直方图
  geom_vline(aes(xintercept=median(two_model$delt_R2, na.rm=TRUE)),  # 添加中位数线
             color="grey50", linetype="dashed", size=1) +
  labs(title="Frequency Distribution", x="Variable", y="Frequency") +
  theme_light()+
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14,family = "Helvetica"),
    axis.title = element_text(size = 14,family = "Helvetica"),
    axis.text.x = element_text(family = "Helvetica") , # 如果组名较长，可以倾斜x轴标签
    axis.ticks.length = unit(0.25, "cm"),  # 增加刻度线长度
    axis.ticks = element_line(linewidth = 0.8)  # 增加刻度线粗细
  ) +xlab("GBDT R²-lasso R²")