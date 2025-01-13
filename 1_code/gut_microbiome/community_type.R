# 计算肠型

library(microbiome)
library(DirichletMultinomial)
library(reshape2)
library(magrittr)
library(dplyr)


library(microbiomedataset)
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

gut_object <-
  gut_object %>%
  summarize_variables(what = "sum_intensity", group_by = "Genus")







## 导入成phyloseq对象



OTU = otu_table(gut_object@expression_data, taxa_are_rows = TRUE)
TAX<-as.matrix(data.frame(gut_object@variable_info))
TAX = tax_table(TAX)
samples<-data.frame(gut_object@sample_info)
rownames(samples)<-samples[,1]
samples = sample_data(samples)

gut_phyloseq<- phyloseq(OTU, TAX, samples)


pseq <- gut_phyloseq

pseq.comp <- microbiome::transform(pseq, "compositional")


taxa <- core_members(pseq.comp, prevalence = 5/100,detection = 1/100)
pseq <- prune_taxa(taxa, pseq)

# Pick the OTU count matrix
# and convert it into samples x taxa format
dat <- abundances(pseq)
count <- as.matrix(t(dat))


fit <- lapply(1:6, dmn, count = count, verbose=TRUE)

lplc <- base::sapply(fit, DirichletMultinomial::laplace) # AIC / BIC / Laplace
aic  <- base::sapply(fit, DirichletMultinomial::AIC) # AIC / BIC / Laplace
bic  <- base::sapply(fit, DirichletMultinomial::BIC) # AIC / BIC / Laplace

plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
#lines(aic, type="b", lty = 2)
lines(bic, type="b", lty = 3)

best <- fit[[which.min(unlist(lplc))]]


ass <- apply(mixture(best), 1, which.max)


for (k in seq(ncol(fitted(best)))) {
  d <- melt(fitted(best))
  colnames(d) <- c("OTU", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    # Arrange OTUs by assignment strength
    arrange(value) %>%
    mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
    # Only show the most important drivers
    filter(abs(value) > quantile(abs(value), 0.8))     
  d<-merge(d,data.frame(gut_object@variable_info),by.x="OTU",by.y="variable_id")
  p <- ggplot(d, aes(x = Genus, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top drivers: community type", k))
  print(p)
}



# 绘制PCoA

gut_phyloseq <- microbiome::transform(gut_phyloseq, "compositional")
taxa <- core_members(gut_phyloseq, prevalence = 5/100,detection = 1/100)
gut_phyloseq <- prune_taxa(taxa, gut_phyloseq)
gut_phyloseq_filt<-gut_phyloseq

sample_data(gut_phyloseq_filt)$CT <- as.character(ass)


dis_bray<- phyloseq::distance(gut_phyloseq_filt, "bray")

## 采用PCoA的方法对距离矩阵进行降维
dis_bray.pcoa = ordinate(gut_phyloseq_filt, method="NMDS", distance=dis_bray)

## 绘制初始图形
bray.pcoa <- plot_ordination(gut_phyloseq_filt, dis_bray.pcoa, color="IRIS" ) + geom_point(size=3)

## 提取图形数据
data<-bray.pcoa$data

## 修改列名
colnames(data)[1:2]<-c("NMDS1","NMDS2")

## 获取主坐标轴1,2的解释度
pc1<-""
pc2<-""



ggplot(data, aes(NMDS1, NMDS2)) +
  #绘制样本点，根据分组匹配颜色和形状，size调整点的大小
  geom_point(aes(colour=CT,fill=CT),size=2.5)+
  #匹配形状、边框和填充的图例
  scale_color_manual(values=c ("#2442B2","#E12441"))+
  #设置标题和横纵坐标label文字
  labs(title="NMDS - The composition of gut microbiome") +
  theme(text=element_text(size=30))+
  #添加横纵两条虚线
  geom_vline(aes(xintercept = 0),linetype="dotted")+
  geom_hline(aes(yintercept = 0),linetype="dotted")+
  #调整背景、坐标轴、图例的格式
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid=element_blank(),
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour 
                                 = "black"),
        axis.title.x=element_text(colour='black', size=14),
        axis.title.y=element_text(colour='black', size=14),
        axis.text=element_text(colour='black',size=14),
        legend.title=element_blank(),
        legend.text=element_text(size=14),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "black"),
        legend.key.height=unit(1.6,"cm"))+
  #设置标题的格式
  theme(plot.title = element_text(size=14,colour = "black",hjust = 0.5,face = "bold"))+stat_ellipse(aes(color = CT),geom = "polygon",level = 0.8,alpha = 0,size=2)

