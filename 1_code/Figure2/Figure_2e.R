setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")
library(ggtree)
library(tidyverse)
library(tidymass)
library(ade4)
library(ggnewscale)
library(ggtreeExtra)
all_cors<-readRDS("3_data_analysis/4_site_merge/all_cors")

all_genus_pc1_cor<-readRDS("3_data_analysis/4_site_merge/all_genus_pc1_cor")


all_genus_pc1_cor<-subset(all_genus_pc1_cor,Kingdom=="Bacteria")

taxonomy_data<-all_genus_pc1_cor[,c("Kingdom","Phylum","Class","Order","Family","Genus")]

taxonomy_data<-subset(taxonomy_data,Phylum%in%c("Actinobacteria","Bacteroidetes","Firmicutes","Proteobacteria"))

taxonomy_data_test<-taxonomy_data[,c("Genus","Family","Order","Class","Phylum")]

taxonomy_data_test <- data.frame(lapply(taxonomy_data_test, as.factor))
rownames(taxonomy_data_test)<-taxonomy_data_test$Genus

taxonomy_data_test<-taxonomy_data_test[,-1]
tax.phy<-taxo2phylog(as.taxo(taxonomy_data_test),abbrev = FALSE)
p<-ggtree(tax.phy, layout="fan", size=0.15, open.angle=90)

tips_data <- data.frame(
  Genus = names(tax.phy$leaves),
  Phylum = taxonomy_data_test$Phylum[match(names(tax.phy$leaves), rownames(taxonomy_data_test))]
)

tips_data$Phylum[135]<-"Proteobacteria"

library(reshape2)

p1<-p %<+% tips_data + 
  geom_tippoint(aes(color = Phylum), size = 2) +
  scale_color_discrete(name = "Phylum") +
  theme(legend.position = "right")+scale_color_manual(values = phylum_color)

abundance_matrix<-subset(all_genus_pc1_cor,Genus%in%tips_data$Genus)
rownames(abundance_matrix)<-abundance_matrix$Genus

abundance_matrix<-abundance_matrix[tips_data$Genus,]
abundance_matrix<-abundance_matrix[,1:5]

abundance_matrix<-melt(abundance_matrix,id.vars = "Genus")

p2<-p1+new_scale_fill() +
  geom_fruit(data=abundance_matrix, geom=geom_tile,
             mapping=aes(y=Genus, x=variable, alpha=value, fill=variable),
             color = "grey70", offset = 0.04,size = 0.02)+
  scale_alpha_continuous(range=c(0, 1.5),
                         guide=guide_legend(keywidth = 0.3, 
                                            keyheight = 0.3, order=5))+
  scale_fill_manual(values = body_site_color)




all_cors$Microbiome<-gsub("Gut_","",all_cors$Microbiome)
all_cors$Microbiome<-gsub("Oral_","",all_cors$Microbiome)
all_cors$Microbiome<-gsub("Skin_","",all_cors$Microbiome)
all_cors$Microbiome<-gsub("Nasal_","",all_cors$Microbiome)




data3<-data.frame(table(all_cors$Microbiome,all_cors$Site))

colnames(data3)<-c("Genus","site","Freq")

data3$site<-tolower(data3$site)

data3<-subset(data3,Genus%in%tips_data$Genus)

plot <-
p2+new_scale_fill()+geom_fruit(data=data3, geom=geom_bar,
                               mapping=aes(y=Genus, x=Freq, fill=site),
                               pwidth=0.38, 
                               orientation="y", 
                               stat="identity",offset = 0.04
)+scale_fill_manual(values = body_site_color)

ggsave(plot, 
       filename = "4_manuscript/Figures/Figure_2/figure_2e.pdf", width = 10, height = 8)
