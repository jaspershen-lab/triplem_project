# For GBDTLASSO.
rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")

library(tidyverse)
library(tidymass)
library(plyr)
library(microbiomedataset)


# Import the data.

GBDT_results<-readRDS("3_data_analysis/gut_microbiome/GBDT/cross_section/gut_GBDT_results")

LASSO_results<-readRDS("3_data_analysis/gut_microbiome/Lasso/cross_section/gut_lasso_results")


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

# Generate the plot.


p1<-ggplot(two_model, aes(x=lasso_R2, y=GBDT_R2)) +
  geom_point(shape=21, size=4, fill="#edd064", color="white") +
  geom_abline(intercept=0, slope=1, color="grey50", linetype="dashed", size=1) +  # Add the diagonal reference line.
  scale_x_continuous(limits=c(0, 0.3)) +
  scale_y_continuous(limits=c(0, 0.3)) +
  theme_light()+
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14,family = "Helvetica"),
    axis.title = element_text(size = 14,family = "Helvetica"),
    axis.text.x = element_text(family = "Helvetica") , # Rotate x-axis labels if group names are long.
    axis.ticks.length = unit(0.25, "cm"),  # Increase tick length.
    axis.ticks = element_line(linewidth = 0.8)  # Increase tick line width.
  ) +xlab("lasso R²")+ylab("GBDT R²")




p2<-ggplot(two_model, aes(x=delt_R2)) +
  geom_histogram(binwidth=0.01, fill="# Edd064",color="white") + Plothistogram.
  geom_vline(aes(xintercept=median(two_model$delt_R2, na.rm=TRUE)),  # Add the median line.
             color="grey50", linetype="dashed", size=1) +
  labs(title="Frequency Distribution", x="Variable", y="Frequency") +
  theme_light()+
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14,family = "Helvetica"),
    axis.title = element_text(size = 14,family = "Helvetica"),
    axis.text.x = element_text(family = "Helvetica") , # Rotate x-axis labels if group names are long.
    axis.ticks.length = unit(0.25, "cm"),  # Increase tick length.
    axis.ticks = element_line(linewidth = 0.8)  # Increase tick line width.
  ) +xlab("GBDT R²-lasso R²")


p3<-p1+p2

ggsave("4_manuscript/Figures/Figure_3/Figure_s3_gut.pdf", 
       p3, width = 8, height = 4)



# For GBDTLASSO.
rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")

library(tidyverse)
library(tidymass)
library(plyr)
library(microbiomedataset)


# Import the data.

GBDT_results<-readRDS("3_data_analysis/oral_microbiome/GBDT/cross_section/oral_GBDT_results")

LASSO_results<-readRDS("3_data_analysis/oral_microbiome/Lasso/cross_section/oral_lasso_results")


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

# Generate the plot.


p1<-ggplot(two_model, aes(x=lasso_R2, y=GBDT_R2)) +
  geom_point(shape=21, size=4, fill="#a1d5b9", color="white") +
  geom_abline(intercept=0, slope=1, color="grey50", linetype="dashed", size=1) +  # Add the diagonal reference line.
  scale_x_continuous(limits=c(0, 0.2)) +
  scale_y_continuous(limits=c(0, 0.2)) +
  theme_light()+
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14,family = "Helvetica"),
    axis.title = element_text(size = 14,family = "Helvetica"),
    axis.text.x = element_text(family = "Helvetica") , # Rotate x-axis labels if group names are long.
    axis.ticks.length = unit(0.25, "cm"),  # Increase tick length.
    axis.ticks = element_line(linewidth = 0.8)  # Increase tick line width.
  ) +xlab("lasso R²")+ylab("GBDT R²")




p2<-ggplot(two_model, aes(x=delt_R2)) +
  geom_histogram(binwidth=0.01, fill="# A1d5b9",color="white") + Plothistogram.
  geom_vline(aes(xintercept=median(two_model$delt_R2, na.rm=TRUE)),  # Add the median line.
             color="grey50", linetype="dashed", size=1) +
  labs(title="Frequency Distribution", x="Variable", y="Frequency") +
  theme_light()+
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14,family = "Helvetica"),
    axis.title = element_text(size = 14,family = "Helvetica"),
    axis.text.x = element_text(family = "Helvetica") , # Rotate x-axis labels if group names are long.
    axis.ticks.length = unit(0.25, "cm"),  # Increase tick length.
    axis.ticks = element_line(linewidth = 0.8)  # Increase tick line width.
  ) +xlab("GBDT R²-lasso R²")


p3<-p1+p2

ggsave("4_manuscript/Figures/Figure_3/Figure_s3_oral.pdf", 
       p3, width = 8, height = 4)




# For GBDTLASSO.
rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")

library(tidyverse)
library(tidymass)
library(plyr)
library(microbiomedataset)


# Import the data.

GBDT_results<-readRDS("3_data_analysis/nasal_microbiome/GBDT/cross_section/nasal_GBDT_results")

LASSO_results<-readRDS("3_data_analysis/nasal_microbiome/Lasso/cross_section/nasal_lasso_results")


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

# Generate the plot.


p1<-ggplot(two_model, aes(x=lasso_R2, y=GBDT_R2)) +
  geom_point(shape=21, size=4, fill="#a17db4", color="white") +
  geom_abline(intercept=0, slope=1, color="grey50", linetype="dashed", size=1) +  # Add the diagonal reference line.
  scale_x_continuous(limits=c(0, 0.2)) +
  scale_y_continuous(limits=c(0, 0.2)) +
  theme_light()+
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14,family = "Helvetica"),
    axis.title = element_text(size = 14,family = "Helvetica"),
    axis.text.x = element_text(family = "Helvetica") , # Rotate x-axis labels if group names are long.
    axis.ticks.length = unit(0.25, "cm"),  # Increase tick length.
    axis.ticks = element_line(linewidth = 0.8)  # Increase tick line width.
  ) +xlab("lasso R²")+ylab("GBDT R²")




p2<-ggplot(two_model, aes(x=delt_R2)) +
  geom_histogram(binwidth=0.01, fill="# A17db4",color="white") + Plothistogram.
  geom_vline(aes(xintercept=median(two_model$delt_R2, na.rm=TRUE)),  # Add the median line.
             color="grey50", linetype="dashed", size=1) +
  labs(title="Frequency Distribution", x="Variable", y="Frequency") +
  theme_light()+
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14,family = "Helvetica"),
    axis.title = element_text(size = 14,family = "Helvetica"),
    axis.text.x = element_text(family = "Helvetica") , # Rotate x-axis labels if group names are long.
    axis.ticks.length = unit(0.25, "cm"),  # Increase tick length.
    axis.ticks = element_line(linewidth = 0.8)  # Increase tick line width.
  ) +xlab("GBDT R²-lasso R²")

p3<-p1+p2

ggsave("4_manuscript/Figures/Figure_3/Figure_s3_nasal.pdf", 
       p3, width = 8, height = 4)


# For GBDTLASSO.
rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")

library(tidyverse)
library(tidymass)
library(plyr)
library(microbiomedataset)


# Import the data.

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

# Generate the plot.


p1<-ggplot(two_model, aes(x=lasso_R2, y=GBDT_R2)) +
  geom_point(shape=21, size=4, fill="#f2ccac", color="white") +
  geom_abline(intercept=0, slope=1, color="grey50", linetype="dashed", size=1) +  # Add the diagonal reference line.
  scale_x_continuous(limits=c(0, 0.2)) +
  scale_y_continuous(limits=c(0, 0.2)) +
  theme_light()+
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14,family = "Helvetica"),
    axis.title = element_text(size = 14,family = "Helvetica"),
    axis.text.x = element_text(family = "Helvetica") , # Rotate x-axis labels if group names are long.
    axis.ticks.length = unit(0.25, "cm"),  # Increase tick length.
    axis.ticks = element_line(linewidth = 0.8)  # Increase tick line width.
  ) +xlab("lasso R²")+ylab("GBDT R²")




p2<-ggplot(two_model, aes(x=delt_R2)) +
  geom_histogram(binwidth=0.01, fill="# F2ccac",color="white") + Plothistogram.
  geom_vline(aes(xintercept=median(two_model$delt_R2, na.rm=TRUE)),  # Add the median line.
             color="grey50", linetype="dashed", size=1) +
  labs(title="Frequency Distribution", x="Variable", y="Frequency") +
  theme_light()+
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14,family = "Helvetica"),
    axis.title = element_text(size = 14,family = "Helvetica"),
    axis.text.x = element_text(family = "Helvetica") , # Rotate x-axis labels if group names are long.
    axis.ticks.length = unit(0.25, "cm"),  # Increase tick length.
    axis.ticks = element_line(linewidth = 0.8)  # Increase tick line width.
  ) +xlab("GBDT R²-lasso R²")

p3<-p1+p2

ggsave("4_manuscript/Figures/Figure_3/Figure_s3_skin.pdf", 
       p3, width = 8, height = 4)