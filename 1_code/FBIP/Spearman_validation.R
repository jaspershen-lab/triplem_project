rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")


## read MetaCardis data
library(readxl)

FBIP_metagenomic<-read.table("2_data/FBIP-main/otu_table/merged_abundance_table_species.txt",header = TRUE)

FBIP_metagenomic<-aggregate(FBIP_metagenomic[,8:669],list(FBIP_metagenomic$Genus),sum)

FBIP_metagenomic$Group.1<-gsub("g__","",FBIP_metagenomic$Group.1)


row.names(FBIP_metagenomic)<-FBIP_metagenomic$Group.1

FBIP_metagenomic<-FBIP_metagenomic[,-1]

FBIP_metagenomic<-data.frame(t(FBIP_metagenomic),check.names = FALSE)



FBIP_metabolome<-read.table("2_data/FBIP-main/metabolites data/plasma metabolites data.txt",header = TRUE,sep = "\t")

FBIP_metadata<-read.table("2_data/FBIP-main/metadata/FBIP_metadata.txt",header = TRUE,sep = "\t")
FBIP_metadata<-subset(FBIP_metadata,time_num%in%c("W0","W4","W16"))

FBIP_metagenomic<-FBIP_metagenomic[FBIP_metadata$SampleID,]


FBIP_metabolome<-subset(FBIP_metabolome,Time%in%c("W0","W4","W16"))

FBIP_metabolome<-merge(FBIP_metabolome,FBIP_metadata[,c("Time","SampleID")],by.x="ST",by.y="Time")
rownames(FBIP_metabolome)<-FBIP_metabolome$SampleID
FBIP_metabolome<-FBIP_metabolome[,c(-1,-2,-3,-197)]

FBIP_metagenomic<-FBIP_metagenomic[rownames(FBIP_metabolome),]

FBIP_metabolome_annotation<-read.table("2_data/FBIP-main/metabolites data/plasma_metabolites_name_group.txt",header = TRUE,sep = "\t")

# Translated comment.

FBIP_metabolome_annotation<-subset(FBIP_metabolome_annotation,HMDB!=c(""))

FBIP_metabolome<-FBIP_metabolome[,FBIP_metabolome_annotation$ID]

colnames(FBIP_metabolome)<-FBIP_metabolome_annotation$HMDB


# Translated comment.


library(tidyverse)
library(tidymass)
library(plyr)
library(microbiomedataset)
###load("data)
load("3_data_analysis/gut_microbiome/data_preparation/object_cross_section")

gut_object<-object_cross_section

load("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section")

metabolomics_object<-object_cross_section

metabolite_annotation<-read_excel("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx")




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



## Translated comment.
com_tax<-intersect(gut_temp_object@variable_info$Genus,colnames(FBIP_metagenomic))


microbiome_data_ipop<-gut_temp_object@expression_data
microbiome_tax_ipop<-gut_temp_object@variable_info

rownames(microbiome_data_ipop)<-microbiome_tax_ipop$Genus

microbiome_data_ipop<-microbiome_data_ipop[com_tax,]

FBIP_metagenomic<-data.frame(t(FBIP_metagenomic),check.names = FALSE)

FBIP_metagenomic<-FBIP_metagenomic[com_tax,]

FBIP_metagenomic <-
  FBIP_metagenomic%>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()


### Translated comment.

metabolome_data_ipop<-metabolomics_temp_object@expression_data
metabolome_ann_ipop<-data.frame(metabolite_annotation)

com_metabolites<-intersect(metabolome_ann_ipop$HMDB,FBIP_metabolome_annotation$HMDB)



# 

metabolome_ann_ipop<-subset(metabolome_ann_ipop,HMDB%in%com_metabolites)

metabolome_data_ipop<-metabolome_data_ipop[metabolome_ann_ipop$variable_id,]


FBIP_metabolome_annotation<-subset(FBIP_metabolome_annotation,HMDB%in%com_metabolites)

FBIP_metabolome<-FBIP_metabolome[,com_metabolites]

common_samples <- intersect(colnames(metabolome_data_ipop), colnames(microbiome_data_ipop))
metabolome_data <- metabolome_data_ipop[, common_samples]
microbiome_data <- microbiome_data_ipop[, common_samples]

# Translated comment.
metabolome_t <- t(metabolome_data)
microbiome_t <- t(microbiome_data)

# Translated comment.
correlation_results <- data.frame()

# Translated comment.
for (i in 1:ncol(metabolome_t)) {
  metabolite_name <- colnames(metabolome_t)[i]
  metabolite_data <- metabolome_t[, i]
  
  for (j in 1:ncol(microbiome_t)) {
    microbe_name <- colnames(microbiome_t)[j]
    microbe_data <- microbiome_t[, j]
    
    # Translated comment.
    valid_indices <- complete.cases(metabolite_data, microbe_data)
    
    if (sum(valid_indices) > 5) {  # translated comment
      # Translated comment.
      cor_test <- cor.test(metabolite_data[valid_indices], 
                           microbe_data[valid_indices], 
                           method = "spearman")
      
      # Translated comment.
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

# Translated comment.
correlation_results$FDR <- p.adjust(correlation_results$P_value, method = "BH")

# Translated comment.
correlation_results$Abs_Correlation <- abs(correlation_results$Correlation)

# Translated comment.
correlation_results <- correlation_results[order(correlation_results$Abs_Correlation, decreasing = TRUE), ]

# Translated comment.
top_100_correlations <- head(correlation_results, 500)


top_100_correlations<-merge(top_100_correlations,metabolite_annotation[,c("variable_id","HMDB","HMDB.Name")],by.x="Metabolite",by.y="variable_id")

top_100_correlations$bac_meta <- paste(top_100_correlations$Microbe, top_100_correlations$HMDB, sep = "_")


top_100_correlations <- top_100_correlations %>%
  group_by(bac_meta) %>%
  dplyr::summarise(across(where(is.numeric), mean, na.rm = TRUE))






library(tidyverse)  # translated comment
library(Hmisc) 

#######FBIP
# Translated comment.
# Translated comment.
metagenomic_data <- FBIP_metagenomic
metabolome_data <- data.frame(t(FBIP_metabolome))

# Translated comment.
cat("代谢物数据维度:", dim(metagenomic_data), "\n")
cat("细菌数据维度:", dim(metabolome_data), "\n")

# Translated comment.
# Translated comment.
common_samples <- intersect(colnames(metagenomic_data), colnames(metabolome_data))

# Translated comment.
if (length(common_samples) == 0) {
  stop("没有共同的样本名在两个数据集中")
}

# Translated comment.
metagenomic_filtered <- metagenomic_data[, common_samples]
metabolome_filtered <- metabolome_data[, common_samples]

# Translated comment.
# Translated comment.
metagenomic_t <- t(metagenomic_filtered)
metabolome_t <- t(metabolome_filtered)

# Translated comment.
correlation_result <- rcorr(as.matrix(metagenomic_t), as.matrix(metabolome_t), type = "spearman")

# Translated comment.
cor_coef <- correlation_result$r
cor_pval <- correlation_result$P

# Translated comment.
# Translated comment.
cor_subset <- cor_coef[1:ncol(metagenomic_t), (ncol(metagenomic_t)+1):(ncol(metagenomic_t)+ncol(metabolome_t))]
pval_subset <- cor_pval[1:ncol(metagenomic_t), (ncol(metagenomic_t)+1):(ncol(metagenomic_t)+ncol(metabolome_t))]

# Translated comment.
correlation_table <- data.frame()

for (i in 1:nrow(cor_subset)) {
  for (j in 1:ncol(cor_subset)) {
    metabolite <- rownames(cor_subset)[i]
    bacteria <- colnames(cor_subset)[j]
    rho <- cor_subset[i, j]
    pvalue <- pval_subset[i, j]
    
    correlation_table <- rbind(correlation_table, 
                               data.frame(Metabolite = metabolite,
                                          Bacteria = bacteria,
                                          Rho = rho,
                                          Pvalue = pvalue))
  }
}

# Translated comment.
correlation_table <- correlation_table[order(correlation_table$Pvalue), ]

# Translated comment.
correlation_table$AdjustedPvalue <- p.adjust(correlation_table$Pvalue, method = "BH")

# Translated comment.
correlation_table$Significance <- ""
correlation_table$Significance[correlation_table$Pvalue < 0.05] <- "*"
correlation_table$Significance[correlation_table$Pvalue < 0.01] <- "**"
correlation_table$Significance[correlation_table$Pvalue < 0.001] <- "***"


correlation_table_FBIP<-correlation_table

correlation_table_FBIP$bac_meta <- paste(correlation_table_FBIP$Metabolite, correlation_table_FBIP$Bacteria, sep = "_")

correlation_table_FBIP<-subset(correlation_table_FBIP,bac_meta%in%top_100_correlations$bac_meta)



correlation_table_FBIP_ipop<-merge(correlation_table_FBIP,top_100_correlations,by="bac_meta")
correlation_table_FBIP_ipop <- correlation_table_FBIP_ipop[which(sign(correlation_table_FBIP_ipop$Rho) * sign(correlation_table_FBIP_ipop$Correlation) > 0), ]

correlation_table_FBIP_ipop<-subset(correlation_table_FBIP_ipop,abs(Rho)>0.1)


ggplot(correlation_table_FBIP_ipop, aes(x=correlation_table_FBIP_ipop$Rho, y= correlation_table_FBIP_ipop$Correlation)) +
  geom_point(shape=21,size=4,fill="#A1D0C7",color="white") +
  geom_smooth(method="lm",colour = "grey50") +theme_light() +stat_cor(method = "spearman")+theme(legend.position="none", # translated comment
                                                                                                 axis.text.x=element_text(colour="black",size=14), # translated comment
                                                                                                 axis.text.y=element_text(size=14,face="plain"), # translated comment
                                                                                                 axis.title.y=element_text(size = 14,face="plain"), # translated comment
                                                                                                 axis.title.x=element_text(size = 14,face="plain"), # translated comment
                                                                                                 plot.title = element_text(size=15,face="bold",hjust = 0.5))+xlab("iPOP_Rho")+ylab("FBIP_Rho")