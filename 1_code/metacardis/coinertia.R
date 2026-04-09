rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")


## read MetaCardis data
library(readxl)

MetaCardis_metagenomic<-data.frame(read_excel("2_data/MetaCardis cohort/European MetaCardis cohort _metagenomic.xlsx"))
row.names(MetaCardis_metagenomic)<-MetaCardis_metagenomic$ID

MetaCardis_metagenomic<-MetaCardis_metagenomic[,-1]


MetaCardis_tax<-data.frame(read_excel("2_data/MetaCardis cohort/taxonomy.xlsx"))

MetaCardis_metabolome<-data.frame(read_excel("2_data/MetaCardis cohort/European MetaCardis cohort_metabolome.xlsx"))

row.names(MetaCardis_metabolome)<-MetaCardis_metabolome$ID

MetaCardis_metabolome<-MetaCardis_metabolome[,-1]

MetaCardis_metabome_annotation<-data.frame(read_excel("2_data/MetaCardis cohort/metabolites_annotation.xlsx"))

MetaCardis_metadata<-data.frame(read_excel("2_data/MetaCardis cohort/metadata.xlsx"))
metabolite_annotation<-read_excel("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx")
#metabolite_annotation<-subset(metabolite_annotation,metabolite_annotation$HMDB.Source.Microbial=="TRUE")

MetaCardis_metabome_annotation_microbial<-subset(MetaCardis_metabome_annotation,HMDB%in%metabolite_annotation$HMDB)

# Remove duplicate participants.

MetaCardis_metadata<-subset(MetaCardis_metadata,Status%in%c("MMC372","HC275","UMCC222","IHD372"))

MetaCardis_metagenomic<-MetaCardis_metagenomic[MetaCardis_metadata$ID,]

MetaCardis_metabolome<-MetaCardis_metabolome[MetaCardis_metadata$ID,]


# Keep participants with both metabolomics and microbiome data.
MetaCardis_metagenomic <- MetaCardis_metagenomic[!apply(MetaCardis_metagenomic == "NA", 1, all), ]
MetaCardis_metabolome <- MetaCardis_metabolome[!apply(MetaCardis_metabolome == "NA", 1, all), ]

com_sample<-intersect(row.names(MetaCardis_metagenomic),row.names(MetaCardis_metabolome))

MetaCardis_metagenomic<-MetaCardis_metagenomic[com_sample,]
MetaCardis_metabolome<-MetaCardis_metabolome[com_sample,]
MetaCardis_metadata<-subset(MetaCardis_metadata,ID%in%com_sample)
MetaCardis_metabolome <- type.convert(MetaCardis_metabolome, as.is = TRUE)

MetaCardis_metabolome<-na.omit(MetaCardis_metabolome)
# Aggregate MetaCardis_metagenomic to the genus level.




MetaCardis_metagenomic <- type.convert(MetaCardis_metagenomic, as.is = TRUE)

MetaCardis_metagenomic<-data.frame(t(MetaCardis_metagenomic))

MetaCardis_metagenomic<-cbind(MetaCardis_metagenomic,MetaCardis_tax[,8:4])


MetaCardis_metagenomic<-aggregate(MetaCardis_metagenomic[,1:1146],list(MetaCardis_metagenomic$genus),sum)

rownames(MetaCardis_metagenomic)<-  MetaCardis_metagenomic$Group.1

MetaCardis_metagenomic<-MetaCardis_metagenomic[,-1]

MetaCardis_metagenomic<-t(MetaCardis_metagenomic)
MetaCardis_metagenomic<-MetaCardis_metagenomic/apply(MetaCardis_metagenomic,1,sum)*100

MetaCardis_metagenomic<-MetaCardis_metagenomic[rownames(MetaCardis_metabolome),]
MetaCardis_metadata<-MetaCardis_metadata%>%filter(ID%in%rownames(MetaCardis_metabolome))



library(ade4)




# Read the data.
metabolome <-MetaCardis_metagenomic
metadata <- MetaCardis_metadata

# Ensuremetabolite data.
metabolome <- as.data.frame(apply(metabolome, 2, as.numeric))

# Extractdata.
age <- metadata$Age..years.

# Createdata framecorrelationresults.
correlation_results <- data.frame(
  Metabolite = colnames(metabolome),
  Correlation = NA,
  P_value = NA
)

# Calculatemetabolitesand correlation.
for(i in 1:ncol(metabolome)) {
  # Usecomplete.casesEnsureUsedatafor.
  complete_data <- complete.cases(metabolome[,i], age)
  
  # Calculatecorrelation.
  cor_test <- cor.test(metabolome[complete_data,i], 
                       age[complete_data], 
                       method = "pearson",  # UseSpearmancorrelation coefficient.
                       exact = FALSE)        # For samplesUseCalculate.
  
  # results.
  correlation_results$Correlation[i] <- cor_test$estimate
  correlation_results$P_value[i] <- cor_test$p.value
}

# AddFDRp-values.
correlation_results$FDR <- p.adjust(correlation_results$P_value, method = "BH")

# correlation coefficientfor sort.
correlation_results <- correlation_results[order(abs(correlation_results$Correlation), decreasing = TRUE),]

# Addsignificance.
correlation_results$Significance <- ifelse(correlation_results$FDR < 0.05, 
                                           "Significant", 
                                           "Not Significant")


correlation_results_MGS<-subset(correlation_results,Correlation>0&Significance=="Significant")



metabolome <-MetaCardis_metabolome
metadata <- MetaCardis_metadata

# Ensuremetabolite data.
metabolome <- as.data.frame(apply(metabolome, 2, as.numeric))

# Extractdata.
age <- metadata$Age..years.

# Createdata framecorrelationresults.
correlation_results <- data.frame(
  Metabolite = colnames(metabolome),
  Correlation = NA,
  P_value = NA
)

# Calculatemetabolitesand correlation.
for(i in 1:ncol(metabolome)) {
  # Usecomplete.casesEnsureUsedatafor.
  complete_data <- complete.cases(metabolome[,i], age)
  
  # Calculatecorrelation.
  cor_test <- cor.test(metabolome[complete_data,i], 
                       age[complete_data], 
                       method = "pearson",  # UseSpearmancorrelation coefficient.
                       exact = FALSE)        # For samplesUseCalculate.
  
  # results.
  correlation_results$Correlation[i] <- cor_test$estimate
  correlation_results$P_value[i] <- cor_test$p.value
}

# AddFDRp-values.
correlation_results$FDR <- p.adjust(correlation_results$P_value, method = "BH")

# correlation coefficientfor sort.
correlation_results <- correlation_results[order(abs(correlation_results$Correlation), decreasing = TRUE),]

# Addsignificance.
correlation_results$Significance <- ifelse(correlation_results$FDR < 0.05, 
                                           "Significant", 
                                           "Not Significant")

correlation_results_Met<-subset(correlation_results,Correlation>0.1&Significance=="Significant")




# 2. for datasetPCA.
pca_microbiome <- dudi.pca(data.frame(MetaCardis_metagenomic), scannf = FALSE, nf = 5)
pca_metabolome <- dudi.pca(data.frame(MetaCardis_metabolome[,correlation_results_Met$Metabolite]), scannf = FALSE, nf = 5)

# 3. CoinertiaAnalyze.
coia <- coinertia(pca_microbiome, pca_metabolome, scannf = FALSE, nf = 2)

# 4. Viewresults.
# RV(correlation).
coia$RV

# Calculatesamples.
distances <- sqrt(rowSums((coia$mX - coia$mY)^2))

# resultsCreatedata frame.
# Viewcoinertiafor .
str(coia)

# Viewfeatures(eigenvalues).
coia$eig

# View.
percent_var <- (coia$eig/sum(coia$eig))*100
print(percent_var)

# Getdatacoordinates.
# Microbiome datacoordinates.
micro_scores <- coia$li  # Samplesmicrobiomecoordinates.
micro_loadings <- coia$c1 # Microbiomevariablescontribution.

# Metabolite datacoordinates.
metab_scores <- coia$li # Samplesmetabolitescoordinates.
metab_loadings <- coia$l1 # Metabolitesvariablescontribution.

# Viewvariablescontribution.
# Microbiomevariablescontribution.
head(micro_loadings)
# Metabolitesvariablescontribution.
head(metab_loadings)





a<-cbind(MetaCardis_metadata,distances)
a<- type.convert(a, as.is = TRUE)

a<-subset(a,a$Status=="UMCC222")
ggplot(a, aes(x=Age..years., y= distances)) +
  geom_point(shape=21,size=4,fill="#A1D0C7",color="white") +
  geom_smooth(method="lm",colour = "grey50") +theme_light() +stat_cor(method = "pearson")+theme(legend.position="none", # Hide the legend.
                                                                                                axis.text.x=element_text(colour="black",size=14), # Set x-axis tick label text properties.
                                                                                                axis.text.y=element_text(size=14,face="plain"), # Set x-axis tick label text properties.
                                                                                                axis.title.y=element_text(size = 14,face="plain"), # Set y-axis title text properties.
                                                                                                axis.title.x=element_text(size = 14,face="plain"), # Set x-axis title text properties.
                                                                                                plot.title = element_text(size=15,face="bold",hjust = 0.5))+xlab("Age")+xlim(c(27,70))



a<-subset(a,Age..years.>27&Age..years.<70)

