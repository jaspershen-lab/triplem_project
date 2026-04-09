rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")
source("1_code/mantel_Procrustes_code.R")
library(tidyverse)
library(tidymass)
library(readxl)


###load("data)
load("3_data_analysis/nasal_microbiome/data_preparation/object_cross_section")

nasal_object<-object_cross_section

load("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section")

metabolomics_object<-object_cross_section

metabolomics_class<-read_excel("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx")

dir.create("3_data_analysis/nasal_microbiome/MMC",recursive = TRUE)



####only remain the genus level
library(microbiomedataset)

nasal_object <-
  nasal_object %>%
  summarize_variables(what = "sum_intensity", group_by = "Genus")

##only remain the genus at least 10% samples
dim(nasal_object)

non_zero_per <-
  apply(nasal_object, 1, function(x) {
    sum(x != 0) / ncol(nasal_object)
  })

idx <-
  which(non_zero_per > 0.1)

nasal_object <-
  nasal_object[idx, ]


nasal_object <-
  nasal_object %>%
  transform2relative_intensity()



##
##adjust BMI, sex, and IRIS, ethnicity
library(tidyverse)
library(ggpubr)
library(rstatix)

nasal_expression_data <-
  extract_expression_data(nasal_object) %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

library(plyr)

nasal_sample_info <-
  nasal_object@sample_info

#######adjust BMI, sex, and IRIS, ethnicity
nasal_expression_data <-
  lm_adjust(expression_data = nasal_expression_data,
            sample_info = nasal_sample_info,
            threads = 3)

nasal_temp_object <- nasal_object
nasal_temp_object@expression_data <- nasal_expression_data



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


########

microbiome_data<-data.frame(nasal_temp_object@expression_data,check.names = FALSE)
metabolomics_data<-data.frame(metabolomics_temp_object@expression_data,check.names = FALSE)

metabolomics_data<-metabolomics_data[,colnames(microbiome_data)]

metabolomics_Microbial<-metabolomics_class%>%filter(HMDB.Source.Microbial=="TRUE")
metabolomics_data<-metabolomics_data[metabolomics_Microbial$variable_id,]

library(ade4)

# data:.
# Microbiome_data: samples,microbiome.
# Metabolome_data: samples,metabolites.

# 1. data preprocessing.
# data,Ensurerow names.
microbiome_scaled <- t(microbiome_data)
metabolome_scaled <- t(metabolomics_data)

# 2. for datasetPCA.
pca_microbiome <- dudi.pca(data.frame(microbiome_scaled), scannf = FALSE, nf = 5)
pca_metabolome <- dudi.pca(data.frame(metabolome_scaled), scannf = FALSE, nf = 5)

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







a<-cbind(nasal_temp_object@sample_info,distances)




ggplot(a, aes(x=adjusted_age, y= distances)) +
  geom_point(shape=21,size=4,fill="#A1D0C7",color="white") +
  geom_smooth(method="lm",colour = "grey50") +theme_light() +stat_cor(method = "pearson")+theme(legend.position="none", # Hide the legend.
                                                                                                axis.text.x=element_text(colour="black",size=14), # Set x-axis tick label text properties.
                                                                                                axis.text.y=element_text(size=14,face="plain"), # Set x-axis tick label text properties.
                                                                                                axis.title.y=element_text(size = 14,face="plain"), # Set y-axis title text properties.
                                                                                                axis.title.x=element_text(size = 14,face="plain"), # Set x-axis title text properties.
                                                                                                plot.title = element_text(size=15,face="bold",hjust = 0.5))+xlab("Age")