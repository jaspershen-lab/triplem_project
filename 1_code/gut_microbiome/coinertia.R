rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")
source("1_code/mantel_Procrustes_code.R")
library(tidyverse)
library(tidymass)
library(readxl)


###load("data)
load("3_data_analysis/gut_microbiome/data_preparation/object_cross_section")

gut_object<-object_cross_section

load("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section")

metabolomics_object<-object_cross_section

metabolomics_class<-read_excel("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx")

dir.create("3_data_analysis/gut_microbiome/MMC",recursive = TRUE)



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


########

microbiome_data<-data.frame(gut_temp_object@expression_data,check.names = FALSE)
metabolomics_data<-data.frame(metabolomics_temp_object@expression_data,check.names = FALSE)

metabolomics_data<-metabolomics_data[,colnames(microbiome_data)]

metabolomics_Microbial<-metabolomics_class%>%filter(HMDB.Source.Microbial=="TRUE")
metabolomics_data<-metabolomics_data[metabolomics_Microbial$variable_id,]

library(ade4)

# Translated comment.
# Translated comment.
# Translated comment.

# Translated comment.
# Translated comment.
microbiome_scaled <- t(microbiome_data)
metabolome_scaled <- t(metabolomics_data)

# Translated comment.
pca_microbiome <- dudi.pca(data.frame(microbiome_scaled), scannf = FALSE, nf = 5)
pca_metabolome <- dudi.pca(data.frame(metabolome_scaled), scannf = FALSE, nf = 5)

# Translated comment.
coia <- coinertia(pca_microbiome, pca_metabolome, scannf = FALSE, nf = 2)

# Translated comment.
# Translated comment.
coia$RV

# Translated comment.
distances <- sqrt(rowSums((coia$mX - coia$mY)^2))

# Translated comment.
# Translated comment.
str(coia)

# Translated comment.
coia$eig

# Translated comment.
percent_var <- (coia$eig/sum(coia$eig))*100
print(percent_var)

# Translated comment.
# Translated comment.
micro_scores <- coia$li  # translated comment
micro_loadings <- coia$c1 # translated comment

# Translated comment.
metab_scores <- coia$li # translated comment
metab_loadings <- coia$l1 # translated comment

# Translated comment.
# Translated comment.
head(micro_loadings)
# Translated comment.
head(metab_loadings)




a<-cbind(gut_temp_object@sample_info,distances)

ggplot(a, aes(x=adjusted_age, y= distances)) +
  geom_point(shape=21,size=4,fill="#A1D0C7",color="white") +
  geom_smooth(method="lm",colour = "grey50") +theme_light() +stat_cor(method = "pearson")+theme(legend.position="none", # translated comment
                                                                                                axis.text.x=element_text(colour="black",size=14), # translated comment
                                                                                                axis.text.y=element_text(size=14,face="plain"), # translated comment
                                                                                                axis.title.y=element_text(size = 14,face="plain"), # translated comment
                                                                                                axis.title.x=element_text(size = 14,face="plain"), # translated comment
                                                                                                plot.title = element_text(size=15,face="bold",hjust = 0.5))+ylab("MMC")+xlab("Age")+xlim(c(33,70))



### Translated comment.


micro_loadings<-subset(micro_loadings,CS1<(-0.1))
metab_loadings<-subset(metab_loadings,RS1<(-0.1))


micro_loadings <- merge(micro_loadings, data.frame(gut_temp_object@variable_info), by = "row.names")

row.names(metabolomics_Microbial)<-metabolomics_Microbial$variable_id


metab_loadings$variable_id<-row.names(metab_loadings)
metab_loadings <- merge(metab_loadings, metabolomics_Microbial, by = "variable_id")

metab_loadings <- arrange(metab_loadings, RS1,decreasing=TRUE)
metab_loadings$HMDB.Name<-factor(metab_loadings$HMDB.Name,levels = unique(metab_loadings$HMDB.Name))
p1<-ggplot(metab_loadings, aes(x = HMDB.Name, y = abs(RS1))) +
  geom_segment(aes(xend = HMDB.Name, yend = 0), color = "gray50") +
  geom_point(size = 3, color = "steelblue") +
  coord_flip() + # translated comment
  theme_bw() +
  labs(
    title = "Top Estimate Metabolites ",
    x = "Metabolites",
    y = "Loadings"
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(face = "bold")
  )


micro_loadings <- arrange(micro_loadings, CS1,decreasing=TRUE)
micro_loadings$Genus<-factor(micro_loadings$Genus,levels = unique(micro_loadings$Genus))

p2<-ggplot(micro_loadings, aes(x = Genus, y = abs(CS1))) +
  geom_segment(aes(xend = Genus, yend = 0), color = "gray50") +
  geom_point(size = 3, color = "steelblue") +
  coord_flip() + # translated comment
  theme_bw() +
  labs(
    title = "Top Estimate Microbiome ",
    x = "Genus",
    y = "Loadings"
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(face = "bold")
  )



