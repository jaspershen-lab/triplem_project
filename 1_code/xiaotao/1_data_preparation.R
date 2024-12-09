no_source()
setwd(r4projects::get_project_wd())
library(tidyverse)
rm(list = ls())

## load data
load("3_data_analysis/plasma_metabolomics/data_preparation/peak/object")

setwd("3_data_analysis/plasma_metabolomics/data_preparation/peak")

dim(object)

colnames(object)

library(tidymass)

extract_process_info(object)

### remove QC and Blank samples from the object
object@sample_info$class

dim(object)

object <-
    object %>%
    dplyr::filter(class != "QC" & class != "Blank")

dim(object)

colnames(object@sample_info)

table(object@sample_info$CL4)

### only remain the healthy samples for each person
object %>%
    dplyr::filter(CL4 == "Healthy") -> object

length(unique(object@sample_info$subject_id))

dim(object)

object@sample_info %>% 
dplyr::count(subject_id) %>% 
dplyr::arrange(desc(n))


###get the cross sectional dataset

library(microbiomedataset)

object_cross_section <-
object %>% 
microbiomedataset::convert2microbiome_dataset() %>% 
microbiomedataset::summarise_samples(
what = "mean_intensity",
group_by = "subject_id") %>% 
microbiomedataset::convert2mass_dataset()

dim(object_cross_section)

save(object_cross_section, file = "object_cross_section.rda")





