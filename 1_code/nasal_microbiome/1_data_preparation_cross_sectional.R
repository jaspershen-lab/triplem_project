setwd(r4projects::get_project_wd())
rm(list = ls())


library(tidyverse)
library(tidymass)
load("2_data/nasal_microbiome/data_preparation/object")
setwd("3_data_analysis/nasal_microbiome/data_preparation/")

massdataset::export_mass_dataset(object = object, 
                                 file_type = "xlsx")

dim(object)

object <-
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(!is.na(adjusted_age))

####only remain the health vist
object <-
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(CL4 == "Healthy")

###for each participant, calculate the median value
object_cross_section <-
  massdataset::summarise_samples(object = object,
                                 what = "mean_intensity",
                                 group_by = "subject_id")

dim(object_cross_section)

save(object_cross_section, file = "object_cross_section")