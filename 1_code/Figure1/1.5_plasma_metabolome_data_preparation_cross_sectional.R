
###
setwd(r4projects::get_project_wd())
rm(list = ls())
# source("1-code/100-tools.R")
library(tidyverse)
library(tidymass)

load("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object")

setwd("3_data_analysis/plasma_metabolomics/data_preparation/metabolite")

massdataset::export_mass_dataset(object = object, 
                                 file_type = "xlsx")

object

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

save(object_cross_section, file = "object_cross_section")