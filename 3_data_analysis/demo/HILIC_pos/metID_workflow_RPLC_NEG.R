###
library(tidyverse)

library(massdataset)
library(metid)

rm(list = ls())

setwd("3_data_analysis/demo/RPLC_neg/")####please change this path

####Positive mode
variable_info <-
  readr::read_csv("metID_nRPLC.csv") %>%
  dplyr::rename(variable_id = name)

expression_data <-
  matrix(1, nrow = nrow(variable_info),
         ncol = 1) %>%
  as.data.frame()

rownames(expression_data) <-
  variable_info$variable_id

colnames(expression_data) <- "sample1"

sample_id = colnames(expression_data)

sample_info <-
  data.frame(sample_id = colnames(expression_data)) %>%
  dplyr::mutate(class = "Subject")

dim(expression_data)
dim(sample_info)
dim(variable_info)

colnames(expression_data) == sample_info$sample_id
rownames(expression_data) == variable_info$variable_id

object <-
  create_mass_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info
  )

###Add MS2 data
object <-
  mutate_ms2(
    object = object,
    column = "rp",
    polarity = "negative",
    path = "."
  )

#####annotation
load("hmdb_ms2.rda")
load("massbank_ms2.rda")
load("mona_ms2.rda")
load("mpsnyder_hilic_ms2.rda")
load("mpsnyder_rplc_ms2.rda")
load("nist_ms2.rda")

###mpsnyder_rplc
object <-
  metid::annotate_metabolites_mass_dataset(
    object = object,
    polarity = "negative",
    rt.match.tol = 30,
    column = "rp",
    candidate.num = 3,
    database = mpsnyder_rplc_ms2,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )

##MassBank
object <-
  metid::annotate_metabolites_mass_dataset(
    object = object,
    polarity = "negative",
    rt.match.tol = 30,
    column = "rp",
    candidate.num = 3,
    database = massbank_ms2,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )

###MONA
object <-
  metid::annotate_metabolites_mass_dataset(
    object = object,
    polarity = "negative",
    rt.match.tol = 30,
    column = "rp",
    candidate.num = 3,
    database = mona_ms2,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )

###NIST
object <-
  metid::annotate_metabolites_mass_dataset(
    object = object,
    polarity = "negative",
    rt.match.tol = 30,
    column = "rp",
    candidate.num = 3,
    database = nist_ms2,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )


###HMDB
object <-
  metid::annotate_metabolites_mass_dataset(
    object = object,
    polarity = "negative",
    rt.match.tol = 30,
    column = "rp",
    candidate.num = 3,
    database = hmdb_ms2,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )


save(object, file = "object")

variable_info <-
  extract_variable_info(object)

write.csv(variable_info, file = "annotation_result.csv", row.names = FALSE)
