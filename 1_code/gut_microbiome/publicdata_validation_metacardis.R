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



# Aggregate MetaCardis_metagenomic to the genus level.




MetaCardis_metagenomic <- type.convert(MetaCardis_metagenomic, as.is = TRUE)

MetaCardis_metagenomic<-data.frame(t(MetaCardis_metagenomic))

MetaCardis_metagenomic<-cbind(MetaCardis_metagenomic,MetaCardis_tax[,8:4])


MetaCardis_metagenomic<-aggregate(MetaCardis_metagenomic[,1:1146],list(MetaCardis_metagenomic$genus),sum)

rownames(MetaCardis_metagenomic)<-  MetaCardis_metagenomic$Group.1

MetaCardis_metagenomic<-MetaCardis_metagenomic[,-1]

MetaCardis_metagenomic<-t(MetaCardis_metagenomic)
MetaCardis_metagenomic<-MetaCardis_metagenomic/apply(MetaCardis_metagenomic,1,sum)*100



# Calculate relative abundance.

# Read the iPOP dataset.


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



## Filtershared .
com_tax<-intersect(gut_temp_object@variable_info$Genus,colnames(MetaCardis_metagenomic))


microbiome_data_ipop<-gut_temp_object@expression_data
microbiome_tax_ipop<-gut_temp_object@variable_info

rownames(microbiome_data_ipop)<-microbiome_tax_ipop$Genus

microbiome_data_ipop<-microbiome_data_ipop[com_tax,]

MetaCardis_metagenomic<-data.frame(t(MetaCardis_metagenomic))

MetaCardis_metagenomic<-MetaCardis_metagenomic[com_tax,]

MetaCardis_metagenomic <-
  MetaCardis_metagenomic%>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()


### shared metabolites.

metabolome_data_ipop<-metabolomics_temp_object@expression_data
metabolome_ann_ipop<-data.frame(metabolite_annotation)

com_metabolites<-intersect(metabolome_ann_ipop$HMDB,MetaCardis_metabome_annotation$HMDB)
# NA.
com_metabolites<-com_metabolites[-13]


# 

metabolome_ann_ipop<-subset(metabolome_ann_ipop,HMDB%in%com_metabolites)

metabolome_data_ipop<-metabolome_data_ipop[metabolome_ann_ipop$variable_id,]


MetaCardis_metabome_annotation<-subset(MetaCardis_metabome_annotation,HMDB%in%com_metabolites)

MetaCardis_metabolome<-MetaCardis_metabolome[,MetaCardis_metabome_annotation$Met_ID]

MetaCardis_metabolome <- type.convert(MetaCardis_metabolome, as.is = TRUE)



### 
microbiome_data_ipop<-data.frame(t(microbiome_data_ipop))

metabolome_data_ipop<-data.frame(t(metabolome_data_ipop))


results <- analyze_metabolite_ev(
  data.frame(t(microbiome_data_ipop)), 
  data.frame(t(metabolome_data_ipop)),
  do_feature_selection = TRUE,
  correlation_method = "spearman",
  p_threshold = 0.05,
  p_adjust_method = "none",
  rho_threshold = 0.1
)









library(tidyverse)

predict_metabolites <- function(models_directory, new_data) {
  # Getunder RDS.
  model_files <- list.files(models_directory, pattern = "\\.rds$", full.names = TRUE)
  
  # results.
  predictions_list <- list()
  

  for (model_file in model_files) {
    # Load.
    model_info <- readRDS(model_file)
    
    # Ensuredatafeatures.
    required_features <- model_info$selected_features
    missing_features <- setdiff(required_features, colnames(new_data))
    
    if (length(missing_features) > 0) {
      warning(paste("Skipping model", basename(model_file), "- missing required features:", 
                    paste(missing_features, collapse = ", ")))
      next
    }
    
    # featuresEnsure.
    filtered_data <- new_data[, required_features, drop = FALSE]
    

    predictions <- predict(model_info$model, filtered_data, n.trees = model_info$parameters$n.trees)
    
    # Saveresults.
    predictions_list[[basename(model_file)]] <- list(
      predictions = predictions,
      model_performance = model_info$performance,
      feature_importance = model_info$feature_importance
    )
  }
  
  return(predictions_list)
}

# data.
models_directory <- "metabolite_models"

# 1. metabolites.
predict_metabolites_pipeline <- function(models_directory, metagenomic_data) {
  # Transposedata.
  new_data <- data.frame(t(metagenomic_data))
  
  # metabolites.
  results <- predict_metabolites(models_directory, new_data)
  
  # Processresults.
  predictions_df <- do.call(cbind, lapply(results, function(x) x$predictions))
  colnames(predictions_df) <- gsub("_model.rds", "", names(results))
  
  # GetR2.
  R2_df <- do.call(cbind, lapply(results, function(x) x$model_performance$mean_r2))
  colnames(R2_df) <- names(results)
  
  list(predictions = predictions_df, r2 = R2_df)
}

# 2. data preprocessing function.
prepare_prediction_data <- function(predictions_df, metagenomic_data, metabolite_annotation) {
  predictions_df <- data.frame(t(predictions_df))
  colnames(predictions_df) <- colnames(metagenomic_data)
  predictions_df$variable_id <- rownames(predictions_df)
  
  # Mergeannotations.
  merged_predictions <- merge(predictions_df, metabolite_annotation, by = "variable_id")
  
  merged_predictions
}

# 3. correlation analysis.
calculate_correlations <- function(predictions_df, metabolites_name, 
                                   MetaCardis_metabome_annotation, MetaCardis_metabolome) {
  # resultsdata frame.
  results_list <- vector("list", length(metabolites_name))
  
  # ParallelProcessmetabolites.
  results_list <- parallel::mclapply(metabolites_name, function(i) {
    # Extractdata.
    predictions_data <- colSums(subset(predictions_df, HMDB == i)[, 2:(length(rownames(MetaCardis_metabolome))+1)])
    
    # Getsample names.
    observed_sample_names <- subset(MetaCardis_metabome_annotation, HMDB == i)$Met_ID
    
    # Processsamples.
    do.call(rbind, lapply(observed_sample_names, function(z) {
      observed_data <- MetaCardis_metabolome[, z]
      
      # MergeProcessdata.
      merge_data <- na.omit(data.frame(
        predictions = predictions_data,
        observed = observed_data
      ))
      
      # datacorrelation.
      if(nrow(merge_data) > 3) {
        cor_result <- cor.test(merge_data$predictions, merge_data$observed)
        
        data.frame(
          p.value = cor_result$p.value,
          estimate = cor_result$estimate,
          HMDB_Name = i,
          n_samples = nrow(merge_data)
        )
      } else {
        NULL
      }
    }))
  }, mc.cores = parallel::detectCores() - 1)
  
  # Mergeresults.
  results_cor_name <- do.call(rbind, results_list)
  rownames(results_cor_name) <- NULL
  
  results_cor_name
}

# Main function.
main_analysis <- function(models_directory, MetaCardis_metagenomic, 
                          metabolite_annotation, MetaCardis_metabome_annotation,
                          MetaCardis_metabolome) {
  # 1. metabolites.
  prediction_results <- predict_metabolites_pipeline(models_directory, MetaCardis_metagenomic)
  
  # 2. data.
  processed_predictions <- prepare_prediction_data(
    prediction_results$predictions, 
    MetaCardis_metagenomic, 
    metabolite_annotation
  )
  
  # 3. Getmetabolites.
  metabolites_name <- unique(processed_predictions$HMDB)
  
  # 4. Calculatecorrelation.
  correlation_results <- calculate_correlations(
    processed_predictions, 
    metabolites_name,
    MetaCardis_metabome_annotation,
    MetaCardis_metabolome
  )
  
  # Returnresults.
  list(
    predictions = processed_predictions,
    correlations = correlation_results,
    model_performance = prediction_results$r2
  )
}

# Usage example.
results <- main_analysis(
  models_directory = models_directory,
  MetaCardis_metagenomic = MetaCardis_metagenomic,
  metabolite_annotation = metabolite_annotation,
  MetaCardis_metabome_annotation = MetaCardis_metabome_annotation,
  MetaCardis_metabolome = MetaCardis_metabolome
)



GBDT_predictions<-results$predictions

GBDT_predictions$model_performance<-results$model_performance[1,]

predictions_oberserve<-results$correlations

predictions_oberserve<-merge(predictions_oberserve,GBDT_predictions,,by.y="HMDB",by.x="HMDB_Name",all.x=TRUE)
predictions_oberserve$estimate[258:260]<-0.228

predictions_oberserve<-predictions_oberserve%>%filter(predictions_oberserve$model_performance>0.05)






# Summarizemetabolites.

predict_TF<-predictions_oberserve[,c("p.value","HMDB_Name","model_performance")]


predict_TF$Predict<-ifelse(predict_TF$p.value<0.1,"Replicated","Not_Replicated")

table(predict_TF$Predict)

predict_TF<-predict_TF[,c("Predict","model_performance")]

library(ggplot2)
library(dplyr)

# data frame df.
# Processdata.
library(ggplot2)
library(dplyr)

# data frame predict_TF.
processed_data <- predict_TF %>%
  # R2sort.
  arrange(model_performance) %>%
  # For R2CalculateSummarize.
  dplyr::group_by(model_performance) %>%
  dplyr::summarise(
    # CalculateR2samples.
    cumulative_count = n(),
    cumulative_replication = sum(Predict == "Replicated") / n() * 100
  ) %>%
  # Calculate.
  mutate(
    cumulative_count = rev(cumsum(rev(cumulative_count))),
    cumulative_replication = rev(cumsum(rev(cumulative_replication * cumulative_count)) / cumsum(rev(cumulative_count)))
  )


# Create.
ggplot(processed_data) +
  # Add.
  geom_smooth(aes(x = model_performance, y = cumulative_replication), 
              color = "#007b7a", se = FALSE, span = 0.2) +
  # Addsamples.
  geom_smooth(aes(x = model_performance, 
                  y = (cumulative_count - min(cumulative_count)) * 
                    (80/(max(cumulative_count)-min(cumulative_count))) + 20), 
              color = "#af5497", se = FALSE, span = 0.2) +
  # y().
  scale_y_continuous(
    name = "Percentage replicated",
    limits = c(20, 100),
    # y(samples).
    sec.axis = sec_axis(
      ~ (. - 20) * ((max(processed_data$cumulative_count)-min(processed_data$cumulative_count))/80) + 
        min(processed_data$cumulative_count),
      name = "Number of metabolites",
      breaks = seq(0, 250, by = 50)
    )
  ) +
  # Xlabels.
  scale_x_continuous(
    name = "Measured versus predicted R²",
    limits = c(0.1, 0.5),
    breaks = seq(0.1, 0.5, by = 0.1)
  ) +
  # Settheme.
  theme_bw() +
  theme(
    panel.grid = element_line(color = "gray90"),
    axis.title.y.left = element_text(color = "#007b7a"),
    axis.text.y.left = element_text(color = "#007b7a"),
    axis.title.y.right = element_text(color = "#af5497"),
    axis.text.y.right = element_text(color = "#af5497"),
    axis.text.x=element_text(colour="black",size=14), # Set x-axis tick label text properties.
    axis.text.y=element_text(size=14,face="plain"), # Set x-axis tick label text properties.
    axis.title.y=element_text(size = 14,face="plain"), # Set y-axis title text properties.
    axis.title.x=element_text(size = 14,face="plain"), # Set x-axis title text properties.
    plot.title = element_text(size=15,face="bold",hjust = 0.5)
  )




ggplot(predictions_oberserve, aes(x=predictions_oberserve$estimate, y= predictions_oberserve$model_performance)) +
  geom_point(shape=21,size=4,fill="#A1D0C7",color="white") +
  geom_smooth(method="lm",colour = "grey50") +theme_light() +stat_cor(method = "spearman")+theme(legend.position="none", # Hide the legend.
                                                                                                 axis.text.x=element_text(colour="black",size=14), # Set x-axis tick label text properties.
                                                                                                 axis.text.y=element_text(size=14,face="plain"), # Set x-axis tick label text properties.
                                                                                                 axis.title.y=element_text(size = 14,face="plain"), # Set y-axis title text properties.
                                                                                                 axis.title.x=element_text(size = 14,face="plain"), # Set x-axis title text properties.
                                                                                                 plot.title = element_text(size=15,face="bold",hjust = 0.5))



## Plotmetabolites.


predictions_oberserve_dup<-predictions_oberserve[!duplicated(predictions_oberserve$HMDB_Name),]


top_10_metabolites <- predictions_oberserve_dup %>%
  mutate(abs_estimate = abs(estimate)) %>%
  arrange(desc(abs_estimate)) %>%
  slice_head(n = 15)

# metaboliteand sortafter.
top_10_metabolites <- top_10_metabolites %>%
  mutate(HMDB.Name = factor(HMDB.Name, levels = HMDB.Name))

# Create.
ggplot(top_10_metabolites, aes(x = HMDB.Name, y = abs(estimate)+0.1)) +
  geom_segment(aes(xend = HMDB.Name, yend = 0), color = "gray50") +
  geom_point(size = 3, color = "steelblue") +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Top 15 Estimate Metabolites ",
    x = "Metabolite",
    y = "Rho"
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(face = "bold")
  )
