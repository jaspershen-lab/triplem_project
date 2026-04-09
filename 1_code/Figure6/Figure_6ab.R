setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")
library(tidyverse)
library(tidymass)
###load("data)
load("3_data_analysis/gut_microbiome/data_preparation/object_cross_section")

gut_object <- object_cross_section
gut_microbiome_metadata<-data.frame(gut_object@sample_info)

gut_microbiome_metadata <- gut_microbiome_metadata %>% filter(!is.na(IRIS))

rownames(gut_microbiome_metadata)<-gut_microbiome_metadata$sample_id

library(phyloseq)

OTU = otu_table(gut_object@expression_data, taxa_are_rows = TRUE)
samples = sample_data(gut_microbiome_metadata)
TAX = tax_table(as.matrix(gut_object@variable_info))
gut_phyloseq<- phyloseq(OTU, TAX, samples)

load("3_data_analysis/oral_microbiome/data_preparation/object_cross_section")

oral_object <- object_cross_section
oral_microbiome_metadata<-data.frame(oral_object@sample_info)

oral_microbiome_metadata <- oral_microbiome_metadata %>% filter(!is.na(IRIS))

rownames(oral_microbiome_metadata)<-oral_microbiome_metadata$sample_id

library(phyloseq)

OTU = otu_table(oral_object@expression_data, taxa_are_rows = TRUE)
samples = sample_data(oral_microbiome_metadata)
TAX = tax_table(as.matrix(oral_object@variable_info))
oral_phyloseq<- phyloseq(OTU, TAX, samples)




## Calculatealpha.


alpha_meas = c("Shannon")


gut_phyloseq <- subset_samples(gut_phyloseq, !(IRIS =="Unknown"))
gut_phyloseq_filt <- filter_taxa(gut_phyloseq,function(x) sum(x) >100 , TRUE )
my_comparisons <- list(c("IS", "IR"))

p<-plot_richness(gut_phyloseq_filt, "IRIS", "IRIS",nrow = 1,measures=alpha_meas)

library(ggpubr)

p$data$IRIS <- factor(p$data$IRIS, levels = c("IS", "IR"))

gut_alpha<-
  ggplot(data=p$data,aes(x=IRIS,y=value,fill=IRIS))+ # "fill="Setcolors.
  stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+ # no no ,Use.
  geom_boxplot(size=0.5,fill="white",outlier.fill="white",outlier.color="white")+ # SizeSetboxplot,fillSetcolors,outlier.filloutlier.colorSet.
  geom_jitter(aes(fill=IRIS),width =0.2,shape = 21,size=2.5)+ # Setscatter plot,width,.
  scale_fill_manual(values = c("# E69F00","0072B2"))+ Setcolors.
  scale_color_manual(values=c("black","black","black"))+ # Setscatter plotcolors.
  # Ggtitle("Gut microbiome")+ Set.
  theme_bw()+ # Use a white background.
  theme(legend.position="none", # Hide the legend.
        axis.text.x=element_text(colour="black",size=14), # Set x-axis tick label text properties.
        axis.text.y=element_text(size=14,face="plain"), # Set x-axis tick label text properties.
        axis.title.y=element_text(size = 14,face="plain"), # Set y-axis title text properties.
        axis.title.x=element_text(size = 14,face="plain"), # Set x-axis title text properties.
        plot.title = element_text(size=15,face="bold",hjust = 0.5), # Set.
        panel.grid.minor = element_blank())+
  ylab("Alpha Diversity Index")+xlab("")+ # Set x- and y-axis titles.
  stat_compare_means(comparisons=my_comparisons,label = "p.value")+
  facet_wrap(~variable,scales = "free")

gut_alpha

ggsave(
  gut_alpha,
  filename = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_6/figure_s6a.pdf"
  ),
  width = 3,
  height = 4
)

alpha_meas = c("Shannon")


oral_phyloseq <- subset_samples(oral_phyloseq, !(IRIS =="Unknown"))
oral_phyloseq_filt <- filter_taxa(oral_phyloseq,function(x) sum(x) >100 , TRUE )
my_comparisons <- list(c("IS", "IR"))

p<-plot_richness(oral_phyloseq_filt, "IRIS", "IRIS",nrow = 1,measures=alpha_meas)

p$data$IRIS <- factor(p$data$IRIS, levels = c("IS", "IR"))

oral_alpha<-ggplot(data=p$data,aes(x=IRIS,y=value,fill=IRIS))+ # "fill="Setcolors.
  stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+ # no no ,Use.
  geom_boxplot(size=0.5,fill="white",outlier.fill="white",outlier.color="white")+ # SizeSetboxplot,fillSetcolors,outlier.filloutlier.colorSet.
  geom_jitter(aes(fill=IRIS),width =0.2,shape = 21,size=2.5)+ # Setscatter plot,width,.
  scale_fill_manual(values = c("# E69F00","0072B2"))+ Setcolors.
  scale_color_manual(values=c("black","black","black"))+ # Setscatter plotcolors.
  # Ggtitle("oral microbiome")+ Set.
  theme_bw()+ # Use a white background.
  theme(legend.position="none", # Hide the legend.
        axis.text.x=element_text(colour="black",size=14), # Set x-axis tick label text properties.
        axis.text.y=element_text(size=14,face="plain"), # Set x-axis tick label text properties.
        axis.title.y=element_text(size = 14,face="plain"), # Set y-axis title text properties.
        axis.title.x=element_text(size = 14,face="plain"), # Set x-axis title text properties.
        plot.title = element_text(size=15,face="bold",hjust = 0.5), # Set.
        panel.grid.minor = element_blank())+
  ylab("Alpha Diversity Index")+xlab("")+ # Set x- and y-axis titles.
  stat_compare_means(comparisons=my_comparisons,label = "p.value")+facet_wrap(~variable,scales = "free")
oral_alpha
ggsave(
  oral_alpha,
  filename = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_6/figure_s6b.pdf"
  ),
  width = 3,
  height = 4
)


####################


dis_bray<- phyloseq::distance(gut_phyloseq_filt, "bray")

## PCoAfor matrix.
dis_bray.pcoa = ordinate(gut_phyloseq_filt, method="NMDS", distance=dis_bray)

## Plot.
bray.pcoa <- plot_ordination(gut_phyloseq_filt, dis_bray.pcoa, color="IRIS" ) + geom_point(size=3)

## Extractdata.
data<-bray.pcoa$data

## column names.
colnames(data)[1:2]<-c("NMDS1","NMDS2")

## Getcoordinates1,2explained variance.
pc1<-""
pc2<-""

data$IRIS <- factor(data$IRIS, levels = c("IS", "IR"))

gut_beta<-
  ggplot(data, aes(NMDS1, NMDS2)) +
  # Plotsamples,colors,size.
  geom_point(aes(fill=IRIS),size=2.5, shape = 21)+
  # Match the legends for shape, border, and fill.
  scale_fill_manual(values=c("#E69F00", "#0072B2"))+
  scale_color_manual(values=c("#E69F00", "#0072B2"))+
  # Set the title and axis label text.
  labs(title="NMDS - The composition of gut microbiome") +
  theme(text=element_text(size=30))+
  # Add horizontal and vertical dashed lines.
  geom_vline(aes(xintercept = 0),linetype="dotted")+
  geom_hline(aes(yintercept = 0),linetype="dotted")+
  # Adjust the background, axes, and legend styling.
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid=element_blank(),
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour 
                                 = "black"),
        axis.title.x=element_text(colour='black', size=14),
        axis.title.y=element_text(colour='black', size=14),
        axis.text=element_text(colour='black',size=14),
        legend.title=element_blank(),
        legend.text=element_text(size=14),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "black"),
        legend.key.height=unit(1.6,"cm"))+
  # Set the title formatting.
  theme(plot.title = element_text(size=14,colour = "black",hjust = 0.5,face = "bold"))+stat_ellipse(aes(color = IRIS),geom = "polygon",level = 0.5,alpha = 0,size=2)
gut_beta
ggsave(
  gut_beta,
  filename = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_6/figure_s6c.pdf"
  ),
  width = 6,
  height = 5
)


dis_bray<- phyloseq::distance(oral_phyloseq_filt, "bray")

## PCoAfor matrix.
dis_bray.pcoa = ordinate(oral_phyloseq_filt, method="NMDS", distance=dis_bray)

## Plot.
bray.pcoa <- plot_ordination(oral_phyloseq_filt, dis_bray.pcoa, color="IRIS" ) + geom_point(size=3)

## Extractdata.
data<-bray.pcoa$data

## column names.
colnames(data)[1:2]<-c("NMDS1","NMDS2")

## Getcoordinates1,2explained variance.
pc1<-""
pc2<-""

data$IRIS <- factor(data$IRIS, levels = c("IS", "IR"))

gut_beta<-
  ggplot(data, aes(NMDS1, NMDS2)) +
  # Plotsamples,colors,size.
  geom_point(aes(fill=IRIS),size=2.5, shape = 21)+
  # Match the legends for shape, border, and fill.
  scale_fill_manual(values=c("#E69F00", "#0072B2"))+
  scale_color_manual(values=c("#E69F00", "#0072B2"))+
  # Set the title and axis label text.
  labs(title="NMDS - The composition of gut microbiome") +
  theme(text=element_text(size=30))+
  # Add horizontal and vertical dashed lines.
  geom_vline(aes(xintercept = 0),linetype="dotted")+
  geom_hline(aes(yintercept = 0),linetype="dotted")+
  # Adjust the background, axes, and legend styling.
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid=element_blank(),
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour 
                                 = "black"),
        axis.title.x=element_text(colour='black', size=14),
        axis.title.y=element_text(colour='black', size=14),
        axis.text=element_text(colour='black',size=14),
        legend.title=element_blank(),
        legend.text=element_text(size=14),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "black"),
        legend.key.height=unit(1.6,"cm"))+
  # Set the title formatting.
  theme(plot.title = element_text(size=14,colour = "black",hjust = 0.5,face = "bold"))+stat_ellipse(aes(color = IRIS),geom = "polygon",level = 0.5,alpha = 0,size=2)
gut_beta
ggsave(
  gut_beta,
  filename = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_6/figure_s6d.pdf"
  ),
  width = 6,
  height = 5
)


# Generate the genus-level table.
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

####only remain the genus level
library(microbiomedataset)

oral_object <-
  oral_object %>%
  summarize_variables(what = "sum_intensity", group_by = "Genus")

##only remain the genus at least 10% samples
dim(oral_object)

non_zero_per <-
  apply(oral_object, 1, function(x) {
    sum(x != 0) / ncol(oral_object)
  })

idx <-
  which(non_zero_per > 0.1)

oral_object <-
  oral_object[idx, ]

oral_object <-
  oral_object %>%
  transform2relative_intensity()

##
##adjust BMI, sex, and IRIS, ethnicity
library(tidyverse)
library(ggpubr)
library(rstatix)

oral_expression_data <-
  extract_expression_data(oral_object) %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

library(plyr)

oral_sample_info <-
  oral_object@sample_info

#######adjust BMI, sex, and IRIS, ethnicity
oral_expression_data <-
  lm_adjust(expression_data = oral_expression_data,
            sample_info = oral_sample_info,
            threads = 3)

oral_temp_object <- oral_object
oral_temp_object@expression_data <- oral_expression_data

rownames(oral_expression_data)<-oral_temp_object@variable_info$Genus
rownames(gut_expression_data)<-gut_temp_object@variable_info$Genus


rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")
setwd("1_code/4_site_merge/")
library(tidymass)
library(tidyverse)
library(readxl)
###load("data)
load("../../3_data_analysis/gut_microbiome/data_preparation/object_cross_section")

gut_object<-object_cross_section

load("../../3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section")

metabolomics_object<-object_cross_section


load("../../3_data_analysis/oral_microbiome/data_preparation/object_cross_section")

oral_object<-object_cross_section

load("../../3_data_analysis/skin_microbiome/data_preparation/object_cross_section")

skin_object<-object_cross_section

load("../../3_data_analysis/nasal_microbiome/data_preparation/object_cross_section")
nasal_object<-object_cross_section
metabolite_annotation<-read_excel("../../3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx")
my_comparisons <- list(c("IS", "IR"))
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

####only remain the genus level
library(microbiomedataset)

oral_object <-
  oral_object %>%
  summarize_variables(what = "sum_intensity", group_by = "Genus")

##only remain the genus at least 10% samples
dim(oral_object)

non_zero_per <-
  apply(oral_object, 1, function(x) {
    sum(x != 0) / ncol(oral_object)
  })

idx <-
  which(non_zero_per > 0.1)

oral_object <-
  oral_object[idx, ]

oral_object <-
  oral_object %>%
  transform2relative_intensity()

####only remain the genus level
library(microbiomedataset)

skin_object <-
  skin_object %>%
  summarize_variables(what = "sum_intensity", group_by = "Genus")

##only remain the genus at least 10% samples
dim(skin_object)

non_zero_per <-
  apply(skin_object, 1, function(x) {
    sum(x != 0) / ncol(skin_object)
  })

idx <-
  which(non_zero_per > 0.1)

skin_object <-
  skin_object[idx, ]

skin_object <-
  skin_object %>%
  transform2relative_intensity()

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

# Combine microbiome data from four body sites for a PCoA plot.

# Read microbiome data from four body regions.
# Assume the files are in the current working directory.
gut_genus<-gut_object@expression_data
rownames(gut_genus)<-gut_object@variable_info$Genus

oral_genus<-oral_object@expression_data
rownames(oral_genus)<-oral_object@variable_info$Genus

skin_genus<-skin_object@expression_data
rownames(skin_genus)<-skin_object@variable_info$Genus

nasal_genus<-nasal_object@expression_data
rownames(nasal_genus)<-nasal_object@variable_info$Genus

# Load required R packages.
library(vegan)      # CalculateNMDS.
library(ggplot2)    # Plot.
library(readr)      # ReadCSV.
library(dplyr)      # dataProcess.
library(tidyr)      # data.


# Add a source label for each dataset.
gut_samples <- colnames(gut_genus)
oral_samples <- colnames(oral_genus)
skin_samples <- colnames(skin_genus)
nasal_samples <- colnames(nasal_genus)

# Combine all data.
# Transpose the matrix so rows are samples and columns are taxa.
gut_t <- t(gut_genus)
oral_t <- t(oral_genus)
skin_t <- t(skin_genus)
nasal_t <- t(nasal_genus)

# Rename samples to avoid duplicates while keeping the original information.
# Assume identical sample names represent different body sites from the same person.
rownames(gut_t) <- paste0(rownames(gut_t), "_gut")
rownames(oral_t) <- paste0(rownames(oral_t), "_oral")
rownames(skin_t) <- paste0(rownames(skin_t), "_skin")
rownames(nasal_t) <- paste0(rownames(nasal_t), "_nasal")

# Create sample type labels.
gut_labels <- data.frame(Sample = rownames(gut_t), Site = "Gut", 
                         Subject = sub("_gut$", "", rownames(gut_t)))
oral_labels <- data.frame(Sample = rownames(oral_t), Site = "Oral", 
                          Subject = sub("_oral$", "", rownames(oral_t)))
skin_labels <- data.frame(Sample = rownames(skin_t), Site = "Skin", 
                          Subject = sub("_skin$", "", rownames(skin_t)))
nasal_labels <- data.frame(Sample = rownames(nasal_t), Site = "Nasal", 
                           Subject = sub("_nasal$", "", rownames(nasal_t)))

# Merge all taxa.
# First ensure all tables share the same taxa columns.
all_species <- unique(c(colnames(gut_t), colnames(oral_t), colnames(skin_t), colnames(nasal_t)))

# Update the helper that fills missing taxa to avoid index errors.
fill_missing_species <- function(df, all_species) {
  # Create a new data frame containing all possible taxa.
  result <- matrix(0, nrow = nrow(df), ncol = length(all_species))
  rownames(result) <- rownames(df)
  colnames(result) <- all_species
  
  # Fill in the existing data.
  common_species <- intersect(colnames(df), all_species)
  for (sp in common_species) {
    result[, sp] <- df[, sp]
  }
  
  # Convert to a data frame and return it.
  return(as.data.frame(result))
}

# Apply the updated function.
gut_complete <- fill_missing_species(gut_t, all_species)
oral_complete <- fill_missing_species(oral_t, all_species)
skin_complete <- fill_missing_species(skin_t, all_species)
nasal_complete <- fill_missing_species(nasal_t, all_species)

# Merge all sample data.
all_data <- rbind(gut_complete, oral_complete, skin_complete, nasal_complete)

# Merge sample labels.
sample_metadata <- rbind(gut_labels, oral_labels, skin_labels, nasal_labels)
rownames(sample_metadata) <- sample_metadata$Sample

# Ensure sample order matches.
sample_metadata <- sample_metadata[rownames(all_data), ]





# Filter gut and oral samples.

sample_metadata<-subset(sample_metadata,Site%in%c("Gut","Oral"))
sample_metadata<-merge(sample_metadata,oral_object@sample_info,by.x="Subject",by.y="sample_id")

sample_metadata<-subset(sample_metadata,IRIS%in%c("IR","IS"))

sample_metadata$IRIS[5]<-"IR"
sample_metadata$IRIS[6]<-"IR"
sample_metadata$IRIS[9]<-"IR"
sample_metadata$IRIS[10]<-"IR"
## Filter paired oral and gut samples.
com_samples<-intersect(gut_object@sample_info$sample_id, oral_object@sample_info$sample_id)


sample_metadata<-subset(sample_metadata,Subject%in%com_samples)

sample_metadata<-subset(sample_metadata,!(Subject%in%c("69-087","70-1005","70-1004")))

all_data<-all_data[sample_metadata$Sample,]



# Calculate Bray-Curtis distances.
bray_dist <- vegdist(all_data, method = "bray")

# Run NMDS analysis.
set.seed(123)  # Set a random seed so results are reproducible.
nmds_result <- metaMDS(bray_dist, k = 2, trymax = 100, autotransform = FALSE)

# Check whether NMDS converged and report the stress value.
cat("NMDS Stress:", nmds_result$stress, "\n")
if(nmds_result$stress > 0.2) {
  warning("NMDS stress > 0.2, indicating poor ordination quality")
} else if(nmds_result$stress > 0.1) {
  cat("NMDS stress is between 0.1 and 0.2, indicating moderate ordination quality\n")
} else {
  cat("NMDS stress < 0.1, indicating good ordination quality\n")
}

# Extract NMDS coordinates.
nmds_df <- as.data.frame(nmds_result$points)
colnames(nmds_df) <- c("NMDS1", "NMDS2")

# Add sample information to the NMDS data.
nmds_df$Sample <- rownames(nmds_df)
nmds_df <- merge(nmds_df, sample_metadata, by = "Sample")






# Plot the NMDS results.
nmds_plot<-
  ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, fill = IRIS, shape = Site),
         color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_point(size = 4, alpha = 0.8,colour="black") +
  scale_fill_manual(values=c("IR"="#E69F00", "IS"="#0072B2")) +
  #geom_path(aes(group = Subject), linetype = "dashed", color = "grey70")+
  scale_shape_manual(values = c(21,22)) +
  labs(
       subtitle = paste("Stress =", round(nmds_result$stress, 3)),
       x = "NMDS1",
       y = "NMDS2") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 12),
    panel.grid.minor = element_blank()
  )
nmds_plot


ggsave(
  nmds_plot,
  filename = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_6/figure_6a.pdf"
  ),
  width = 6,
  height = 5
)


###

subjects <- unique(sample_metadata$Subject)
paired_distances <- data.frame(Subject=character(), 
                               Distance=numeric(), 
                               IRIS=character(),
                               stringsAsFactors=FALSE)

# 2. for ,CalculateoralgutsamplesbetweenBray-Curtis.
for(subj in subjects) {
  # Getsamples.
  subj_samples <- sample_metadata[sample_metadata$Subject == subj, ]
  
  # Ensure both oral and gut samples are present.
  if(nrow(subj_samples) == 2 && all(c("Gut", "Oral") %in% subj_samples$Site)) {
    # Get sample names.
    gut_sample <- subj_samples$Sample[subj_samples$Site == "Gut"]
    oral_sample <- subj_samples$Sample[subj_samples$Site == "Oral"]
    
    # Extract sample data.
    gut_data <- all_data[gut_sample, , drop=FALSE]
    oral_data <- all_data[oral_sample, , drop=FALSE]
    
    # Calculate Bray-Curtis distances.
    bc_dist <- vegdist(rbind(gut_data, oral_data), method="bray")[1]
    
    # Save the results, including IRIS group information.
    iris_group <- unique(subj_samples$IRIS)
    paired_distances <- rbind(paired_distances, 
                              data.frame(Subject=subj, 
                                         Distance=as.numeric(bc_dist), 
                                         IRIS=iris_group))
  }
}

# 3. PlotboxplotIRIS.
library(ggplot2)

# Createboxplot.
paired_distances$IRIS <- factor(paired_distances$IRIS, levels = c("IS", "IR"))
box_plot <- 
  ggplot(paired_distances, aes(x=IRIS, y=Distance, fill=IRIS)) +
  geom_boxplot(alpha=0.7) +
  geom_jitter(width=0.2, size=2, alpha=0.6) +
  scale_fill_manual(values=c("IR"="#E69F00", "IS"="#0072B2")) +
  labs(
       subtitle="",
       x="Group",
       y="Bray-Curtis Distance") +
  theme_bw() +
  theme(
    plot.title = element_text(size=14, face="bold", hjust=0.5),
    plot.subtitle = element_text(size=12, hjust=0.5),
    axis.title = element_text(size=12),
    axis.text = element_text(size=10),
    legend.position = "none")+stat_compare_means(comparisons=my_comparisons,label = "p.value",method = 't.test')

box_plot

ggsave(
  box_plot,
  filename = file.path(
    r4projects::get_project_wd(),
    "4_manuscript/Figures/Figure_6/figure_6b.pdf"
  ),
  width = 2.5,
  height = 4
)
