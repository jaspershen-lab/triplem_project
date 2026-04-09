rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")

library(tidyverse)
library(tidymass)
library(readxl)


load(
  "3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section"
)

metabolomics_object <- object_cross_section


### Calculatebody sitemicrobiomealpha.

metabolite_annotation <- read_excel(
  "1_code/revise_code/metabolite_annotation.xlsx"
)

load("3_data_analysis/gut_microbiome/data_preparation/object_cross_section")

### Calculatebody sitemicrobiomealpha.

load("3_data_analysis/gut_microbiome/data_preparation/object_cross_section")

gut_object <- object_cross_section


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

library(vegan)
shannon_div <- diversity(t(gut_object@expression_data), index = "shannon")

# Create the results data frame.
results_gut <- data.frame(Sample = names(shannon_div), Shannon = shannon_div)


load("3_data_analysis/oral_microbiome/data_preparation/object_cross_section")

oral_object <- object_cross_section


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


shannon_div <- diversity(t(oral_object@expression_data), index = "shannon")

# Create the results data frame.
results_oral <- data.frame(Sample = names(shannon_div), Shannon = shannon_div)


load("3_data_analysis/skin_microbiome/data_preparation/object_cross_section")

skin_object <- object_cross_section


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







shannon_div <- diversity(t(skin_object@expression_data), index = "shannon")

# Create the results data frame.
results_skin <- data.frame(Sample = names(shannon_div), Shannon = shannon_div)


load("3_data_analysis/nasal_microbiome/data_preparation/object_cross_section")

nasal_object <- object_cross_section


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







shannon_div <- diversity(t(nasal_object@expression_data), index = "shannon")

# Create the results data frame.
results_nasal <- data.frame(Sample = names(shannon_div), Shannon = shannon_div)

demographic_data <- data.frame(metabolomics_object@sample_info)
# Mergealpha diversity.
Sample_ID <- demographic_data[, 1:2]
colnames(Sample_ID)[1] <- "Sample"

alpha_diversity <- Sample_ID %>%
  full_join(results_gut, by = "Sample") %>%
  full_join(results_oral, by = "Sample") %>%
  full_join(results_skin, by = "Sample") %>%
  full_join(results_nasal, by = "Sample")

rownames(alpha_diversity) <- alpha_diversity$Sample

alpha_diversity <- alpha_diversity[, -1:-2]

colnames(alpha_diversity) <- c("gut", "oral", "skin", "nasal")






## Extractmetabolitedata.

#######adjust BMI, sex, and IRIS, ethnicity
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

metabolome_data <- metabolomics_temp_object@expression_data


metabolome_data <- data.frame(t(metabolome_data))

# Filtershared samples.

metabolome_data <- metabolome_data[rownames(alpha_diversity), ]


# Load required R packages.
library(tidyverse)



# Ensure sample IDs.
common_samples <- intersect(rownames(alpha_diversity), rownames(metabolome_data))
if (length(common_samples) == 0) {
  stop("The two datasets do not share common sample IDs")
}
cat("Shared ", length(common_samples), " samples are available for analysis\n")

# Use shared samples to filter the data.
alpha_diversity_filtered <- alpha_diversity[common_samples, ]
metabolome_data_filtered <- metabolome_data[common_samples, ]

# Calculatesitealphaand metabolitesSpearmancorrelation.
sites <- colnames(alpha_diversity_filtered)
metabolites <- colnames(metabolome_data_filtered)

# resultsdata frame.
results <- data.frame(
  Site = character(),
  Metabolite = character(),
  Rho = numeric(),
  P_value = numeric(),
  Adjusted_P_value = numeric(),
  stringsAsFactors = FALSE
)

# Calculatecorrelation.
for (site in sites) {
  alpha_values <- alpha_diversity_filtered[[site]]
  
  for (metabolite in metabolites) {
    metabolite_values <- metabolome_data_filtered[[metabolite]]
    
    # CalculateSpearmancorrelation.
    cor_test <- cor.test(alpha_values,
                         metabolite_values,
                         method = "pearson",
                         exact = FALSE)
    
    # Addresultsdata frame.
    results <- rbind(
      results,
      data.frame(
        Site = site,
        Metabolite = metabolite,
        Rho = cor_test$estimate,
        P_value = cor_test$p.value,
        Adjusted_P_value = NA,
        # NA,after.
        stringsAsFactors = FALSE
      )
    )
  }
}

# For p-valuesBH.
results$Adjusted_P_value <- p.adjust(results$P_value, method = "BH")

results <- subset(results, P_value < 0.05)

results <- merge(results,
                 metabolite_annotation[, c("variable_id",
                                           "HMDB.Name",
                                           "HMDB.Class",
                                           "HMDB.Source.Microbial")],
                 by.x = "Metabolite",
                 by.y = "variable_id")


results <- subset(results, !(HMDB.Class == "NA"))


# Define.
keep_classes <- c(
  "Benzene and substituted derivatives",
  "Carboxylic acids and derivatives",
  "Fatty Acyls",
  "Glycerophospholipids",
  "Organic sulfuric acids and derivatives",
  "Organooxygen compounds",
  "Piperidines",
  "Steroids and steroid derivatives"
)

# in HMDB.Class"Others".
results$HMDB.Class <- ifelse(results$HMDB.Class %in% keep_classes,
                             results$HMDB.Class,
                             "Others")



# Updated visualization for metabolite correlations by HMDB Class using original Rho values
library(ggplot2)
library(dplyr)
library(scales)
library(ggrepel)
library(patchwork) # For combining plots

# Data preprocessing
plot_data <- results %>%
  # Use original Rho values
  mutate(
    # Transform p-values for better visualization
    neg_log_p = -log10(P_value),
    # Size based on p-value significance
    p_size = case_when(P_value < 0.001 ~ 4, P_value < 0.01 ~ 3, P_value < 0.05 ~ 2, TRUE ~ 1)
  ) %>%
  # Calculate the number of metabolites per class and filter
  group_by(HMDB.Class) %>%
  dplyr::mutate(
    class_size = n(),
    class_median_rho = median(Rho, na.rm = TRUE),
    class_pct_significant = mean(Adjusted_P_value < 0.05, na.rm = TRUE) * 100,
    # Rank within each class by absolute correlation (for labeling purposes)
    class_rank = rank(-abs(Rho))
  ) %>%
  ungroup() %>%
  # Filter classes with at least 10 members (adjust this threshold as needed)
  filter(class_size >= 10) %>%
  # Sort by median Rho within class
  arrange(HMDB.Class, desc(Rho)) %>%
  # Add significance markers and labels
  mutate(# Label top 3 metabolites in each class
    label = ifelse(class_rank <= 3 &
                     Adjusted_P_value < 0.05, HMDB.Name, ""))

# Assign sequential numbers for x-axis
plot_data$metabolite_sort <- 1:nrow(plot_data)

# Calculate class ranges for background shading
class_summary <- plot_data %>%
  group_by(HMDB.Class) %>%
  dplyr::summarise(
    start = min(metabolite_sort),
    end = max(metabolite_sort),
    mid = mean(c(
      min(metabolite_sort), max(metabolite_sort)
    )),
    median_rho = median(Rho, na.rm = TRUE),
    pct_significant = mean(Adjusted_P_value < 0.05, na.rm = TRUE) * 100,
    n = n()
  )

# Custom color palette for sites
site_colors <- c(
  "gut" = "#edd064",
  "oral" = "#a1d5b9",
  "skin" = "#f2ccac",
  "nasal" = "#a17db4"
)

# Main plot
p1 <- ggplot(plot_data, aes(x = metabolite_sort, y = Rho)) +
  # Add alternating backgrounds for classes
  geom_rect(
    data = class_summary,
    aes(
      xmin = start,
      xmax = end,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = "gray95",
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  # Add zero line to distinguish positive from negative correlations
  geom_hline(
    yintercept = 0,
    linetype = "solid",
    color = "gray50",
    size = 0.7
  ) +
  # Add points with site colors
  geom_point(aes(color = Site, size = p_size), alpha = 0.85) +
  # Add FDR significance thresholds (both for positive and negative correlations)
  geom_hline(
    yintercept = median(subset(plot_data, Adjusted_P_value == 0.05 &
                                 Rho > 0)$Rho, na.rm = TRUE),
    linetype = "dashed",
    color = "darkred",
    size = 0.5
  ) +
  geom_hline(
    yintercept = median(subset(plot_data, Adjusted_P_value == 0.05 &
                                 Rho < 0)$Rho, na.rm = TRUE),
    linetype = "dashed",
    color = "darkblue",
    size = 0.5
  ) +
  # Add labels for top metabolites
  geom_text_repel(
    data = subset(plot_data, label != ""),
    aes(label = label),
    size = 3.5,
    box.padding = 0.4,
    point.padding = 0.3,
    force = 8,
    max.overlaps = 30,
    segment.color = "grey50",
    segment.size = 0.2,
    min.segment.length = 0.1
  ) +
  # Customize scales
  scale_color_manual(values = site_colors, name = "Site") +
  scale_size_continuous(
    name = "P-value",
    breaks = c(1, 2, 3, 4),
    labels = c("ns", "P < 0.05", "P < 0.01", "P < 0.001"),
    range = c(2, 5)
  ) +
  scale_y_continuous(
    breaks = seq(-1, 1, by = 0.1),
    minor_breaks = seq(-1, 1, by = 0.05),
    limits = c(min(plot_data$Rho) * 1.05, max(plot_data$Rho) * 1.05),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  # Customize theme
  theme_bw() +
  theme(
    panel.grid.minor = element_line(color = "gray95"),
    panel.grid.major = element_line(color = "gray90"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    # axis.text.y = element_text(size = 14),
    # axis.title = element_text(size = 16),
    legend.position = "right",
    legend.box = "vertical",
    # plot.title = element_text(size = 12, face = "bold"),
    # plot.subtitle = element_text(size = 10),
    plot.caption = element_text(hjust = 0)
  ) +
  # Add class names as x-axis labels
  scale_x_continuous(
    breaks = class_summary$mid,
    labels = class_summary$HMDB.Class,
    expand = c(0.01, 0.01)
  ) +
  # Add informative labels
  labs(
    x = "HMDB Class",
    y = "Correlation (Rho)",
    title = "Metabolite Correlations by HMDB Class",
    subtitle = paste0(
      "Showing ",
      nrow(plot_data),
      " metabolites across ",
      length(unique(plot_data$HMDB.Class)),
      " HMDB classes"
    )
  )

# Display the plot
p1

ggsave(p1,
       filename = "4_manuscript/Figures/Figure_1/figure_1e.pdf",
       width = 8,
       height = 6)

library(ggpubr)

# Plotwhether microbesmetabolitesRho.

plot_microbial_source_data <- plot_data

plot_microbial_source_data$HMDB.Source.Microbial <- gsub("NA",
                                                         "FALSE",
                                                         plot_microbial_source_data$HMDB.Source.Microbial)

HMDB.Source.Microbial_color <- c("", "")

library(gghalves)
library(ggsignif)

p <-
  ggplot(data = plot_microbial_source_data,
         aes(x = HMDB.Source.Microbial, y = abs(Rho), fill = HMDB.Source.Microbial)) +
  # stat_boxplot(geom = "errorbar", width = 0.15, aes(color = "black")) +
  geom_boxplot(
    size = 0.5,
    outlier.fill = "white",
    outlier.color = "white"
  ) +
  ggsignif::geom_signif(
    comparisons = list(c("FALSE", "TRUE")),
    map_signif_level = TRUE,
    textsize = 5,
    vjust = 0.5,
    y_position = c(0.5, 0.6),
    step_increase = 0.05
  ) +
  geom_dotplot(
    binaxis = "y",
    stackdir = "center",
    dotsize = 0.9,
    alpha = 0.7
  ) +
  scale_fill_manual(values = c("#3d95d2", "#f16147")) +
  scale_color_manual(values = c("black", "black", "black")) + # Setscatter plotcolors.
  ggtitle("") + # Set the overall title.
  theme_bw() + theme(
    legend.position = "none",
    # Hide the legend.
    axis.text.x = element_text(colour = "black", size =
                                 14, "Helvetica"),
    # Set x-axis tick label text properties.
    axis.text.y = element_text(size = 14, family = "Helvetica"),
    # Set x-axis tick label text properties.
    axis.title.y = element_text(size = 14, family = "Helvetica"),
    # Set y-axis title text properties.
    axis.title.x = element_text(size = 14, family = "Helvetica"),
    # Set x-axis title text properties.
    plot.title = element_text(
      size = 15,
      face = "bold",
      "Helvetica",
      hjust = 0.5
    )
  ) +
  ylab("Absolute correlation") + xlab("Microbial source") + # Set x- and y-axis titles.
  stat_compare_means()
p
ggsave(p,
       filename = "4_manuscript/submission/Microbiome/revision_Figure/figure_1f.pdf",
       width = 8,
       height = 6)


## Plotbody siterho.

# Ridge plot (mountain plot) of Rho values by microbiome site
library(ggplot2)
library(dplyr)
library(ggridges)

# Custom color palette for sites


# Data preparation
# Data preparation
ridge_data <- results %>%
  # Ensure Site is a factor with desired order
  mutate(Site = factor(Site, levels = c("gut", "oral", "skin", "nasal")),
         # Add significance flag
         is_significant = Adjusted_P_value < 0.05)

# Perform statistical tests comparing gut to each other site
# Store p-values from comparisons
site_comparisons <- data.frame(Site = character(),
                               p_value = numeric(),
                               significance = character())

# Get unique sites excluding gut
other_sites <- c("oral", "skin", "nasal")

# Run Wilcoxon tests comparing gut vs each other site
for (site in other_sites) {
  gut_data <- ridge_data$Rho[ridge_data$Site == "gut"]
  site_data <- ridge_data$Rho[ridge_data$Site == site]
  
  # Perform the test on absolute values
  test_result <- wilcox.test(abs(gut_data), abs(site_data))
  
  # Determine significance symbol
  sig_symbol <- if (test_result$p.value < 0.001)
    "***"
  else if (test_result$p.value < 0.01)
    "**"
  else if (test_result$p.value < 0.05)
    "*"
  else
    "ns"
  
  # Add to results
  site_comparisons <- rbind(
    site_comparisons,
    data.frame(
      Site = site,
      p_value = test_result$p.value,
      significance = sig_symbol
    )
  )
}

# Find the maximum x value for positioning the annotations
max_x <- max(abs(ridge_data$Rho))

# Create ridge plot with significance markers
p1 <- ggplot(ridge_data, aes(x = abs(Rho), y = Site, fill = Site)) +
  # Add density ridges
  geom_density_ridges(
    aes(height = ..density..),
    alpha = 0.8,
    scale = 3,
    rel_min_height = 0.01,
    quantile_lines = TRUE,
    quantiles = 1
  ) +
  # Add significance markers
  geom_text(
    data = site_comparisons,
    aes(x = max_x * 0.95, y = Site, label = significance),
    hjust = 0,
    size = 5,
    fontface = "bold"
  ) +
  # Customize colors
  scale_fill_manual(values = site_colors) +
  labs(
    title = "",
    subtitle = "",
    x = "Spearman Correlation Coefficient (Rho)",
    y = NULL,
    caption = "Significance vs gut: * p < 0.05, ** p < 0.01, *** p < 0.001, ns: not significant"
  ) +
  # Customize theme
  theme_ridges(center_axis_labels = TRUE) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.title.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    plot.caption = element_text(hjust = 0, size = 8)
  )



p1

ggsave(p1,
       filename = "4_manuscript/Figures/Figure_1/figure_1g.pdf",
       width = 8,
       height = 6)

# Summarizemetabolitesbar chart.


# Read the data.
# data"metabolite_annotation.csv"in.
# If data,.
data <- metabolite_annotation



# Process(NA).
data <- subset(data, !(HMDB.Class == "NA"))

# Load.
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(forcats)

# data frame:.
# 1. metabolite_annotation - metabolite data.
# 2. plot_data - sitemetabolite data.

# Step 1: Processmetabolite data.
total_data <- metabolite_annotation
total_data <- subset(total_data, !(HMDB.Class == "NA"))

# Summarize.
class_summary <- total_data %>%
  dplyr::count(HMDB.Class) %>%  # dpcountcount.
  mutate(total = sum(n), percentage = n / total * 100) %>%
  arrange(desc(percentage))

# 7,"Others".
top_n <- 7
if (nrow(class_summary) > top_n) {
  top_classes <- class_summary[1:top_n, ]
  others <- data.frame(
    HMDB.Class = "Others",
    n = sum(class_summary$n[(top_n + 1):nrow(class_summary)]),
    total = class_summary$total[1],
    percentage = sum(class_summary$percentage[(top_n + 1):nrow(class_summary)])
  )
  class_summary <- rbind(top_classes, others)
}

# Add"Site","Total".
class_summary <- class_summary %>%
  mutate(Site = "Total")

# EnsureSiteAdd.
print(names(class_summary))  # Checkcolumn nameswhether Site.

# Step 2: Processsitemetabolite data.
site_data <- plot_data
site_data$Site <- tolower(site_data$Site)
expected_sites <- c("gut", "oral", "skin", "nasal")
site_data <- site_data %>% filter(Site %in% expected_sites)

# Extractdatain ,sitedatain .
important_classes <- class_summary$HMDB.Class
if ("Others" %in% important_classes) {
  important_classes <- important_classes[important_classes != "Others"]
}

# Summarizesite.
site_class_summary <- site_data %>%
  dplyr::group_by(Site, HMDB.Class) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  dplyr::group_by(Site) %>%
  dplyr::mutate(total = sum(count), percentage = count / total * 100)

# Processsitedata,Ensure.
site_class_final <- site_class_summary %>%
  # whether .
  mutate(is_important = HMDB.Class %in% important_classes) %>%
  group_by(Site) %>%
  # For site,"Others".
  mutate(HMDB.Class = if_else(is_important, HMDB.Class, "Others")) %>%
  # Summarize"Others".
  group_by(Site, HMDB.Class) %>%
  dplyr::summarise(percentage = sum(percentage), .groups = "drop")

# Checkresults.
print("site_class_final structure:")
print(str(site_class_final))
print("site_class_final first few rows:")
print(head(site_class_final))

# If results,.
if (ncol(site_class_final) < 3 ||
    !"Site" %in% names(site_class_final)) {
  print("Recompute site_class_final with an alternative method")
  
  # : Checkin results.
  site_data_with_important <- site_class_summary %>%
    mutate(is_important = HMDB.Class %in% important_classes)
  
  print("Intermediate data:")
  print(head(site_data_with_important))
  
  # Remove incomplete cases.
  site_data_reclassified <- site_data_with_important %>%
    mutate(HMDB.Class_new = ifelse(is_important, as.character(HMDB.Class), "Others"))
  
  # Remove incomplete cases.
  site_class_final <- site_data_reclassified %>%
    group_by(Site, HMDB.Class_new) %>%
    summarise(percentage = sum(percentage), .groups = "drop") %>%
    rename(HMDB.Class = HMDB.Class_new)
  
  print("New site_class_final:")
  print(head(site_class_final))
}

# Step 3: Mergesitedata.
# Checkdata frame.
print("class_summary columns:")
print(names(class_summary))
print("site_class_final columns:")
print(names(site_class_final))

# Ensuredata framecolumn names.
class_summary_selected <- class_summary %>%
  dplyr::select(Site, HMDB.Class, percentage)
site_class_selected <- site_class_final %>%
  dplyr::select(Site, HMDB.Class, percentage)

# Mergedata.
combined_data <- rbind(class_summary_selected, site_class_selected)

# Ensure"Others"legendin after.
combined_data$HMDB.Class <- fct_relevel(as.factor(combined_data$HMDB.Class), "Others", after = Inf)

# SetSite,"Total".
combined_data$Site <- factor(combined_data$Site, levels = c("Total", expected_sites))

# colors.
unique_classes <- levels(combined_data$HMDB.Class)
num_colors <- length(unique_classes)

# Use.
if (num_colors <= 8) {
  colors <- brewer.pal(max(3, num_colors), "Set1")
} else if (num_colors <= 12) {
  colors <- brewer.pal(max(3, num_colors), "Paired")
} else {
  # colors.
  colors <- c(brewer.pal(8, "Set1"),
              brewer.pal(8, "Set2"),
              brewer.pal(8, "Set3"))
  # If ,Usecolors.
  if (length(colors) < num_colors) {
    colors <- colorRampPalette(colors)(num_colors)
  } else {
    colors <- colors[1:num_colors]
  }
}

# If 6colors.
if (length(colors) >= 6) {
  colors <- colors[-6]
}

# colors.
print(paste("Length of colors:", length(colors)))
print(paste("Length of unique_classes:", length(unique_classes)))

# Ensurecolors.
if (length(colors) != length(unique_classes)) {
  # If colors,Addcolors.
  if (length(colors) < length(unique_classes)) {
    additional_colors_needed <- length(unique_classes) - length(colors)
    additional_colors <- colorRampPalette(colors)(additional_colors_needed)
    colors <- c(colors, additional_colors)
  }
  # If colors,colors.
  else {
    colors <- colors[1:length(unique_classes)]
  }
}

# Check.
print(paste("Adjusted length of colors:", length(colors)))

# Ensure"Others"Use.
if ("Others" %in% unique_classes) {
  names(colors) <- unique_classes  # Remove incomplete cases.
  colors["Others"] <- "# 999999" .
}
combined_data$Site <- factor(combined_data$Site,
                             levels = c("Total", "gut", "oral", "skin", "nasal"))
# Createbar chart.
p <- ggplot(combined_data, aes(x = Site, y = percentage, fill = HMDB.Class)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colors) +
  theme_bw() +
  labs(title = "",
       x = "Site",
       y = "Percent (%)",
       fill = "HMDB.Class") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right") +
  scale_y_continuous(breaks = seq(0, 100, 10))

# Addlabels().
# Calculate.
label_data <- combined_data %>%
  dplyr::group_by(Site) %>%
  dplyr::arrange(desc(HMDB.Class)) %>%
  dplyr::mutate(
    ymin = lag(cumsum(percentage), default = 0),
    ymax = cumsum(percentage),
    pos = (ymin + ymax) / 2
  ) %>%
  dplyr::ungroup()

# Addlabels(>=5%).
label_threshold <- 2
p <- p +
  geom_text(
    data = subset(label_data, percentage >= label_threshold),
    aes(
      x = Site,
      y = pos,
      label = sprintf("%.1f%%", percentage)
    ),
    size = 3,
    color = "white"
  )

# Settheme.
p <-
p + theme(
  axis.text.x = element_text(
    colour = "black",
    size = 14,
    family = "Helvetica"
  ),
  axis.text.y = element_text(size = 14, family = "Helvetica"),
  axis.title.y = element_text(size = 14, family = "Helvetica"),
  axis.title.x = element_text(size = 14, family = "Helvetica"),
  plot.title = element_text(
    size = 15,
    face = "bold",
    family = "Helvetica",
    hjust = 0.5
  )
) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

p

ggsave(p,
       filename = "4_manuscript/Figures/Figure_1/figure_1h.pdf",
       width = 10,
       height = 6)

# Load.
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)

# Usesite_class_finalAnalyze.
# ,EnsureTotalsitedatain.
# If ,Add.

# Checkdatain site.
print(unique(site_class_final$Site))

# 1. data.
# Createsite.
chi_square_test <- function(site1, site2, data) {
  # Extractsitedata.
  site_data <- data %>%
    filter(Site %in% c(site1, site2))
  
  # dataConvert,.
  # Usecount.
  contingency_table <- site_data %>%
    dplyr::select(Site, HMDB.Class, percentage) %>%
    spread(key = Site,
           value = percentage,
           fill = 0)  # Use0.
  
  # CheckConvertaftercounts.
  print(paste("Contingency table for", site1, "vs", site2))
  print(contingency_table)
  
  # HMDB.Class,counts.
  freq_table <- contingency_table %>%
    select(-HMDB.Class) %>%
    as.matrix()
  
  # Remove incomplete cases.
  chi_test <- tryCatch({
    chisq.test(freq_table)
  }, error = function(e) {
    # If (counts),ReturnNA.
    warning(paste(
      "Error in chi-square test for",
      site1,
      "vs",
      site2,
      ":",
      e$message
    ))
    return(list(
      p.value = NA,
      parameter = NA,
      statistic = NA
    ))
  })
  
  # Returnp-values.
  return(
    list(
      site1 = site1,
      site2 = site2,
      p_value = chi_test$p.value,
      df = chi_test$parameter,
      statistic = chi_test$statistic
    )
  )
}

# 2. Getsite.
sites <- unique(site_class_final$Site)
site_pairs <- expand.grid(site1 = sites,
                          site2 = sites,
                          stringsAsFactors = FALSE)
# and .
site_pairs <- site_pairs %>% filter(site1 != site2)

# 3. for .
results <- list()
for (i in 1:nrow(site_pairs)) {
  pair <- site_pairs[i, ]
  cat("Testing", pair$site1, "vs", pair$site2, "\n")
  result <- chi_square_test(pair$site1, pair$site2, site_class_final)
  results[[i]] <- result
}

# 4. resultsConvertdata frame.
result_df <- do.call(rbind, lapply(results, function(x) {
  data.frame(
    site1 = x$site1,
    site2 = x$site2,
    p_value = x$p_value,
    statistic = x$statistic,
    df = x$df,
    stringsAsFactors = FALSE
  )
}))

# results.
print(result_df)

# 5. Createp-valuesmatrixheatmap.
p_value_matrix <- matrix(NA, length(sites), length(sites))
rownames(p_value_matrix) <- sites
colnames(p_value_matrix) <- sites

# matrix.
for (i in 1:nrow(result_df)) {
  row_idx <- which(sites == result_df$site1[i])
  col_idx <- which(sites == result_df$site2[i])
  p_value_matrix[row_idx, col_idx] <- result_df$p_value[i]
}

# For 1(and ).
diag(p_value_matrix) <- 1

# 6. Createheatmap.
# Convertmatrixggplot.
p_value_long <- melt(p_value_matrix)
names(p_value_long) <- c("Site1", "Site2", "p_value")

# Addsignificance.
p_value_long$significance <- "NS"
p_value_long$significance[p_value_long$p_value < 0.05] <- "*"
p_value_long$significance[p_value_long$p_value < 0.01] <- "**"
p_value_long$significance[p_value_long$p_value < 0.001] <- "***"
# NACalculate.
p_value_long$significance[is.na(p_value_long$p_value)] <- "NA"

# Createheatmap,-log10Convertp-values.
p_value_long$neg_log_p <- -log10(p_value_long$p_value)
# NAInf0.
p_value_long$neg_log_p[is.na(p_value_long$neg_log_p) |
                         is.infinite(p_value_long$neg_log_p)] <- 0

# Createheatmap.
p <- ggplot(p_value_long, aes(x = Site2, y = Site1, fill = neg_log_p)) +
  geom_tile(color = "white") +
  geom_text(aes(label = significance),
            color = "black",
            size = 5) +
  scale_fill_gradient2(
    low = "white",
    high = "red",
    mid = "pink",
    midpoint = 1,
    limit = c(0, max(p_value_long$neg_log_p, na.rm = TRUE)),
    name = "-log10(p-value)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1,
      size = 12
    ),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  labs(title = "", x = "", y = "") +
  coord_fixed()

# 7. heatmap.


p

ggsave(p,
       filename = "4_manuscript/Figures/Figure_1/figure_s1_chi.test.pdf",
       width = 5,
       height = 5)

# 8. Createsignificance.
significance_summary <- result_df %>%
  mutate(
    significance = case_when(
      is.na(p_value) ~ "NA",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "NS"
    ),
    comparison = paste(site1, "vs", site2)
  ) %>%
  select(comparison, p_value, statistic, df, significance) %>%
  arrange(p_value)

# Remove incomplete cases.
print(significance_summary)



######revise

# Remove incomplete cases.
library(vegan)
library(tidyverse)

rare_list <- rarecurve(
  otu_table,
  step   = 200,
  sample = min(rowSums(otu_table)),
  label  = FALSE,
  draw   = FALSE   # ⭐ .
)


rare_df <- purrr::imap_dfr(rare_list, function(x, sam) {
  data.frame(
    Sample  = sam,
    Reads   = attr(x, "Subsample"),
    Species = as.numeric(x)
  )
})

ggplot(rare_df, aes(Reads, Species, group = Sample)) +
  geom_line(alpha = 0.3, color = "#2C7BB6", linewidth = 0.6) +
  geom_vline(
    xintercept = min(rowSums(otu_t)),
    linetype = "dashed",
    color = "grey40"
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "Sequencing depth (reads)",
    y = "Observed species"
  )


# Load required R packages.
library(tidyverse)
library(boot)

# Bootstrap function for correlation
boot_cor <- function(data, indices, x_col, y_col) {
  d <- data[indices, ]
  cor(d[[x_col]], d[[y_col]], method = "pearson", use = "complete.obs")
}

# Ensure sample IDs.
common_samples <- intersect(rownames(alpha_diversity), rownames(metabolome_data))
if (length(common_samples) == 0) {
  stop("The two datasets do not share common sample IDs")
}
cat("Shared ", length(common_samples), " samples are available for analysis\n")

# Use shared samples to filter the data.
alpha_diversity_filtered <- alpha_diversity[common_samples, ]
metabolome_data_filtered <- metabolome_data[common_samples, ]

# Calculatesitealphaand metabolitesPearsoncorrelationBootstrap 95% CI.
sites <- colnames(alpha_diversity_filtered)
metabolites <- colnames(metabolome_data_filtered)

# resultsdata frame.
results <- data.frame(
  Site = character(),
  Metabolite = character(),
  Rho = numeric(),
  P_value = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  Adjusted_P_value = numeric(),
  stringsAsFactors = FALSE
)

# SetBootstrap.
n_bootstrap <- 1000  # Bootstrap.
set.seed(123)  # Setrandom seedresultsreproducible.

# Calculate correlations and bootstrap confidence intervals.
cat("Starting correlation analysis and bootstrap 95% CI estimation...\n")

for (site in sites) {
  cat("Processing site:", site, "\n")
  
  for (metabolite in metabolites) {
    # data.
    temp_data <- data.frame(
      alpha = alpha_diversity_filtered[[site]],
      metabolite = metabolome_data_filtered[[metabolite]]
    )
    
    # Remove incomplete cases.
    temp_data <- temp_data[complete.cases(temp_data), ]
    
    if (nrow(temp_data) < 3) {
      warning(paste("Insufficient sample size:", site, "-", metabolite))
      next
    }
    
    # CalculatePearsoncorrelation.
    cor_test <- cor.test(temp_data$alpha,
                         temp_data$metabolite,
                         method = "pearson")
    
    # BootstrapCalculate95% CI.
    boot_results <- boot(
      data = temp_data,
      statistic = boot_cor,
      R = n_bootstrap,
      x_col = "alpha",
      y_col = "metabolite"
    )
    
    # Calculate95% CI(Usepercentile).
    boot_ci <- boot.ci(boot_results, type = "perc", conf = 0.95)
    
    # ExtractCI.
    ci_lower <- boot_ci$percent[4]
    ci_upper <- boot_ci$percent[5]
    
    # Addresultsdata frame.
    results <- rbind(
      results,
      data.frame(
        Site = site,
        Metabolite = metabolite,
        Rho = cor_test$estimate,
        P_value = cor_test$p.value,
        CI_lower = ci_lower,
        CI_upper = ci_upper,
        Adjusted_P_value = NA,
        stringsAsFactors = FALSE
      )
    )
  }
}

cat("Calculation completed!\n")

# For p-valuesBH.
results$Adjusted_P_value <- p.adjust(results$P_value, method = "BH")

# Filtersignificantresults.
results_significant <- subset(results, P_value < 0.05)

# Addmetabolitesannotations.
results_significant <- merge(
  results_significant,
  metabolite_annotation[, c("variable_id",
                            "HMDB.Name",
                            "HMDB.Class",
                            "HMDB.Source.Microbial")],
  by.x = "Metabolite",
  by.y = "variable_id",
  all.x = TRUE
)

# HMDB.ClassNA.
results_significant <- subset(results_significant, !(HMDB.Class == "NA"))

# Saveresults.
write.csv(results, 
          "correlation_results_with_bootstrap_CI.csv", 
          row.names = FALSE)

write.csv(results_significant, 
          "correlation_results_significant_with_bootstrap_CI.csv", 
          row.names = FALSE)

# results.
cat("\n=== Result summary ===\n")
cat("Total number of tested correlations:", nrow(results), "\n")
cat("Number of significant correlations (P < 0.05):", nrow(results_significant), "\n")
cat("\nNumber of significant correlations by site:\n")
print(table(results_significant$Site))

# results.
cat("\nTop 10 significant results:\n")
print(head(results_significant[order(results_significant$P_value), ], 10))

# : correlationforest plot.
library(ggplot2)

# sitetop 10significantresults.
top_results <- results_significant %>%
  group_by(Site) %>%
  arrange(P_value) %>%
  slice_head(n = 10) %>%
  ungroup()

# Createforest plot.
p_forest <- ggplot(top_results, 
                   aes(x = Rho, 
                       y = reorder(paste(Site, HMDB.Name, sep = " - "), Rho),
                       color = Site)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), 
                 height = 0.3, 
                 size = 0.8) +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "gut" = "#edd064",
    "oral" = "#a1d5b9",
    "skin" = "#f2ccac",
    "nasal" = "#a17db4"
  )) +
  labs(
    title = "Top Correlations with 95% Bootstrap CI",
    x = "Pearson Correlation Coefficient",
    y = NULL
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )

print(p_forest)

ggsave(p_forest,
       filename = "4_manuscript/Figures/Figure_1/correlation_forest_plot_with_CI.pdf",
       width = 10,
       height = 8)

cat("\nAnalysis complete. Results have been saved.\n")

# ===== UseBootstrap CIresultsPlot =====.

# Define.
keep_classes <- c(
  "Benzene and substituted derivatives",
  "Carboxylic acids and derivatives",
  "Fatty Acyls",
  "Glycerophospholipids",
  "Organic sulfuric acids and derivatives",
  "Organooxygen compounds",
  "Piperidines",
  "Steroids and steroid derivatives"
)

# Plotdata.
plot_data <- results_significant %>%
  # in HMDB.Class"Others".
  dplyr::mutate(
    HMDB.Class = ifelse(HMDB.Class %in% keep_classes, HMDB.Class, "Others"),
    # Transform p-values for better visualization
    neg_log_p = -log10(P_value),
    # Size based on p-value significance
    p_size = case_when(
      P_value < 0.001 ~ 4,
      P_value < 0.01 ~ 3,
      P_value < 0.05 ~ 2,
      TRUE ~ 1
    )
  ) %>%
  # Calculate the number of metabolites per class and filter
  dplyr::group_by(HMDB.Class) %>%
  dplyr::mutate(
    class_size = n(),
    class_median_rho = median(Rho, na.rm = TRUE),
    class_pct_significant = mean(Adjusted_P_value < 0.05, na.rm = TRUE) * 100,
    # Rank within each class by absolute correlation (for labeling purposes)
    class_rank = rank(-abs(Rho))
  ) %>%
  ungroup() %>%
  # Filter classes with at least 10 members
  filter(class_size >= 10) %>%
  # Sort by median Rho within class
  arrange(HMDB.Class, desc(Rho)) %>%
  # Add significance markers and labels
  dplyr::mutate(
    # Label top 3 metabolites in each class
    label = ifelse(class_rank <= 3 & Adjusted_P_value < 0.05, HMDB.Name, "")
  )

# Assign sequential numbers for x-axis
plot_data$metabolite_sort <- 1:nrow(plot_data)

# Calculate class ranges for background shading
class_summary <- plot_data %>%
  dplyr::group_by(HMDB.Class) %>%
  dplyr::summarise(
    start = min(metabolite_sort),
    end = max(metabolite_sort),
    mid = mean(c(min(metabolite_sort), max(metabolite_sort))),
    median_rho = median(Rho, na.rm = TRUE),
    pct_significant = mean(Adjusted_P_value < 0.05, na.rm = TRUE) * 100,
    n = n()
  )

# Custom color palette for sites
site_colors <- c(
  "gut" = "#edd064",
  "oral" = "#a1d5b9",
  "skin" = "#f2ccac",
  "nasal" = "#a17db4"
)

# Main plot with Bootstrap CI
p1 <- ggplot(plot_data, aes(x = metabolite_sort, y = Rho)) +
  # Add alternating backgrounds for classes
  geom_rect(
    data = class_summary,
    aes(
      xmin = start,
      xmax = end,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = "gray95",
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  # Add zero line to distinguish positive from negative correlations
  geom_hline(
    yintercept = 0,
    linetype = "solid",
    color = "gray50",
    size = 0.7
  ) +
  # Add 95% CI error bars (50% smaller)
  geom_errorbar(
    aes(
      ymin = Rho - (Rho - CI_lower) * 0.5, 
      ymax = Rho + (CI_upper - Rho) * 0.5, 
      color = Site
    ),
    width = 0.5,
    alpha = 0.3,
    size = 0.7
  ) +
  # Add points with site colors
  geom_point(aes(color = Site, size = p_size), alpha = 0.85) +
  # Add FDR significance thresholds (both for positive and negative correlations)
  geom_hline(
    yintercept = median(subset(plot_data, Adjusted_P_value == 0.05 & Rho > 0)$Rho, na.rm = TRUE),
    linetype = "dashed",
    color = "darkred",
    size = 0.5
  ) +
  geom_hline(
    yintercept = median(subset(plot_data, Adjusted_P_value == 0.05 & Rho < 0)$Rho, na.rm = TRUE),
    linetype = "dashed",
    color = "darkblue",
    size = 0.5
  ) +
  # Add labels for top metabolites
  geom_text_repel(
    data = subset(plot_data, label != ""),
    aes(label = label),
    size = 3.5,
    box.padding = 0.4,
    point.padding = 0.3,
    force = 8,
    max.overlaps = 30,
    segment.color = "grey50",
    segment.size = 0.2,
    min.segment.length = 0.1
  ) +
  # Customize scales
  scale_color_manual(values = site_colors, name = "Site") +
  scale_size_continuous(
    name = "P-value",
    breaks = c(1, 2, 3, 4),
    labels = c("ns", "P < 0.05", "P < 0.01", "P < 0.001"),
    range = c(2, 5)
  ) +
  scale_y_continuous(
    breaks = seq(-1, 1, by = 0.1),
    minor_breaks = seq(-1, 1, by = 0.05),
    limits = c(min(plot_data$CI_lower, na.rm = TRUE) * 1.05, 
               max(plot_data$CI_upper, na.rm = TRUE) * 1.05),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  # Customize theme
  theme_bw() +
  theme(
    panel.grid.minor = element_line(color = "gray95"),
    panel.grid.major = element_line(color = "gray90"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.box = "vertical",
    plot.caption = element_text(hjust = 0)
  ) +
  # Add class names as x-axis labels
  scale_x_continuous(
    breaks = class_summary$mid,
    labels = class_summary$HMDB.Class,
    expand = c(0.01, 0.01)
  ) +
  # Add informative labels
  labs(
    x = "HMDB Class",
    y = "Correlation (Rho) with 95% Bootstrap CI",
    title = "Metabolite Correlations by HMDB Class",
    subtitle = paste0(
      "Showing ",
      nrow(plot_data),
      " metabolites across ",
      length(unique(plot_data$HMDB.Class)),
      " HMDB classes (with 95% Bootstrap CI)"
    )
  )

# Display the plot
print(p1)

ggsave(p1,
       filename = "4_manuscript/submission/Microbiome/revision_Figure//figure_1e_with_bootstrap_CI.pdf",
       width = 10,
       height = 6)

cat("\nThe correlation plot with bootstrap CI has been saved!\n")




rownames(sample_info)<-sample_info$sample_id
# Load required R packages.
library(tidyverse)
library(boot)
library(broom)

# ===== Part 1: Spearmancorrelation analysis(for ) =====.

cat("=== Starting Spearman correlation analysis ===\n")

# Ensure sample IDs.
common_samples <- intersect(rownames(alpha_diversity), rownames(metabolome_data))
if (length(common_samples) == 0) {
  stop("The two datasets do not share common sample IDs")
}
cat("Shared ", length(common_samples), " samples are available for analysis\n")

# Use shared samples to filter the data.
alpha_diversity_filtered <- alpha_diversity[common_samples, ]
metabolome_data_filtered <- metabolome_data[common_samples, ]

# Getvariables.
sample_info_filtered <- sample_info[common_samples, ]

# Calculatesitealphaand metabolitesSpearmancorrelation.
sites <- colnames(alpha_diversity_filtered)
metabolites <- colnames(metabolome_data_filtered)

# Spearmanresults.
results_spearman <- data.frame(
  Site = character(),
  Metabolite = character(),
  Rho_spearman = numeric(),
  P_value_spearman = numeric(),
  stringsAsFactors = FALSE
)

cat("Calculating Spearman correlations...\n")
for (site in sites) {
  for (metabolite in metabolites) {
    alpha_values <- alpha_diversity_filtered[[site]]
    metabolite_values <- metabolome_data_filtered[[metabolite]]
    
    # CalculateSpearmancorrelation.
    cor_test <- cor.test(alpha_values, metabolite_values, 
                         method = "spearman", exact = FALSE)
    
    results_spearman <- rbind(
      results_spearman,
      data.frame(
        Site = site,
        Metabolite = metabolite,
        Rho_spearman = cor_test$estimate,
        P_value_spearman = cor_test$p.value,
        stringsAsFactors = FALSE
      )
    )
  }
}

# BH.
results_spearman$Adjusted_P_spearman <- p.adjust(results_spearman$P_value_spearman, 
                                                 method = "BH")

cat("Spearman analysis completed!\n\n")

# ===== Part 2: regression modelAnalyze(variables) =====.

cat("=== Starting regression model analysis (with covariate adjustment) ===\n")

# Bootstrap function for regression coefficient
boot_reg <- function(data, indices, formula_str) {
  d <- data[indices, ]
  model <- lm(as.formula(formula_str), data = d)
  # Returnalpha diversity.
  coef(model)["alpha_div"]
}

# Regressionresults.
results_regression <- data.frame(
  Site = character(),
  Metabolite = character(),
  Beta = numeric(),
  SE = numeric(),
  P_value_regression = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  stringsAsFactors = FALSE
)

# SetBootstrap.
n_bootstrap <- 1000
set.seed(123)

cat("Starting regression analysis and bootstrap 95% CI estimation...\n")

for (site in sites) {
  cat("Processing site:", site, "\n")
  
  for (metabolite in metabolites) {
    # regressiondata.
    reg_data <- data.frame(
      metabolite_value = metabolome_data_filtered[[metabolite]],
      alpha_div = alpha_diversity_filtered[[site]],
      age = sample_info_filtered$adjusted_age,
      sex = sample_info_filtered$Gender,
      bmi = sample_info_filtered$BMI,
      ethnicity = sample_info_filtered$Ethnicity
    )
    
    # Remove incomplete cases.
    reg_data <- reg_data[complete.cases(reg_data), ]
    
    if (nrow(reg_data) < 10) {
      warning(paste("Insufficient sample size:", site, "-", metabolite))
      next
    }
    
    # regression model(age,sex,bmi,ethnicity).
    model <- lm(metabolite_value ~ alpha_div + age + sex + bmi + ethnicity, 
                data = reg_data)
    
    # Extractresults.
    model_summary <- summary(model)
    coef_table <- coef(model_summary)
    
    # Getalpha_div.
    if ("alpha_div" %in% rownames(coef_table)) {
      beta <- coef_table["alpha_div", "Estimate"]
      se <- coef_table["alpha_div", "Std. Error"]
      p_value <- coef_table["alpha_div", "Pr(>|t|)"]
      
      # BootstrapCalculate95% CI.
      formula_str <- "metabolite_value ~ alpha_div + age + sex + bmi + ethnicity"
      
      boot_results <- tryCatch({
        boot(data = reg_data, 
             statistic = boot_reg, 
             R = n_bootstrap,
             formula_str = formula_str)
      }, error = function(e) {
        warning(paste("Bootstrap failed for", site, "-", metabolite, ":", e$message))
        return(NULL)
      })
      
      if (!is.null(boot_results)) {
        # Calculate95% CI.
        boot_ci <- boot.ci(boot_results, type = "perc", conf = 0.95)
        ci_lower <- boot_ci$percent[4]
        ci_upper <- boot_ci$percent[5]
      } else {
        ci_lower <- NA
        ci_upper <- NA
      }
      
      # Addresults.
      results_regression <- rbind(
        results_regression,
        data.frame(
          Site = site,
          Metabolite = metabolite,
          Beta = beta,
          SE = se,
          P_value_regression = p_value,
          CI_lower = ci_lower,
          CI_upper = ci_upper,
          stringsAsFactors = FALSE
        )
      )
    }
  }
}

# BH.
results_regression$Adjusted_P_regression <- p.adjust(results_regression$P_value_regression, 
                                                     method = "BH")

cat("Regression analysis completed!\n\n")

# ===== Part 3: Mergefor results =====.

cat("=== Merging and comparing results ===\n")

# Mergeresults.
results_combined <- full_join(
  results_spearman,
  results_regression,
  by = c("Site", "Metabolite")
)

# Addmetabolitesannotations.
results_combined <- merge(
  results_combined,
  metabolite_annotation[, c("variable_id", "HMDB.Name", "HMDB.Class", 
                            "HMDB.Source.Microbial")],
  by.x = "Metabolite",
  by.y = "variable_id",
  all.x = TRUE
)

# Filterin at least significantresults.
results_combined_sig <- results_combined %>%
  filter(P_value_spearman < 0.05 | P_value_regression < 0.05) %>%
  filter(!is.na(HMDB.Class) & HMDB.Class != "NA")

# Add.
results_combined_sig <- results_combined_sig %>%
  mutate(
    # (whether ).
    direction_consistent = sign(Rho_spearman) == sign(Beta),
    # significant.
    both_significant = (P_value_spearman < 0.05) & (P_value_regression < 0.05),
    # Spearmansignificant.
    only_spearman_sig = (P_value_spearman < 0.05) & (P_value_regression >= 0.05),
    # regressionsignificant.
    only_regression_sig = (P_value_spearman >= 0.05) & (P_value_regression < 0.05),
    # Remove incomplete cases.
    consistency_type = case_when(
      both_significant & direction_consistent ~ "Both significant & consistent",
      both_significant & !direction_consistent ~ "Both significant but inconsistent",
      only_spearman_sig ~ "Only Spearman significant",
      only_regression_sig ~ "Only Regression significant",
      TRUE ~ "Neither significant"
    )
  )

# Saveresults.
write.csv(results_combined_sig, 
          "correlation_comparison_spearman_vs_regression.csv", 
          row.names = FALSE)

# Summarize.
cat("\n=== Result statistics ===\n")
cat("Total number of associations compared:", nrow(results_combined_sig), "\n\n")

consistency_table <- table(results_combined_sig$consistency_type)
print(consistency_table)

cat("\nDirectional consistency ratio:", 
    round(mean(results_combined_sig$direction_consistent, na.rm = TRUE) * 100, 2), "%\n")

# ===== Part 4: for  =====.

library(ggplot2)
library(ggrepel)

# 1. Spearman Rho vs Regression Betascatter plot.
p_scatter <- ggplot(results_combined_sig, 
                    aes(x = Rho_spearman, y = Beta, color = consistency_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  scale_color_manual(
    values = c(
      "Both significant & consistent" = "#2E7D32",
      "Both significant but inconsistent" = "#D32F2F",
      "Only Spearman significant" = "#1976D2",
      "Only Regression significant" = "#F57C00"
    )
  ) +
  labs(
    title = "Spearman Correlation vs Linear Regression",
    subtitle = "Comparing unadjusted vs covariate-adjusted associations",
    x = "Spearman Rho (unadjusted)",
    y = "Regression Beta (adjusted for age, sex, BMI, ethnicity)",
    color = "Consistency"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

print(p_scatter)
ggsave(p_scatter,
       filename = "4_manuscript/Figures/Figure_1/spearman_vs_regression_scatter.pdf",
       width = 10,
       height = 8)

# 2. sitefor .
p_facet <- ggplot(results_combined_sig, 
                  aes(x = Rho_spearman, y = Beta, color = both_significant)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Site, scales = "free") +
  scale_color_manual(
    values = c("TRUE" = "#2E7D32", "FALSE" = "#999999"),
    labels = c("TRUE" = "Both significant", "FALSE" = "At least one NS")
  ) +
  labs(
    title = "Method Comparison by Body Site",
    x = "Spearman Rho",
    y = "Regression Beta",
    color = NULL
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

print(p_facet)
ggsave(p_facet,
       filename = "4_manuscript/Figures/Figure_1/spearman_vs_regression_by_site.pdf",
       width = 10,
       height = 8)

# 3. bar chart.
consistency_summary <- results_combined_sig %>%
  group_by(Site, consistency_type) %>%
  summarise(count = n(), .groups = "drop")

p_bar <- ggplot(consistency_summary, 
                aes(x = Site, y = count, fill = consistency_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    values = c(
      "Both significant & consistent" = "#2E7D32",
      "Both significant but inconsistent" = "#D32F2F",
      "Only Spearman significant" = "#1976D2",
      "Only Regression significant" = "#F57C00"
    )
  ) +
  labs(
    title = "Consistency of Spearman vs Regression Results",
    x = "Body Site",
    y = "Number of Associations",
    fill = "Consistency Type"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

print(p_bar)
ggsave(p_bar,
       filename = "4_manuscript/Figures/Figure_1/consistency_by_site_barplot.pdf",
       width = 10,
       height = 6)

# 4. p-valuesfor .
p_pvalue <- ggplot(results_combined_sig,
                   aes(x = -log10(P_value_spearman), 
                       y = -log10(P_value_regression),
                       color = direction_consistent)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(alpha = 0.6, size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
  scale_color_manual(
    values = c("TRUE" = "#2E7D32", "FALSE" = "#D32F2F"),
    labels = c("TRUE" = "Same direction", "FALSE" = "Opposite direction")
  ) +
  labs(
    title = "P-value Comparison",
    x = "-log10(P-value) Spearman",
    y = "-log10(P-value) Regression",
    color = "Direction"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

print(p_pvalue)
ggsave(p_pvalue,
       filename = "4_manuscript/Figures/Figure_1/pvalue_comparison.pdf",
       width = 8,
       height = 8)

# 5. correlationheatmap - Spearmancorrelation coefficient.
cor_methods <- cor(results_combined_sig$Rho_spearman, 
                   results_combined_sig$Beta, 
                   use = "complete.obs", 
                   method = "spearman")

cat("\nSpearman correlation between the two methods:", round(cor_methods, 3), "\n")

cat("\nAnalysis complete. All results and figures have been saved.\n")
