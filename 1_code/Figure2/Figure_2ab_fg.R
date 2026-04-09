



### adonis
rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")
adonis_r2_all <- readRDS(file = "3_data_analysis/4_site_merge/adonis_r2_all.rds")
adonis_r2_all$bodysite<-as.character(adonis_r2_all$bodysite)


adonis_r2_all[9,]<-c("33.25","All","Total")
adonis_r2_all[10,]<-c("37.22","All","Microbial")



library(ggpattern)
adonis_r2_all$adonis_r2 <- as.numeric(adonis_r2_all$adonis_r2)

adonis_r2_all$bodysite<-factor(adonis_r2_all$bodysite,levels=c("All","gut","oral","skin","nasal"))
adonis_plot <-
  ggplot(adonis_r2_all,
         aes(
           x = bodysite,
           y = adonis_r2,
           fill = bodysite,
           pattern = group
         )) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = 0.7),
    width = 0.6,
    pattern_fill = "black",
    pattern_angle = 45,
    pattern_density = 0.1,
    pattern_spacing = 0.025
  ) +
  geom_text(
    aes(label = paste(adonis_r2, "%", sep = "")),
    position = position_dodge(width = 0.7),
    vjust = -0.5,
    size = 5
  ) +
  scale_fill_manual(values = body_site_color) +
  scale_pattern_manual(values = c("stripe", "none")) +
  labs(x = NULL, y = "Estimated variance (%)") +
  theme_bw() +
  # scale_y_continuous(expand = expansion(0, 0.5)) +
  theme(axis.text.x = element_text(family = "Helvetica", size = 15),
        legend.position = "right")

adonis_plot
ggsave(adonis_plot,
       filename = "4_manuscript/Figures/Figure_2/figure_2a.pdf",
       width = 8,
       height = 6)

rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")

library(tidyverse)
library(tidymass)

###load("data)
load("3_data_analysis/gut_microbiome/data_preparation/object_cross_section")

gut_object <- object_cross_section

load(
  "3_data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section"
)

metabolomics_object <- object_cross_section

dir.create("3_data_analysis/gut_microbiome/spearman/cross_section/",
           recursive = TRUE)

setwd("3_data_analysis/gut_microbiome/spearman/cross_section/")

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


##### run pca for metabolomics and gut microbiome


metabolomics_pca_object <-
  metabolomics_temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  massstat::run_pca()

## Translated comment.
metabolomics_pc <- t(metabolomics_temp_object@expression_data)
gut_pc <- t(gut_temp_object@expression_data)


share_index <- intersect(row.names(metabolomics_pc), row.names(gut_pc))

metabolomics_pc <- metabolomics_pc[share_index, ]
gut_pc <- gut_pc[share_index, ]


# Translated comment.
# Translated comment.

# Translated comment.
metabolomics_matrix <- as.matrix(metabolomics_pc)
gut_matrix <- as.matrix(gut_pc)

# Translated comment.
correlation_results <- apply(metabolomics_matrix, 2, function(metabolite) {
  apply(gut_matrix, 2, function(species) {
    cor.test(metabolite, species, method = "spearman")
  })
})

# Translated comment.
cor_values <- sapply(correlation_results, function(x)
  sapply(x, function(y)
    y$estimate))
p_values <- sapply(correlation_results, function(x)
  sapply(x, function(y)
    y$p.value))

# Translated comment.


cor_values[p_values > 0.05] <- 0

cor_values_filtered <- cor_values[rowSums(cor_values) != 0, ]

cor_values_filtered <- cor_values_filtered[, colSums(cor_values_filtered) != 0]



# Translated comment.
num_count <- apply(cor_values_filtered, 2, function(x)
  length(which(x != 0))) %>% sort(decreasing = TRUE) %>% data.frame()
metabolomics_top100 <- rownames(num_count)[1:100]



metabolomics_pca_object <-
  metabolomics_temp_object[metabolomics_top100, ] %>%
  activate_mass_dataset(what = "variable_info") %>%
  massstat::run_pca()


metabolomics_top100_PC1 <- metabolomics_pca_object$x[, 1:5]
metabolomics_top100_PC1 <- metabolomics_top100_PC1[share_index, ]

##cor(gut_pc,metabolomics_top100_PC1)%>%View()


metabolomics_matrix <- as.matrix(metabolomics_top100_PC1)
gut_matrix <- as.matrix(gut_pc)

# Translated comment.
correlation_results <- apply(metabolomics_matrix, 2, function(metabolite) {
  apply(gut_matrix, 2, function(species) {
    cor.test(metabolite, species, method = "spearman")
  })
})

# Translated comment.
cor_values <- sapply(correlation_results, function(x)
  sapply(x, function(y)
    y$estimate))
p_values <- sapply(correlation_results, function(x)
  sapply(x, function(y)
    y$p.value))

cor_values[p_values > 0.05] <- 0


# top20 genus correlation with metabolomics PC1


genus_top20_mete <- data.frame(cbind(cor_values[, 1], p_values[, 1]))

rownames(genus_top20_mete) <- rownames(p_values)

colnames(genus_top20_mete) <- c("rho", "p_val")

genus_top20_mete <- genus_top20_mete %>% rownames_to_column("row_name")
genus_top30_mete <- genus_top20_mete %>%
  mutate(abs_rho = abs(rho)) %>% # translated comment
  arrange(desc(abs_rho)) %>%       # translated comment
  slice_head(n = 20)           # translated comment

variable_info <- gut_temp_object@variable_info

genus_top30_mete <- merge(genus_top30_mete,
                          variable_info,
                          by.x = "row_name",
                          by.y = "variable_id")

genus_top30_mete <- genus_top30_mete %>%
  arrange(rho)

# Translated comment.
if (!exists("phylum_color")) {
  phylum_color <- setNames(
    c(
      "#E41A1C",
      "#377EB8",
      "#4DAF4A",
      "#984EA3",
      "#FF7F00",
      "#FFFF33",
      "#A65628",
      "#F781BF"
    ),
    unique(genus_top30_mete$Phylum)
  )
}

ggplot(genus_top30_mete, aes(
  x = reorder(Genus, rho),
  y = rho,
  group = Genus
)) +
  geom_segment(aes(y = 0, yend = rho, xend = Genus),
               color = "#00A1D5FF",
               size = 1.5) +  # translated comment
  geom_point(aes(color = Phylum), size = 3) +  # translated comment
  scale_color_manual(values = phylum_color) +  # translated comment
  labs(x = "Bacteria", y = "Rho Value") +    # translated comment
  theme_light(base_size = 14) +  # translated comment
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 10,
      color = "grey30",
      face = "italic"
    ),
    axis.text.y = element_text(size = 12, color = "grey30"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )


## Translated comment.

cor_plot <- cbind(metabolomics_top100_PC1[, "PC1"], gut_matrix[, c("ASV5", "ASV87")])
colnames(cor_plot) <- c("PC1", "ASV5", "ASV87")

Oscillibacter_PC1 <- cor_plot %>%
  ggplot(aes(ASV87, PC1)) +
  geom_point(color = "#edd064",
             show.legend = FALSE,
             size = 6) +
  geom_smooth(color = "#edd064",
              method = "lm",
              show.legend = FALSE) +
  labs(x = "Normalized abundance of Oscillibacter", y = "PC1") + theme_bw() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 16),
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + stat_cor(method = "spearman")
Oscillibacter_PC1
ggsave(
  Oscillibacter_PC1,
  filename = file.path(r4projects::get_project_wd(), "4_manuscript/Figures/Figure_2/figure_2f.pdf"),
  width = 7,
  height = 6
)

Phocaeicola_PC1 <- cor_plot %>%
  ggplot(aes(ASV5, PC1)) +
  geom_point(color = "#edd064",
             show.legend = FALSE,
             size = 6) +
  geom_smooth(color = "#edd064",
              method = "lm",
              show.legend = FALSE) +
  labs(x = "Normalized abundance of Phocaeicola", y = "PC1") + theme_bw() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 16),
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + stat_cor(method = "spearman")
Phocaeicola_PC1
ggsave(
  Phocaeicola_PC1,
  filename = file.path(r4projects::get_project_wd(), "4_manuscript/Figures/Figure_2/figure_2g.pdf"),
  width = 7,
  height = 6
)

# Translated comment.

# Translated comment.
# Translated comment.

# Translated comment.
pca_data <- as.data.frame(metabolomics_pca_object$x[, 1:2])

# Translated comment.
pc1_var <- round(metabolomics_pca_object$sdev[1]^2 / sum(metabolomics_pca_object$sdev^2) * 100,
                 1)
pc2_var <- round(metabolomics_pca_object$sdev[2]^2 / sum(metabolomics_pca_object$sdev^2) * 100,
                 1)

# Translated comment.
top5_bacteria <- genus_top30_mete %>%
  arrange(desc(abs(rho))) %>%
  slice_head(n = 5)

# Translated comment.
top5_bacteria_data <- gut_matrix[share_index, top5_bacteria$row_name, drop = FALSE]
colnames(top5_bacteria_data) <- top5_bacteria$Genus

# Translated comment.
cat("代谢物PCA数据维度:", dim(pca_data), "\n")
cat("前5微生物数据维度:", dim(top5_bacteria_data), "\n")
cat("共有样本数:", length(share_index), "\n")

# Translated comment.

pca_data <- pca_data[rownames(top5_bacteria_data), ]
# Translated comment.
library(vegan)
env_fit <- envfit(pca_data, top5_bacteria_data, perm = 999)

# Translated comment.
env_arrows <- scores(env_fit, "vectors") * sqrt(env_fit$vectors$r) * 2 # translated comment
env_arrows_df <- as.data.frame(env_arrows)
env_arrows_df$Genus <- rownames(env_arrows_df)

# Translated comment.
env_arrows_df <- merge(env_arrows_df,
                       top5_bacteria[, c("row_name", "Genus", "Phylum")],
                       by = "Genus",
                       all.x = TRUE)

# Translated comment.
library(ggplot2)
library(ggrepel)



# Translated comment.
PCA_plot <- 
  ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point(color = "#00BFC4",
             size = 6,
             alpha = 0.7) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray50",
    alpha = 0.5
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "gray50",
    alpha = 0.5
  ) +
  # Translated comment.
  geom_segment(
    data = env_arrows_df,
    aes(
      x = 0,
      y = 0,
      xend = PC1 * 8,
      yend = PC2 * 8
    ),
    arrow = arrow(length = unit(0.25, "cm")),
    size = 1
  ) +
  # Translated comment.
  geom_text_repel(
    data = env_arrows_df,
    aes(x = PC1 * 8, y = PC2 * 8, label = Genus),
    fontface = "italic",
    size = 4,
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = "grey50"
  ) +
  # Translated comment.
  labs(
    title = "Metabolome PCA",
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)"),
    color = "门 (Phylum)"
  ) +
  # Translated comment.
  scale_color_manual(values = phylum_color) +
  # Translated comment.
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.position = "right"
  )

PCA_plot


ggsave(
  PCA_plot,
  filename = file.path(r4projects::get_project_wd(), "4_manuscript/Figures/Figure_2/figure_2b.pdf"),
  width = 6,
  height = 6
)
