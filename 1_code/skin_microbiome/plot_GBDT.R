rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")



results<-readRDS("3_data_analysis/skin_microbiome/GBDT/cross_section/gut_GBDT_results")

GBDT_R2<-results$summary

library(readxl)
metabolite_annotation<-read_excel("3_data_analysis/plasma_metabolomics/data_preparation/metabolite/variable_info_metabolome_HMDB_class.xlsx")

GBDT_R2<-cbind(GBDT_R2,metabolite_annotation[,-1])


GBDT_R2[,]
#GBDT_R2<-subset(GBDT_R2,r2_mean>0.1)


# Translated comment.
library(ggplot2)
library(dplyr)
library(scales)
library(ggrepel)

# Translated comment.
compress_trans <- function(from = 0, to = 0.1, factor = 0.2) {
  trans_new(
    name = "compress",
    transform = function(x) {
      ifelse(x <= to,
             x * factor,
             (x - to) + (to * factor))
    },
    inverse = function(x) {
      ifelse(x <= to * factor,
             x / factor,
             (x - (to * factor)) + to)
    },
    domain = c(0, Inf)
  )
}

# Translated comment.
plot_data <- GBDT_R2 %>%
  # Translated comment.
  group_by(HMDB.Class) %>%
  mutate(
    class_size = n(),
    # Translated comment.
    label = ifelse(HMDB.Class == "Indoles and derivatives" & r2_mean >= 0.05, 
                   HMDB.Name, "")
  ) %>%
  ungroup() %>%
  filter(class_size >= 20) %>%
  # Translated comment.
  arrange(HMDB.Class) %>%
  # Translated comment.
  mutate(significant = ifelse(r2_mean >= 0.05, "significant", "not_significant"))

# Translated comment.
plot_data$metabolite_sort <- 1:nrow(plot_data)

# Translated comment.
class_num <- table(plot_data$HMDB.Class)
class_range <- c(0)
class_name <- numeric(length(class_num))

for (i in 1:length(class_num)) {
  class_range[i + 1] <- class_range[i] + class_num[i]
  class_name[i] <- class_range[i] + class_num[i] / 2
}

# Translated comment.
p <- ggplot(plot_data, aes(x = metabolite_sort, y = r2_mean)) +
  # Translated comment.
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = 'black'),
    panel.background = element_rect(fill = 'transparent'),
    legend.key = element_rect(fill = 'transparent'),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  ) +
  # Translated comment.
  scale_x_continuous(
    breaks = class_name,
    labels = names(class_num),
    expand = c(0, 0)
  ) +
  # Translated comment.
  scale_y_continuous(
    trans = compress_trans(),
    breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
  ) +
  # Translated comment.
  labs(
    x = NULL,
    y = "R² Value",
    size = "R² value"
  )

# Translated comment.
for (i in 1:(length(class_range) - 1)) {
  p <- p + annotate('rect',
                    xmin = class_range[i],
                    xmax = class_range[i + 1],
                    ymin = -Inf,
                    ymax = Inf,
                    fill = ifelse(i %% 2 == 0, 'gray95', 'gray85'),
                    alpha = 0.5
  )
}

# Translated comment.
p <- p +
  # Translated comment.
  geom_point(
    data = subset(plot_data, significant == "not_significant"),
    aes(size = 2),
    color = "grey70",
    alpha = 0.6
  ) +
  # Translated comment.
  geom_point(
    data = subset(plot_data, significant == "significant"),
    aes(size = 2, color = HMDB.Class, 
        shape = HMDB.Source.Microbial),
    alpha = 0.8
  ) +
  # Translated comment.
  geom_text_repel(
    data = subset(plot_data, label != ""),
    aes(label = label),
    size = 2.5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    force = 10,
    max.overlaps = Inf
  ) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17,"NA"=16)) +
  scale_size(range = c(2, 2)) +
  # Translated comment.
  geom_hline(
    yintercept = 0.05,
    color = 'gray',
    linetype = 2,
    size = 1
  )

# Translated comment.
print(p)

# Translated comment.
labeled_metabolites <- plot_data %>%
  filter(label != "") %>%
  select(HMDB.Name, r2_mean, HMDB.Source.Microbial)

print("标记的代谢物信息：")
print(labeled_metabolites)