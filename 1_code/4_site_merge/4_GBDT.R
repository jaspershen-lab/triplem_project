rm(list = ls())
setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")
library(ggbeeswarm)
setwd("1_code/4_site_merge/")



gut_GBDT_results<-readRDS("../../3_data_analysis/gut_microbiome/GBDT/cross_section/gut_GBDT_results")
oral_GBDT_results<-readRDS("../../3_data_analysis/oral_microbiome/GBDT/cross_section/oral_GBDT_results")
skin_GBDT_results<-readRDS("../../3_data_analysis/skin_microbiome/GBDT/cross_section/skin_GBDT_results")
nasal_GBDT_results<-readRDS("../../3_data_analysis/nasal_microbiome/GBDT/cross_section/nasal_GBDT_results")


gut_GBDT_results_R2<-gut_GBDT_results$summary[,c(1,2)]
colnames(gut_GBDT_results_R2)<-c("metabolite","gut")
gut_GBDT_results_R2<-subset(gut_GBDT_results_R2,gut>0.1)

oral_GBDT_results_R2<-oral_GBDT_results$summary[,c(1,2)]
colnames(oral_GBDT_results_R2)<-c("metabolite","oral")
oral_GBDT_results_R2<-subset(oral_GBDT_results_R2,oral>0.1)

skin_GBDT_results_R2<-skin_GBDT_results$summary[,c(1,2)]
colnames(skin_GBDT_results_R2)<-c("metabolite","skin")
skin_GBDT_results_R2<-subset(skin_GBDT_results_R2,skin>0.1)

nasal_GBDT_results_R2<-nasal_GBDT_results$summary[,c(1,2)]
colnames(nasal_GBDT_results_R2)<-c("metabolite","nasal")
nasal_GBDT_results_R2<-subset(nasal_GBDT_results_R2,nasal>0.1)


four_site_GBDT_R2<-merge(gut_GBDT_results_R2,oral_GBDT_results_R2,by = "metabolite",all = TRUE)
four_site_GBDT_R2<-merge(four_site_GBDT_R2,skin_GBDT_results_R2,by = "metabolite",all = TRUE)

four_site_GBDT_R2<-merge(four_site_GBDT_R2,nasal_GBDT_results_R2,by = "metabolite",all = TRUE)
colnames(four_site_GBDT_R2)<-c("metabolite","gut","oral","skin","nasal")

colnames(four_site_GBDT_R2)<-c("metabolite","gut","oral","skin","nasal")

four_site_GBDT_R2<- data.frame(lapply(four_site_GBDT_R2, function(x) ifelse(is.na(x), 0, x)))

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)



# metabolitesin .
data <- four_site_GBDT_R2 %>%
  mutate(Dominant_Factor = case_when(
    gut >= oral & gut >= skin & gut>=nasal ~ "gut",
    oral >= gut & oral>= skin & oral>=nasal~ "oral",
    skin >= oral & skin >= gut & skin>=nasal ~ "skin",
    nasal >= oral & nasal >= skin & nasal>=oral ~ "nasal"
  ))

# Calculatemetabolites R-squared ,main after R-squared descending ordersort.
data <- data %>%
  mutate(Total_R2 = gut + oral + skin + nasal) %>%
  arrange(Dominant_Factor, desc(Total_R2))

# Convert ggplot Plot.
data_long <- data %>%
  pivot_longer(cols = c("gut", "oral", "skin","nasal"),
               names_to = "Factor",
               values_to = "R_squared")

# Plotbar chart.
g1 <- ggplot(data_long, aes(x = factor(metabolite, levels = data$metabolite), y = R_squared, fill = Factor)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = body_site_color, name = "Factor") +
  theme_minimal() +
  labs(x = "256 metabolites (adjusted R² > 10%)", y = "Adjusted R²") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "top",
        panel.spacing = unit(0, "lines"),
        panel.grid = element_blank()) 




## PlotsiteR250samples.
gut_GBDT_results<-readRDS("../../3_data_analysis/gut_microbiome/GBDT/cross_section/gut_GBDT_results")
oral_GBDT_results<-readRDS("../../3_data_analysis/oral_microbiome/GBDT/cross_section/oral_GBDT_results")
skin_GBDT_results<-readRDS("../../3_data_analysis/skin_microbiome/GBDT/cross_section/skin_GBDT_results")
nasal_GBDT_results<-readRDS("../../3_data_analysis/nasal_microbiome/GBDT/cross_section/nasal_GBDT_results")
gut_GBDT_results_R2<-gut_GBDT_results$summary[,c(1,2)]
colnames(gut_GBDT_results_R2)<-c("metabolite","gut")


oral_GBDT_results_R2<-oral_GBDT_results$summary[,c(1,2)]
colnames(oral_GBDT_results_R2)<-c("metabolite","oral")


skin_GBDT_results_R2<-skin_GBDT_results$summary[,c(1,2)]
colnames(skin_GBDT_results_R2)<-c("metabolite","skin")


nasal_GBDT_results_R2<-nasal_GBDT_results$summary[,c(1,2)]
colnames(nasal_GBDT_results_R2)<-c("metabolite","nasal")



four_site_GBDT_R2<-cbind(gut_GBDT_results_R2,oral_GBDT_results_R2$oral,skin_GBDT_results_R2$skin,nasal_GBDT_results_R2$nasal)
colnames(four_site_GBDT_R2)<-c("metabolite","gut","oral","skin","nasal")




library(tidyverse)
library(ggplot2)

# Getbody site50metabolites.
metabolite_analysis <- function(four_site_GBDT_R2) {
  # Getbody site50data.
  gut_data <- four_site_GBDT_R2 %>% 
    top_n(50, gut) %>%
    select(metabolite, gut) %>%
    mutate(Site = "gut",
           Value = gut) %>%
    select(metabolite, Site, Value)
  
  oral_data <- four_site_GBDT_R2 %>% 
    top_n(50, oral) %>%
    select(metabolite, oral) %>%
    mutate(Site = "oral",
           Value = oral) %>%
    select(metabolite, Site, Value)
  
  skin_data <- four_site_GBDT_R2 %>% 
    top_n(50, skin) %>%
    select(metabolite, skin) %>%
    mutate(Site = "skin",
           Value = skin) %>%
    select(metabolite, Site, Value)
  
  nasal_data <- four_site_GBDT_R2 %>% 
    top_n(50, nasal) %>%
    select(metabolite, nasal) %>%
    mutate(Site = "nasal",
           Value = nasal) %>%
    select(metabolite, Site, Value)
  
  # Mergedata.
  combined_data <- bind_rows(gut_data, oral_data, skin_data, nasal_data)
  

  
  combined_data<-data.frame(combined_data)
  combined_data$Site<-as.factor(combined_data$Site)
  summary_stats <- combined_data %>%
    group_by(Site) %>%
    dplyr:: summarise(
      mean = mean(Value),
      sd = sd(Value),
      count = length(Value),  # Use length()  n().
      se = sd(Value)/sqrt(length(Value))
    )
  
  
  ggplot() +
    # Addbar chart.
    geom_bar(data = summary_stats, 
             aes(x = Site, y = mean, fill = Site),
             stat = "identity",
             width = 0.6,
             alpha = 1) +
    # Add.
    geom_errorbar(data = summary_stats,
                  aes(x = Site, 
                      ymin = mean - se, 
                      ymax = mean + se),
                  width = 0.2) +
    # Add.
    geom_quasirandom(data = combined_data,
                     aes(x = Site, y = Value),
                     alpha = 0.8,
                     width = 0.2) +
    # Setcolors.
    scale_fill_manual(values = body_site_color) +
    # Sety.
    scale_y_continuous(expand = c(0, 0)) +  
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 14,family = "Helvetica"),
      axis.title = element_text(size = 14,family = "Helvetica"),
      axis.text.x = element_text(angle = 30, hjust = 1,family = "Helvetica") , # Rotate x-axis labels if group names are long.
      axis.ticks.length = unit(0.25, "cm"),  # Increase tick length.
      axis.ticks = element_line(linewidth = 0.8)  # Increase tick line width.
    ) +
    # Setaxis labels.
    xlab("") +
    ylab("R2")
  
  
}

# Use.
 result <- metabolite_analysis(four_site_GBDT_R2)
 print(result)




  
  
  
  
  
### UpSet
  
  library(UpSetR)
  four_site_GBDT_R2

  metabolite_annotation_mainclass<-subset(metabolite_annotation,HMDB.Class%in%c("Benzene and substituted derivatives","Carboxylic acids and derivatives","Fatty Acyls","Glycerophospholipids","Indoles and derivatives","Organooxygen compounds","Steroids and steroid derivatives"))
  
  
  df<-four_site_GBDT_R2
  
  df <- df %>% mutate(across(everything(), ~ifelse(. > 0.05, 1, 0)))
  
  df$metabolite<-four_site_GBDT_R2$metabolite
  
  
  upset(df, sets = c("gut", "oral", "skin", "nasal"), keep.order = TRUE,sets.bar.color =c("#edd064" , "#a1d5b9" ,"#f2ccac","#a17db4"))
  
  
  
 
  
  
  
  
  
  
  
  
