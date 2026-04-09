# run



# Translated comment.
selection_results_list <- list(
  "Gut" = gut_result,
  "Oral" = oral_result,
  "Skin" = skin_result,
  "Nasal" = nasal_result
)

# Translated comment.
sites <- c("Gut", "Oral", "Skin", "Nasal")

# Translated comment.
p <- plot_quadrant_metabolite_selection(
  selection_results_list = selection_results_list,
  metabolomics_class = metabolomics_class,
  sites = sites
)


Coinertia_RV<-c(gut_cor_results$coinertia$RV,oral_cor_results$coinertia$RV,nasal_cor_results$coinertia$RV,skin_cor_results$coinertia$RV)
Species_num<-c(length(gut_result$significant_metabolites),length(oral_result$significant_metabolites),length(nasal_result$significant_metabolites),length(skin_result$significant_metabolites))
site<-c("Gut","Oral","Nasal","Skin")
summary_data<-cbind(Coinertia_RV,Species_num,site)



# Translated comment.



p1<-ggplot(summary_data,aes(x=site,y=as.numeric(Coinertia_RV)))+geom_bar(aes(fill=site),stat = "identity", width = 0.8)+scale_fill_manual(values = body_site_color)+theme_minimal() +theme(
  axis.text.x = element_text( size = 10, color = "grey30",face = "bold"),
  axis.text.y = element_text(size = 12, color = "grey30"),
  axis.title = element_text(size = 14, face = "bold"),
  plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
  legend.position = "top",
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10)
)+labs(x = "Site", y = "Coinertia_RV")

p2<-ggplot(summary_data,aes(x=site,y=as.numeric(Species_num)))+geom_bar(aes(fill=site),stat = "identity", width = 0.8)+scale_fill_manual(values = body_site_color)+theme_minimal() +theme(
  axis.text.x = element_text( size = 10, color = "grey30",face = "bold"),
  axis.text.y = element_text(size = 12, color = "grey30"),
  axis.title = element_text(size = 14, face = "bold"),
  plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
  legend.position = "top",
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10)
)+labs(x = "Site", y = "Species num")