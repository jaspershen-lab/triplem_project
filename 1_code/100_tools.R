lm_adjust <-
  function(expression_data,
           sample_info,
           threads = 5) {
    library(future)
    library(furrr)
    # plan(strategy = multisession(workers = threads))
    new_expression_data <-
      rownames(expression_data) %>%
      purrr::map(function(name) {
        # cat(name, " ")
        x = as.numeric(expression_data[name,])
        temp_data =
          data.frame(x = x,
                     sample_info)
        
        
        temp_data$Gender[temp_data$Gender == 'Female'] = 0
        temp_data$Gender[temp_data$Gender == 'Male'] = 1
        temp_data$Gender = as.numeric(temp_data$Gender)
        
        temp_data$IRIS[temp_data$IRIS == 'IR'] = 1
        temp_data$IRIS[temp_data$IRIS == 'IS'] = 2
        temp_data$IRIS[temp_data$IRIS == 'Unknown'] = 0
        temp_data$IRIS = as.numeric(temp_data$IRIS)
        
        adjusted_x <-
          lm(x ~ Gender + BMI + IRIS, data = temp_data) %>%
          residuals()
        adjusted_x
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    colnames(new_expression_data) <-
      colnames(expression_data)
    
    rownames(new_expression_data) <-
      rownames(expression_data)
    new_expression_data
  }





body_site_color = c(
  "gut" = "#edd064",
  "skin" = "#f2ccac",
  "oral" = "#a1d5b9",
  "nasal" = "#a17db4"
)


iris_color = c(
  "IR" = "#bf3f7c",
  "IS" = "#46c1be",
  "Unknown" = "#546672"
)


sex_color <-
  c(
    "Female" = wesanderson::wes_palettes$Rushmore1[3],
    "Male" = wesanderson::wes_palettes$Rushmore1[1]
  )

ethnicity_color <-
  c(
    "Caucasian" = wesanderson::wes_palettes$Darjeeling2[1],
    "Asian" = wesanderson::wes_palettes$Darjeeling2[2],
    "Hispanics" = wesanderson::wes_palettes$Darjeeling2[3],
    "Black" = wesanderson::wes_palettes$Darjeeling2[4]
  )



microbiome_genus_filt <- function(gut_microbiome_table, gut_microbiome_tax, prevalence) {
  # 计算原始表中的列数
  sample_num <- length(colnames(gut_microbiome_table))
  
  # 将tax数据的前6列绑定到microbiome表上
  gut_microbiome_table <- cbind(gut_microbiome_table, gut_microbiome_tax[, 1:6])
  
  # 使用aggregate函数按Genus对数据进行聚合求和
  gut_microbiome_table <- aggregate(gut_microbiome_table[, 1:sample_num], 
                                    by = list(gut_microbiome_table$Genus), 
                                    sum)
  
  # 设置行名为Genus名
  row.names(gut_microbiome_table) <- gut_microbiome_table$Group.1
  
  # 移除由aggregate自动创建的'Group.1'列
  gut_microbiome_table <- gut_microbiome_table[, -1]
  
  # 转置微生物组表，以便进行OTU存在性分析
  gut_microbiome_table <- t(gut_microbiome_table)
  
  # 计算每个OTU的存在性
  otu_presence <- colSums(gut_microbiome_table > 0)
  
  # 计算存在性阈值
  threshold <- prevalence * sample_num
  
  # 基于阈值过滤OTU
  gut_microbiome_table <- gut_microbiome_table[, otu_presence >= threshold]
  
  # 将处理后的表转换为数据框
  gut_microbiome_table <- data.frame(gut_microbiome_table)
  
  # 返回处理后的表
  return(gut_microbiome_table)
}


phylum_name =
  c(
    "Actinobacteria",
    "Bacteroidetes",
    "Cyanobacteria/Chloroplast"  ,
    "Firmicutes",
    "Lentisphaerae",
    "Proteobacteria",
    "Synergistetes",
    "Verrucomicrobia",
    "Campilobacterota",
    "Candidatus_Saccharibacteria",
    "Fusobacteria",
    "Plantae",
    "Tenericutes",
    "Spirochaetes"
  )

phylum_color =
  ggsci::pal_simpsons()(n = length(phylum_name))