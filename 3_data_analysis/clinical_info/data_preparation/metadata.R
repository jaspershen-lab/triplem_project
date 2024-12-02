# 加载必要的包
library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(tidyr)
library(dplyr)

# 创建100个物种的树
set.seed(123)
n_species <- 100
tree <- ape::rtree(n_species)

# 生成更真实的物种名称
generate_species_names <- function(n) {
  genera <- c("Bacteroides", "Prevotella", "Faecalibacterium", "Roseburia", 
              "Ruminococcus", "Clostridium", "Eubacterium", "Blautia", 
              "Bifidobacterium", "Streptococcus", "Lactobacillus", "Escherichia",
              "Akkermansia", "Methanobrevibacter", "Coprococcus")
  
  species <- paste0(sample(genera, n, replace = TRUE), "_",
                    "sp_", sprintf("%03d", 1:n))
  return(species)
}

tree$tip.label <- generate_species_names(n_species)

# 修改树的风格
# 创建分类级别信息（模拟数据）
tax_data <- data.frame(
  label = tree$tip.label,
  Class = sample(c("Bacteroidia", "Clostridia", "Actinobacteria", "Bacilli"), n_species, replace = TRUE),
  Order = sample(c("Bacteroidales", "Clostridiales", "Bifidobacteriales", "Lactobacillales"), n_species, replace = TRUE),
  Family = sample(c("Bacteroidaceae", "Lachnospiraceae", "Ruminococcaceae", "Bifidobacteriaceae"), n_species, replace = TRUE),
  Genus = sapply(strsplit(tree$tip.label, "_"), `[`, 1),
  stringsAsFactors = FALSE
)

# 创建四个身体部位的相关性数据
body_sites <- c("Gut", "Oral", "Skin", "Vaginal")
correlation_data <- data.frame(
  label = rep(tree$tip.label, each = 4),
  site = rep(body_sites, times = n_species),
  correlation = runif(n_species * 4, -1, 1)
) %>%
  pivot_wider(names_from = site, values_from = correlation)

# 计算每个身体部位的显著物种数（假设相关系数绝对值>0.5为显著）
significance_counts <- data.frame(
  label = tree$tip.label,
  Gut = sample(0:100, n_species, replace = TRUE),
  Oral = sample(0:100, n_species, replace = TRUE),
  Skin = sample(0:100, n_species, replace = TRUE),
  Vaginal = sample(0:100, n_species, replace = TRUE)
)

# 将显著性计数转换为长格式
significance_long <- significance_counts %>%
  pivot_longer(cols = -label, names_to = "site", values_to = "count")

# 绘图代码基本保持不变，但需要调整一些参数
# 前面的代码保持不变，直到绘图部分

# 绘图部分的修改
# 首先创建基础树
# 前面的代码保持不变，直到绘图部分

# 绘图部分的修改
# 首先创建基础树
p <- ggtree(tree, layout = "fan", open.angle = 180,branch.length="none") %<+% tax_data + 
  geom_tree(size = 1)

# 添加热图层
p1 <- p + 
  geom_fruit(data = correlation_long,
             geom = geom_tile,
             mapping = aes(y = label, x = site, fill = correlation),
             offset = 0.2,
             width = 0.8) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "Correlation"
  )

# 添加显著性计数柱状图
p2 <- p1 +
  new_scale_fill() +
  geom_fruit(data = significance_long,
             geom = geom_bar,
             mapping = aes(y = label, x = count, fill = site),
             stat = "identity",
             orientation = "y",
             offset = 0.05,
             width = 0.8) +
  scale_fill_manual(
    values = c("Gut" = "#E41A1C", 
               "Oral" = "#4DAF4A", 
               "Skin" = "#377EB8", 
               "Vaginal" = "#984EA3"),
    name = "Body site"
  ) +
  scale_x_continuous(
    breaks = seq(0, 25, by = 2),
    limits = c(0, 25)
  )

# 添加标签并设置整体布局
p3 <- p2 + 
  geom_tiplab(size = 2, offset = 0.5) +
  theme(legend.position = "right",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)) 

# 显示图形
# 注意：在保存图片时，使用特定的宽度和高度以确保图形显示正确
# ggsave("phylogenetic_heatmap.pdf", p3, width = 15, height = 15)
print(p3)