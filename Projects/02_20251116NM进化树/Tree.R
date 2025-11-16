## ---------------------------
## 脚本目的: 复现带有多层注释环的组合式进化树
## 核心包: ggtree, ggtreeExtra, ggnewscale
## ---------------------------

library(ape)
library(ggnewscale)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(RColorBrewer)
library(tidyverse)

# --- Step 1: 数据加载与预处理 ---

# 加载数据
BacTree	<- read.tree("data/pMAGs_bact_gtdtk_midroot.tree")
dat <- read_tsv("data/pMAGS_tax.tsv")
dataset <- read_tsv("data/pMAGS_Presence_Datasets.txt")
covmax <- read_tsv("data/pMAGs_cov_sum.txt")

# 清洗与整理分类信息
dat$p_c <- if_else(dat$p == "p__Proteobacteria", dat$c, dat$p)
dat$p_c <- gsub(".__", "", dat$p_c, perl = T)
dat$p_c <- gsub("_.$", "", dat$p_c, perl = T)

# 筛选细菌域的数据，并找出主要分类单元
bacDat <- dat %>% filter(d == "d__Bacteria")
mostAbundantTax	<- bacDat %>%
  group_by(p_c) %>% summarise(total	= n()) %>%
  arrange(desc(total)) %>% slice(1:16)
list <- c(mostAbundantTax$p_c, "Nitrospirota")
bacDat	<- bacDat %>%
  mutate(p_c	= if_else(p_c %in% list,  p_c, "Other"))

# 准备数据集来源的注释数据
bacDatset <- dataset %>% filter(MAGs %in% bacDat$MAGs)
bacDatset <- as.data.frame(bacDatset); rownames(bacDatset) <- bacDatset$MAGs; bacDatset$MAGs <- NULL

# 准备丰度的注释数据
bactcov <- covmax %>% filter(MAGs %in% bacDat$MAGs) %>% select(MAGs, count)
bactcov <- as.data.frame(bactcov); rownames(bactcov) <- bactcov$MAGs; bactcov$MAGs <- NULL
bactcov$count <- log10(bactcov$count)

# --- Step 2: 核心绘图：层层叠加 ---

# 按照分类对树进行分组，以便后续按组上色
a	<- split(bacDat$MAGs, bacDat$p_c)
tree	<- groupOTU(BacTree, a)

# 创建一个颜色调色板
getPaletteBact = colorRampPalette(brewer.pal(9, "Set1"))
bactColor	<- getPaletteBact(length(unique(bacDat$p_c)) + 1)
bactColor[1]	<- "black"

# p1: 绘制核心进化树
p1	<-
  ggtree(tree, layout = 'circular', aes(color = group)) +
  geom_tree() + theme_tree() +
  geom_treescale(width	= 0.1) +
  scale_color_manual(values	= bactColor, guide = "none") +
  theme(legend.position = "right")

# p2: 添加第一层注释环（分类）并为下一层做准备
p2 <- p1 +
  new_scale_colour() + new_scale_fill() +
  geom_fruit(
    data = bacDat, pwidth = 0.01, geom = geom_bar,
    mapping = aes(y = MAGs, fill = p_c, x = 1), stat = "identity"
  ) +
  scale_fill_manual(values = bactColor[-1], name = "Taxa") +
  new_scale_colour() + new_scale_fill() # 【关键】为p3准备新的标度

# p3: 添加第二层注释环（数据集来源）并为下一层做准备
p3 <-
  gheatmap(
    p2, bacDatset, width = 0.2, offset = 0.1,
    colnames = FALSE, color = NULL
  ) +
  scale_fill_discrete(na.translate = F, name = "Datasets") +
  new_scale_colour() + new_scale_fill() # 【关键】为p4准备新的标度

# p4: 添加第三层注释环（丰度）
p4 <- gheatmap(
  p3, bactcov, offset = 0.6, width = 0.05,
  colnames = FALSE, color = NULL
) +
  scale_fill_viridis_c(name = "Normalized\nabundance")

# 显示最终图像
print(p4)


# 保存图片
ggsave("Figures/Fig_1c_bacterialTree.pdf", p4, width=10, height=10)
ggsave("Figures/Fig_1c_bacterialTree.png", p4, width=10, height=10)
