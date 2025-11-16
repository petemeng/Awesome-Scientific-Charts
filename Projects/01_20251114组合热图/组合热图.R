## ---------------------------
## 脚本目的: 复现组合式热图 (Composite Heatmap)
## 核心包: tidyverse, patchwork
## ---------------------------

library(tidyverse)
library(scales)
library(RColorBrewer)
library(patchwork) 

# --- Step 1: 数据加载与清洗 ---
# (假设本地已有相关数据文件，此处为数据预处理标准流程)

dat <- read_tsv("data/pMAGs_listMetabolism.gz")
taxAll <- read_tsv("data/pMAGs_tax.tsv")
marker_list <- read_tsv("data/Other/MetabolishMM_MarkerGenes_Selected.txt")

marker_list$Gene_name <- factor(marker_list$Gene_name, levels = unique(marker_list$Gene_name))

# 处理分类信息，添加MAGs数量统计
taxAll_t <- taxAll %>% group_by(c) %>% summarise(n = n()) %>% filter(n >= 10)
taxAll <- merge(taxAll, taxAll_t, by = "c")
taxAll$c <- paste(taxAll$c, " (", taxAll$n, ")", sep = "")
taxAll$c <- gsub("c__", "", taxAll$c)

# 合并主数据表
dat_merge <- dat %>%
  distinct() %>%
  mutate(value = 1) %>%
  select(-c(Gene_name, Gene_description, Metabolism, Pathway)) %>%
  complete(MAGs, KO_ID, fill = list(value = 0)) %>%
  left_join(marker_list, by = join_by(KO_ID)) %>%
  left_join(taxAll, by = join_by(MAGs)) %>%
  select(-n) %>%
  filter(!(KO_ID == "K14469"))

# 设定Y轴（物种）的排序
tax <- dat_merge %>% select(p, c) %>% distinct()
tax$c <- factor(tax$c)
tax$c <- fct_reorder(tax$c, tax$p)

# 计算绘图用的百分比数据
dat_Pivot <- dat_merge %>%
  group_by(MAGs, KO_ID, Gene_name, Pathway_small, Metabolism, c) %>%
  summarise(n = sum(value)) %>%
  ungroup() %>%
  mutate(number = if_else(n > 1, 1, n)) %>%
  group_by(Gene_name, Pathway_small, Metabolism, c) %>%
  summarise(mean_tax = mean(number) * 100, .groups = 'drop') %>%
  na.omit()

dat_Pivot$c <- factor(dat_Pivot$c, levels = levels(tax$c))

# --- Step 2: 核心绘图循环 ---

# 1. 定义子图的排列顺序，确保与论文原图一致
pathway_order <- c(
  "3-HP/4-HB pathway", "Calvin cycle", "Formaldehyde assimilation", "Reductive citrate cycle",
  "Ammonia oxidation", "Ass. nitrate reduction", "Dissimilatory nitrate reduction",
  "Denitrification", "Urea metabolism", "Ass./diss. sulfate reduction",
  "Sulfite oxidation", "Thiosulfate oxidation (SOX)", "CO oxidation"
)

# 2. 为每个 'Metabolism' (代谢大类) 创建一个颜色映射字典
metabolism_list <- dat_Pivot %>% select(Metabolism, Pathway_small) %>% distinct()
unique_metabolisms <- unique(metabolism_list$Metabolism)
metabolism_colors <- setNames(brewer.pal(length(unique_metabolisms), "Dark2"), unique_metabolisms)

# 3. 按照预设顺序循环遍历每个子通路，创建图表并存入列表
plot_list <- list() # 初始化一个空列表，用于存储所有子图
plots_created <- 0  # 初始化一个成功创建的图表计数器

for (current_pathway in pathway_order) { 
  
  # 筛选出当前子通路的数据
  dat_lol <- dat_Pivot %>%
    filter(Pathway_small == current_pathway)
  
  # 如果该通路没有数据，则跳过，不执行后续代码
  if(nrow(dat_lol) == 0) {
    next
  }
  
  # 如果代码运行到这里，说明有数据，可以创建图表
  plots_created <- plots_created + 1 # 计数器加 1
  
  # 根据当前通路的代谢大类，从颜色字典中查找对应的颜色
  current_metabolism <- dat_lol$Metabolism[1]
  plot_color <- metabolism_colors[current_metabolism]
  
  # 创建 ggplot 热图
  p1 <- ggplot(dat_lol, aes(Gene_name, c)) +
    geom_tile(aes(fill = mean_tax), colour = "black", linewidth = 0.25) +
    scale_fill_gradient(
      low = "white",
      high = plot_color,
      limits = c(0, 100),
      name = "Proportion of\nMAGs (%)"
    ) +
    coord_fixed(ratio = 1) +
    labs(x = NULL, y = NULL, title = current_pathway) +
    theme_minimal() +
    theme(
      text = element_text(size = 12, colour = "black"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=10),
      axis.text.y = element_text(size=10),
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position = "right",
      plot.margin = unit(c(5.5, 1, 5.5, 1), "pt")
    )
  
  # 使用新的计数器来判断
  # 只有在创建的图超过1个时，才移除Y轴标签
  if (plots_created > 1) {
    p1 <- p1 + theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  }
  
  # 将生成好的图 p1 添加到列表中
  plot_list[[current_pathway]] <- p1
}

# 4. 使用 patchwork 将列表中的所有图拼合成一个大图
if (length(plot_list) > 0) {
  # 注意：此版本未实现自适应宽度，所有面板宽度将基本一致
  combined_plot <- wrap_plots(plot_list, nrow = 1) +
    plot_layout(guides = 'collect') &
    theme(legend.position = 'right')

  # 5. 保存最终的组合图
  ggsave("Figures/Fig3a_Combined_Corrected.pdf", combined_plot, width = 18, height = 7, units = "in")
  ggsave("Figures/Fig3a_Combined_Corrected.png", combined_plot, width = 18, height = 7, units = "in", dpi = 300)

  # 在 RStudio 的绘图窗口中显示最终图像
  print(combined_plot)
} else {
  print("没有找到任何通路的数据，无法生成图表。")
}
