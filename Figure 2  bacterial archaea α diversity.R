# ============================================================
# Step 1: 加载必要的包
# ============================================================
library(vegan)
library(ggplot2)
library(dplyr)
library(ggsci)  # NPG配色
library(ggpubr) # 统计检验
library(cowplot) # 拼图

# ============================================================
# Step 2: 定义通用分析与绘图函数
# ============================================================

analyze_and_plot_alpha <- function(data_file, group_file, domain_name) {
  
  # 1. 读取抽平后的数据
  # 注意：check.names=FALSE 防止列名被修改
  otu <- read.csv(data_file, header=T, row.names=1, check.names=FALSE)
  
  # 2. 转置 (vegan要求: 行=样本)
  # 如果你的输入文件已经是 行=ASV, 列=Sample，则需要转置
  otu_t <- as.data.frame(t(otu))
  
  # 3. 读取并构建元数据 (Metadata)
  group_raw <- read.csv(group_file, header=T, stringsAsFactors = FALSE)
  
  # --- 核心：定义 Clade 和 Genus (直接复用你 PCoA 的逻辑) ---
  metadata <- group_raw %>%
    mutate(
      Genus = case_when(
        group %in% c("J17") ~ "Astreopora sp.",  
        group %in% c("J19") ~ "Goniopora sp.", 
        group %in% c("J37") ~ "Porites sp.",    
        group %in% c("J43") ~ "Acropora sp.",     
        group %in% c("J54") ~ "Platygyra sp.",   
        group %in% c("J61") ~ "Favites sp.", 
        group %in% c("J66") ~ "Dipsastraea sp.",    
        TRUE ~ "Unknown"
      ),
      Clade = case_when(
        group %in% c("J54", "J61", "J66") ~ "Robusta",
        group %in% c("J17", "J19", "J37", "J43") ~ "Complexa",
        TRUE ~ "Unknown"
      )
    )
  
  # 4. 对齐样品
  common <- intersect(rownames(otu_t), metadata$sample)
  otu_t <- otu_t[common, ]
  metadata <- metadata[match(common, metadata$sample), ]
  
  # 5. 计算 Alpha 多样性
  # Shannon
  shannon <- diversity(otu_t, index = "shannon")
  # Chao1
  chao1 <- estimateR(otu_t)[2, ]
  
  # 合并到绘图数据
  plot_data <- metadata
  plot_data$Shannon <- shannon
  plot_data$Chao1 <- chao1
  
  # 6. 统计检验 (计算大类 Robusta vs Complexa 的 P值)
  # 我们把这个 P 值作为副标题展示，说明整体差异
  stat_shannon <- compare_means(Shannon ~ Clade, data = plot_data, method = "wilcox.test")
  stat_chao1 <- compare_means(Chao1 ~ Clade, data = plot_data, method = "wilcox.test")
  
  p_val_shannon <- paste0("Clade Diff: P = ", stat_shannon$p.format)
  p_val_chao1 <- paste0("Clade Diff: P = ", stat_chao1$p.format)
  
  # ============================================================
  # 绘图部分：分面展示 (Facet Plot)
  # ============================================================
  
  # 内部函数：画单个指标
  plot_metric <- function(df, metric, title_p_val, y_lab) {
    ggplot(df, aes(x = Genus, y = .data[[metric]], fill = Clade)) +
      # 箱线图
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      # 散点图 (展示每个样品，因为n=3，散点很有必要)
      geom_jitter(width = 0.2, size = 2, alpha = 0.8, color="grey30") +
      
      # 核心设计：分面 (Facet)
      # 按照 Clade 分左右两边，scales="free_x" 让每边只显示自己的属
      facet_grid(~Clade, scales = "free_x", space = "free_x") +
      
      # 配色 (NPG)
      scale_fill_npg() +
      
      # 标签
      labs(x = NULL, y = y_lab, 
           title = paste0(domain_name, " - ", metric),
           subtitle = title_p_val) +
      
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, face = "italic"), # 属名斜体
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = "grey95"), # 分面标题背景
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "none", # 不需要图例，分面标题已经说明了
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10, color = "red")
      )
  }
  
  # 画两张图
  p1 <- plot_metric(plot_data, "Shannon", p_val_shannon, "Shannon Index")
  p2 <- plot_metric(plot_data, "Chao1", p_val_chao1, "Chao1 Index")
  
  # 组合
  p_combined <- plot_grid(p1, p2, ncol = 2, align = "h")
  return(p_combined)
}

# ============================================================
# Step 3: 执行分析
# ============================================================

# 1. 细菌分析
# 确保文件名和你硬盘里的一致
p_bacteria <- analyze_and_plot_alpha(
  data_file = "bacteria_data_rarefied_9480.csv", 
  group_file = "group.csv", 
  domain_name = "Bacteria"
)
print(p_bacteria)
# 2. 古菌分析
p_archaea <- analyze_and_plot_alpha(
  data_file = "archaea_data_rarefied_2913.csv", 
  group_file = "group.csv", 
  domain_name = "Archaea"
)
print(p_archaea)
# ============================================================
# Step 4: 保存与查看
# ============================================================

# 查看细菌结果
print(p_bacteria)
ggsave("Figure_Alpha_Bacteria_Genus_Level.pdf", p_bacteria, width = 10, height = 5)

# 查看古菌结果
print(p_archaea)
ggsave("Figure_Alpha_Archaea_Genus_Level.pdf", p_archaea, width = 10, height = 5)
