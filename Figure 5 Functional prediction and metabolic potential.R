# ============================================================
# Step 1: 加载必要的包
# ============================================================
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(ggsci)
library(cowplot)

# ============================================================
# Step 2: 定义基因字典 (Gene Dictionary)
# ============================================================
# 根据你提供的扩展列表构建
ko_dict <- data.frame(
  KO = c(
    # --- 固氮 (Nitrogen Fixation) ---
    "K02588", "K02586", "K02591", "K00531", "K22896", "K22897", "K22898", "K22899",
    # --- 氨氧化 (Nitrification) ---
     "K28504", "K10945", "K10946", "K10535",
    # --- 反硝化 (Denitrification) ---
    "K00370", "K00371", "K00374", # Nar (Nitrate reductase)
    "K02567", "K02568",           # Nap
    "K00368", "K15864",           # Nir (Nitrite reductase) - 关键步骤
    "K04561", "K02305",           # Nor (NO reductase)
    "K00376",                     # NosZ (N2O reductase) - 终极缺氧指示
    # --- 硫代谢 (Sulfur) ---
    "K11180", # dsrA (Sulfite reductase) - 厌氧/硫循环
    "K17224"  # soxB (Sulfur oxidation)
  ),
  Gene_Name = c(
    "nifH", "nifD", "nifK" , "anfG", "vnfD", "vnfK","vnfG","vnfH",# 将固氮基因归为一大类
     "amoA", "amoB", "amoC", "hao",
    "narG","narH","narI",
    "napA","napB",
    "nirK", "nirS",
    "norB", "norC",
    "nosZ",
    "dsrA", "soxB"
  ),
  Pathway = c(
    rep("N_Fixation", 8),
    rep("Nitrification", 4),
    rep("Denitrification", 10),
    rep("S_Metabolism", 2)
  )
)

# 打印检查
print(ko_dict)

# ============================================================
# Step 3: 数据读取与标准化
# ============================================================

# 1. 读取数据
bac_ko <- read.csv("bacteria_kegg_ko.csv", row.names=1, check.names=F)
arc_ko <- read.csv("archaea_kegg_ko.csv", row.names=1, check.names=F)
group <- read.csv("group.csv", stringsAsFactors=F)

# 2. 标准化函数 (转化成 CPM)
# 这一步是为了让细菌和古菌的数据在同一量级上，且不同样品可比
normalize_cpm <- function(table) {
  sweep(table, 2, colSums(table), "/") * 1e6
}

bac_cpm <- normalize_cpm(bac_ko)
arc_cpm <- normalize_cpm(arc_ko)

# 3. 筛选目标 KO
# 只保留字典里有的行
target_kos <- ko_dict$KO
bac_sub <- bac_cpm[rownames(bac_cpm) %in% target_kos, ]
arc_sub <- arc_cpm[rownames(arc_cpm) %in% target_kos, ]

# 4. 合并细菌和古菌的功能潜力 (Total Functional Potential)
# 逻辑：我们将两者相加，代表整个群落(Holobiont)的潜在能力
# 首先要补全行（有些KO可能只在细菌有，古菌没有，反之亦然）
all_kos <- unique(c(rownames(bac_sub), rownames(arc_sub)))
common_samples <- intersect(colnames(bac_sub), colnames(arc_sub))

# 创建合并矩阵
total_cpm <- matrix(0, nrow=length(all_kos), ncol=length(common_samples),
                    dimnames=list(all_kos, common_samples))

# 填入数据
total_cpm[rownames(bac_sub), common_samples] <- total_cpm[rownames(bac_sub), common_samples] + as.matrix(bac_sub[, common_samples])
total_cpm[rownames(arc_sub), common_samples] <- total_cpm[rownames(arc_sub), common_samples] + as.matrix(arc_sub[, common_samples])

# 5. 按基因名聚合 (Aggregate by Gene Name)
# 因为多个 KO 可能对应同一个基因（如 nifH 有很多变体），我们需要加和
current_dict <- ko_dict[match(rownames(total_cpm), ko_dict$KO), ]
total_gene_cpm <- aggregate(total_cpm, by=list(Gene=current_dict$Gene_Name), FUN=sum)
rownames(total_gene_cpm) <- total_gene_cpm$Gene
total_gene_cpm <- total_gene_cpm[, -1]

# 6. 添加分组信息
plot_data <- as.data.frame(t(total_gene_cpm))
plot_data$Sample <- rownames(plot_data)
plot_data <- merge(plot_data, group, by.x="Sample", by.y="sample")

# 添加 Clade
plot_data <- plot_data %>%
  mutate(Clade = case_when(
    group %in% c("J54", "J61", "J66") ~ "Robust",
    group %in% c("J17", "J19", "J37", "J43") ~ "Complex",
    TRUE ~ "Unknown"
  )) %>% filter(Clade != "Unknown")


# ============================================================
# Step 4: Figure 5A - 功能热图 (带显著性标记)
# ============================================================

# 1. 计算均值和 Z-score
heatmap_data <- plot_data %>%
  group_by(Clade) %>%
  summarise(across(where(is.numeric), mean)) %>%
  as.data.frame()
rownames(heatmap_data) <- heatmap_data$Clade
heatmap_mat <- t(heatmap_data[, -1]) # 转置: 行=基因, 列=组

# 计算 Z-score (按行标准化)
heatmap_z <- t(scale(t(heatmap_mat)))

# 2. 计算显著性 (Wilcoxon Test)
# 我们需要知道哪些基因在两组间有显著差异
sig_label <- c()
p_values <- c()

gene_list <- rownames(heatmap_mat)
for(gene in gene_list) {
  res <- wilcox.test(plot_data[[gene]] ~ plot_data$Clade)
  p_values <- c(p_values, res$p.value)
  # 标记：* P<0.05, ** P<0.01
  if(res$p.value < 0.01) {
    sig_label <- c(sig_label, "**")
  } else if(res$p.value < 0.05) {
    sig_label <- c(sig_label, "*")
  } else {
    sig_label <- c(sig_label, "")
  }
}

# 创建注释矩阵用于 pheatmap 显示
# 我们把显著性标记放在 Robusta 这一列展示 (或者两列都放，这里演示放右边)
display_mat <- matrix("", nrow=nrow(heatmap_z), ncol=ncol(heatmap_z))
# 假设 Robusta 是第2列 (按字母顺序 Complexa, Robusta)
display_mat[, 2] <- sig_label 

# 3. 绘制热图
# 为了美观，我们按照代谢通路给基因排个序
gene_order <- c("nifH", "nifD", "nifK" , "anfG", "vnfD", "vnfK","vnfG","vnfH", # 固氮
                "amoA", "amoB", "amoC", "hao", # 氨氧化
                "narG","narH","narI","napA","napB", "nirK", "nirS", "norB", "norC", "nosZ", # 反硝化
                "dsrA", "soxB") # 硫
# 取交集防止某些基因数据里没有
gene_order <- intersect(gene_order, rownames(heatmap_z))
heatmap_z_ordered <- heatmap_z[gene_order, ]
display_mat_ordered <- display_mat[match(gene_order, rownames(heatmap_z)), ]

# 定义行注释 (基因所属通路)
annotation_row <- data.frame(
  Pathway = ko_dict$Pathway[match(rownames(heatmap_z_ordered), ko_dict$Gene_Name)]
)
rownames(annotation_row) <- rownames(heatmap_z_ordered)

# 绘图
p_heatmap <- pheatmap(heatmap_z_ordered,
                      cluster_rows = FALSE, # 不聚类，按我们要的生理顺序排
                      cluster_cols = FALSE,
                      display_numbers = display_mat_ordered, # 显示星号
                      fontsize_number = 15,
                      number_color = "black",
                      annotation_row = annotation_row,
                      color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(50), # 蓝白红配色
                      main = "Metabolic Potential (Z-score)",
                      silent = TRUE # 不直接画，存起来
)$gtable


# ============================================================
# Step 5: Figure 5B - 关键基因箱线图
# ============================================================

# 挑选几个最有故事的基因
# 1. amoA: 验证古菌保守性 (预期无差异或差异小)
# 2. nirS/nirK/nosZ: 验证厌氧环境 (预期 Robusta 高)
# 3. dsrA: 硫循环/厌氧 (预期 Robusta 高)

key_genes <- c( "nifH", "nifD", "nifK","amoB","amoC", "nirK", "nirS","soxB")
# 检查数据里有没有这些基因
key_genes <- intersect(key_genes, colnames(plot_data))

# 转换长格式用于 ggplot
plot_data_long <- melt(plot_data, id.vars = c("Sample", "group", "Clade"), 
                       measure.vars = key_genes, variable.name = "Gene", value.name = "CPM")

p_box <- ggplot(plot_data_long, aes(x = Clade, y = CPM, fill = Clade)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
  facet_wrap(~Gene, scales = "free_y", ncol = 4) + # 分面展示
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y.npc = 0.95) +
  scale_fill_manual(values = c("Complex" = "#E64B35", "Robust" = "#4DBBD5")) +
  labs(y = "Gene Abundance (CPM)", x = NULL, title = "Key Metabolic Genes") +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 12, face = "bold")
  )

print(p_box)
ggsave("Figure5_Functional_box.pdf", p_box, width = 10, height = 4)
# ============================================================
# Step 6: 组合图片
# ============================================================
# Heatmap 是 grid 对象，ggplot 是 ggplot 对象，用 cowplot 拼
final_fig5 <- plot_grid(p_heatmap, p_box, ncol = 1, rel_heights = c(2.5, 2.5), labels = c("A", "B"))

print(final_fig5)
ggsave("Figure5_Functional_Profile.pdf", final_fig5, width = 10, height = 10)


###############################################################################
# ============================================================
# Step 1: 加载包
# ============================================================
library(dplyr)
library(reshape2)
library(ggplot2)
library(scales)

# ============================================================
# Step 2: 数据准备 (基于你已有的 plot_data)
# ============================================================
# 假设你之前的代码已经运行到了 Step 3，拥有了 'plot_data' (包含 Sample, Clade, group, 和各基因CPM)
# 且 group 列是 "J17" 这样的编号

# 1. 定义基因顺序 (按代谢通路排序，这决定了Y轴顺序)
gene_levels <- c(
  "nifH", "nifD", "nifK" , "anfG", "vnfD", "vnfK","vnfG","vnfH", # 固氮
  "amoA", "amoB", "amoC", "hao",                         # 氨氧化
  "narG","narH","narI","napA","napB", "nirK", "nirS", "norB", "norC", "nosZ", # 反硝化
  "dsrA", "soxB"                                  # 硫代谢
)

# 筛选数据中实际存在的基因
genes_in_data <- intersect(gene_levels, colnames(plot_data))

# 2. 按“属 (Genus)”汇总数据
bubble_data <- plot_data %>%
  # 映射属名
  mutate(Genus = case_when(
    group %in% c("J17") ~ "Astreopora",  
    group %in% c("J19") ~ "Goniopora", 
    group %in% c("J37") ~ "Porites",    
    group %in% c("J43") ~ "Acropora",     
    group %in% c("J54") ~ "Platygyra",   
    group %in% c("J61") ~ "Favites", 
    group %in% c("J66") ~ "Dipsastraea",
        TRUE ~ "Unknown"
  )) %>%
  filter(Genus != "Unknown") %>%
  # 按 Clade 和 Genus 分组求均值
  group_by(Clade, Genus) %>%
  summarise(across(all_of(genes_in_data), mean), .groups="drop")

# 3. 转换为长格式
bubble_long <- melt(bubble_data, id.vars = c("Clade", "Genus"), 
                    variable.name = "Gene", value.name = "CPM")

# 4. 计算 Z-score (按基因标准化，用于颜色)
# 这样可以看出同一个基因在哪个属里最高
bubble_long <- bubble_long %>%
  group_by(Gene) %>%
  mutate(Z_score = (CPM - mean(CPM)) / sd(CPM)) %>%
  ungroup()

# 5. 设置因子水平 (控制排序，关键步骤！)
# (1) Clade 顺序
bubble_long$Clade <- factor(bubble_long$Clade, levels = c("Complex", "Robust"))

# (2) Genus 顺序 (Complex 在左，Robust 在右)
genus_order <- c("Acropora", "Astreopora", "Goniopora", "Porites", 
                 "Dipsastraea", "Favites", "Platygyra")
bubble_long$Genus <- factor(bubble_long$Genus, levels = genus_order)

# (3) Gene 顺序 (反转，因为ggplot从下往上画)
bubble_long$Gene <- factor(bubble_long$Gene, levels = rev(genes_in_data))

# ============================================================
# Step 3: 绘制高级气泡图
# ============================================================

p_bubble <- ggplot(bubble_long, aes(x = Genus, y = Gene)) +
  
  # 1. 气泡层
  geom_point(aes(size = CPM, fill = Z_score), shape = 21, color = "grey30", stroke = 0.3) +
  
  # 2. 分面 (将两大类分开展示)
  facet_grid(~ Clade, scales = "free_x", space = "free_x") +
  
  # 3. 颜色映射 (蓝-白-红，经典的表达量配色)
  scale_fill_gradient2(low = "#3B4992", mid = "white", high = "#EE0000", midpoint = 0,
                       name = "Relative Intensity\n(Z-score)") +
  
  # 4. 大小映射 (控制气泡大小范围)
  scale_size_continuous(range = c(2, 8), name = "Abundance\n(Mean CPM)") +
  
  # 5. 添加代谢通路背景线 (可选，增加可读性)
  # 可以在 AI 里加，或者这里用 geom_hline
  
  # 6. 主题美化
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, face = "italic", color="black"),
    axis.text.y = element_text(size = 10, face = "bold", color="black"),
    axis.title = element_blank(),
    panel.grid.major = element_line(color = "grey92", linetype = "dashed"), # 浅色网格
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA), # 分面背景
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold")
  ) +
  
  # 7. Y轴分组标签 (模拟通路分组)
  # 这里可以利用 strip 或者 AI 后期添加 "Nitrogen Fixation" 等大括号
  labs(title = "Functional Gene Profiling across Coral Genera")

print(p_bubble)

# 保存
ggsave("Figure5B_Genus_Bubble.pdf", p_bubble, width = 10, height = 6)
