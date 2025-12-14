library(ggplot2)
library(dplyr)
library(vegan)
library(ggpubr)
library(cowplot)

# ============================================================
# Step 1: 准备孔隙率数据 (字典)
# ============================================================
# 假设实体密度 2.93
density_ref <- data.frame(
  Genus = c("Acropora sp.", "Astreopora sp.", "Goniopora sp.", "Porites sp.", 
            "Dipsastraea sp.", "Favites sp.", "Platygyra sp."),
  Density = c(1.57, 1.48, 1.31, 1.31, 1.26, 1.24, 1.22)
) %>%
  mutate(Porosity = (1 - Density / 2.93) * 100)

# ============================================================
# Step 2: 计算每个样品的 Alpha 多样性 (Raw Data)
# ============================================================

# 1. 读取 Group
group_df <- read.csv("group.csv", stringsAsFactors=F)
# 补充 Genus 和 Clade 信息 (复用之前的逻辑)
metadata <- group_df %>%
  mutate(
    Genus = case_when(
      group %in% c("J17") ~ "Astreopora sp.",  
      group %in% c("J19") ~ "Goniopora sp.", 
      group %in% c("J37") ~ "Porites sp.",    
      group %in% c("J43") ~ "Acropora sp.",     
      group %in% c("J54") ~ "Platygyra sp.",   
      group %in% c("J61") ~ "Favites sp.", 
      group %in% c("J66") ~ "Dipsastraea sp.",    
      TRUE ~ "Unknown" # 处理可能的 J77
    ),
    Clade = case_when(
      group %in% c("J54", "J61", "J66") ~ "Robust",
      group %in% c("J17", "J19", "J37", "J43") ~ "Complex",
      TRUE ~ "Unknown"
    )
  ) %>% filter(Clade != "Unknown")

# 2. 读取微生物数据并计算 Chao1
# 细菌
bac_otu <- read.csv("bacteria_data_rarefied_9480_no_zeros.csv", row.names=1)
bac_chao1 <- estimateR(t(bac_otu))[2, ] # vegan 需要行=Sample

# 古菌
arc_otu <- read.csv("archaea_data_rarefied_2913_no_zeros.csv", row.names=1)
arc_chao1 <- estimateR(t(arc_otu))[2, ]

# 3. 整合数据框
df_bac <- data.frame(Sample = names(bac_chao1), Chao1 = bac_chao1, Domain = "Bacteria")
df_arc <- data.frame(Sample = names(arc_chao1), Chao1 = arc_chao1, Domain = "Archaea")
df_all <- rbind(df_bac, df_arc)

# 4. 合并 Metadata 和 孔隙率
plot_data_raw <- df_all %>%
  left_join(metadata, by = c("Sample" = "sample")) %>%
  left_join(density_ref, by = "Genus") %>%
  filter(!is.na(Porosity)) # 去掉没匹配上的

# ============================================================
# Step 3: 绘图函数 (散点+抖动)
# ============================================================

plot_raw_correlation <- function(data, domain_name, color_palette) {
  
  sub_data <- subset(data, Domain == domain_name)
  
  p <- ggplot(sub_data, aes(x = Porosity, y = Chao1)) +
    
    # 1. 拟合线 (线性回归)
    geom_smooth(method = "lm", color = "grey40", fill = "grey85", alpha = 0.5, linetype="dashed") +
    
    # 2. 散点 (关键：使用 geom_jitter 防止点重叠)
    # width = 0.2 表示在 X 轴方向轻微抖动，避免点连成一条死板的直线
    geom_jitter(aes(color = Clade, shape = Genus), size = 4, alpha = 0.8, width = 0.2, stroke = 1) +
    
    # 3. 相关性统计
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 5) +
    
    # 4. 样式
    scale_color_manual(values = color_palette) +
    scale_shape_manual(values = c(13, 8, 1, 2, 4, 0, 7)) +
    
    labs(x = "Skeletal Porosity (%)", 
         y = "Chao1 Richness (Individual Samples)",
         title = paste(domain_name, "- Porosity Correlation")) +
    
    theme_bw() +
    theme(
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

# ============================================================
# Step 4: 执行
# ============================================================

my_colors <- c("Complex" = "#E64B35", "Robust" = "#4DBBD5")

p1 <- plot_raw_correlation(plot_data_raw, "Bacteria", my_colors)
p2 <- plot_raw_correlation(plot_data_raw, "Archaea", my_colors)

final_fig <- plot_grid(p1, p2, ncol = 2, align = "hv", labels = c("A", "B"))
print(final_fig)

ggsave("FigureS5_Porosity_RawData.pdf", final_fig, width = 12, height = 4)
ggsave("FigureS5_Porosity_RawData.png", final_fig, width = 12, height = 4, dpi = 600)


################################################################################
################################################################################
library(vegan)
library(dplyr)

# ============================================================
# Step 1: 准备数据
# ============================================================

# 1. 读取微生物抽平数据 (Bacteria)
bac_otu <- read.csv("bacteria_data_rarefied_9480_no_zeros.csv", row.names=1, check.names=F)
# 转置为 vegan 格式 (行=Sample)
bac_otu_t <- t(bac_otu)

# 2. 准备 Metadata (含孔隙率)
group_df <- read.csv("group.csv", stringsAsFactors=F)

# 你的孔隙率参考表 (从之前代码里拿来的)
density_ref <- data.frame(
  Genus = c("Acropora sp.", "Astreopora sp.", "Goniopora sp.", "Porites sp.", 
            "Dipsastraea sp.", "Favites sp.", "Platygyra sp."),
  Density = c(1.57, 1.48, 1.31, 1.31, 1.26, 1.24, 1.22)
) %>%
  mutate(Porosity = (1 - Density / 2.93) * 100)

# 合并生成完整的 Metadata
metadata <- group_df %>%
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
      group %in% c("J54", "J61", "J66") ~ "Robust",
      group %in% c("J17", "J19", "J37", "J43") ~ "Complex",
      TRUE ~ "Unknown"
    )
  ) %>%
  left_join(density_ref, by="Genus") %>%
  filter(Clade != "Unknown")

# 3. 对齐数据
common <- intersect(rownames(bac_otu_t), metadata$sample)
bac_otu_t <- bac_otu_t[common, ]
metadata <- metadata[match(common, metadata$sample), ]

# ============================================================
# Step 2: 统计检验
# ============================================================

# --- 方法 A: Mantel Test ---
# 1. 计算微生物 Bray-Curtis 距离矩阵
dist_bac <- vegdist(bac_otu_t, method = "bray")

# 2. 计算孔隙率 Euclidean 距离矩阵
dist_por <- dist(metadata$Porosity, method = "euclidean")

# 3. 检验
mantel_res <- mantel(dist_bac, dist_por, method = "spearman", permutations = 999)
print("--- Mantel Test Result (Bacteria) ---")
print(mantel_res)

# --- 方法 B: Adonis (PERMANOVA) 控制变量 ---
# 模型含义：先看 Clade 能解释多少，剩下的残差里看 Porosity 能解释多少
adonis_res <- adonis2(dist_bac ~ Clade + Porosity, data = metadata, permutations = 999)
print("--- Adonis Result (Bacteria) ---")
print(adonis_res)

# 反过来：先看 Porosity
adonis_res_rev <- adonis2(dist_bac ~ Porosity + Clade, data = metadata, permutations = 999)
print("--- Adonis Result (Porosity First) ---")
print(adonis_res_rev)

################################################################################
# 1. 读取微生物抽平数据 (Archaea)
aac_otu <- read.csv("archaea_data_rarefied_2913_no_zeros.csv", row.names=1, check.names=F)
# 转置为 vegan 格式 (行=Sample)
aac_otu_t <- t(aac_otu)

# 2. 准备 Metadata (含孔隙率)
group_df <- read.csv("group.csv", stringsAsFactors=F)

# 你的孔隙率参考表 (从之前代码里拿来的)
density_ref <- data.frame(
  Genus = c("Acropora sp.", "Astreopora sp.", "Goniopora sp.", "Porites sp.", 
            "Dipsastraea sp.", "Favites sp.", "Platygyra sp."),
  Density = c(1.57, 1.48, 1.31, 1.31, 1.26, 1.24, 1.22)
) %>%
  mutate(Porosity = (1 - Density / 2.93) * 100)

# 合并生成完整的 Metadata
metadata <- group_df %>%
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
      group %in% c("J54", "J61", "J66") ~ "Robust",
      group %in% c("J17", "J19", "J37", "J43") ~ "Complex",
      TRUE ~ "Unknown"
    )
  ) %>%
  left_join(density_ref, by="Genus") %>%
  filter(Clade != "Unknown")

# 3. 对齐数据
common <- intersect(rownames(aac_otu_t), metadata$sample)
aac_otu_t <- aac_otu_t[common, ]
metadata <- metadata[match(common, metadata$sample), ]

# ============================================================
# Step 2: 统计检验
# ============================================================

# --- 方法 A: Mantel Test ---
# 1. 计算微生物 Bray-Curtis 距离矩阵
dist_aac <- vegdist(aac_otu_t, method = "bray")

# 2. 计算孔隙率 Euclidean 距离矩阵
dist_por <- dist(metadata$Porosity, method = "euclidean")

# 3. 检验
mantel_res <- mantel(dist_aac, dist_por, method = "spearman", permutations = 999)
print("--- Mantel Test Result (Archaea) ---")
print(mantel_res)

# --- 方法 B: Adonis (PERMANOVA) 控制变量 ---
# 模型含义：先看 Clade 能解释多少，剩下的残差里看 Porosity 能解释多少
adonis_res <- adonis2(dist_aac ~ Clade + Porosity, data = metadata, permutations = 999)
print("--- Adonis Result (Archaea) ---")
print(adonis_res)

# 反过来：先看 Porosity
adonis_res_rev <- adonis2(dist_aac ~ Porosity + Clade, data = metadata, permutations = 999)
print("--- Adonis Result (Porosity First) ---")
print(adonis_res_rev)














################################################################################
################################################################################
# ============================================================
# Step 1: 加载必要的包
# ============================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr) # 用于添加相关性统计
library(stringr) # 用于处理字符串

# ============================================================
# Step 2: 准备数据
# ============================================================

# 1. 手动构建骨骼密度数据 (基于你提供的文献值)
# 确保 Genus 名字与 Table S2 的列名完全一致
density_data <- data.frame(
  Genus = c("Acropora sp.", "Astreopora sp.", "Goniopora sp.", "Porites sp.", 
            "Dipsastraea sp.", "Favites sp.", "Platygyra sp."),
  Density = c(1.57, 1.48, 1.31, 1.31, 1.60, 1.24, 1.22),
  Clade = c("Complex", "Complex", "Complex", "Complex", 
            "Robust", "Robust", "Robust")
)

# 2. 读取并清洗 Table S2
# 假设你的文件名为 TableS2_Alpha_Diversity_Summary.csv
library(readr)
raw_table <- read_csv("TableS2_Alpha_Diversity_Summary.csv", 
                      locale = locale(encoding = "UTF-8"))

# 我们只需要 "Bacteria Chao1" 和 "Archaea Chao1" 这两行
target_rows <- c("Bacteria Chao1", "Archaea Chao1")
clean_table <- raw_table %>%
  filter(Index %in% target_rows) %>%
  select(-Complex, -Robust) # 去掉大类的汇总列，只保留属的列

# 3. 数据重塑 (宽变长) 并提取数值
diversity_long <- clean_table %>%
  pivot_longer(cols = -Index, names_to = "Genus", values_to = "Value_Str") %>%
  mutate(
    # 提取 "Mean" (±号前面的数字)
    Mean_Value = as.numeric(str_extract(Value_Str, "^[0-9.]+")),
    # 标记域 (Domain)
    Domain = ifelse(grepl("Bacteria", Index), "Bacteria", "Archaea")
  ) %>%
  select(Domain, Genus, Mean_Value)

# 4. 合并 密度 和 多样性
plot_data <- merge(diversity_long, density_data, by="Genus")

# 查看数据检查是否正确
print(plot_data)

# ============================================================
# Step 3: 定义绘图函数 (顶刊风格)
# ============================================================

plot_correlation <- function(data, domain_name, color_palette) {
  
  # 筛选对应的数据
  sub_data <- subset(data, Domain == domain_name)
  
  p <- ggplot(sub_data, aes(x = Density, y = Mean_Value)) +
    
    # 1. 拟合线 (线性回归)
    # se=TRUE 显示置信区间，alpha设置透明度
    geom_smooth(method = "lm", color = "grey40", fill = "grey85", alpha = 0.5, linetype="dashed") +
    
    # 2. 散点 (按 Clade 上色，按 Genus 形状)
    geom_point(aes(color = Clade, shape = Genus), size = 5, alpha = 0.9, stroke = 1) +
    
    # 3. 添加相关性统计 (Pearson)
    # label.x.npc 控制位置，防止挡住点
    stat_cor(method = "pearson", 
             label.x.npc = "left", 
             label.y.npc = "top", 
             size = 5,
             p.accuracy = 0.001,
             r.accuracy = 0.01) +
    
    # 4. 颜色与形状
    scale_color_manual(values = color_palette) +
    # 自动生成足够的形状，或者手动指定
    scale_shape_manual(values = c(13, 8, 1, 2, 4, 0, 7)) +
    
    # 5. 标签与美化
    labs(x = expression(bold(paste("Skeletal Density (g cm"^-3, ")"))), 
         y = "Chao1 Richness (Mean)",
         title = paste(domain_name, "Diversity vs. Skeletal Density")) +
    
    theme_bw() +
    theme(
      axis.title = element_text(size = 13, face = "bold"),
      axis.text = element_text(size = 11, color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank() # 去掉次级网格
    )
  
  return(p)
}

# ============================================================
# Step 4: 执行绘图
# ============================================================

# 设置配色 (Complexa=红, Robusta=蓝)
my_colors <- c("Complex" = "#E64B35", "Robust" = "#4DBBD5")

# 1. 细菌图
p_bac <- plot_correlation(plot_data, "Bacteria", my_colors)

# 2. 古菌图
p_arc <- plot_correlation(plot_data, "Archaea", my_colors)

# 3. 组合图片
library(cowplot)
# 使用 align="hv" 保证坐标轴对齐
final_fig5_corr <- plot_grid(p_bac, p_arc, ncol = 2, align = "hv", labels = c("A", "B"))

print(final_fig5_corr)

# 保存
ggsave("Figure5_Density_Correlation.pdf", final_fig5_corr, width = 12, height = 5)
ggsave("Figure5_Density_Correlation.png", final_fig5_corr, width = 12, height = 5, dpi = 300)