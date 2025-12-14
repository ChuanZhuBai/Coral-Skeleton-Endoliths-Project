# ============================================================
# Step 1: 加载必要的包
# ============================================================
library(vegan)
library(ggplot2)
library(dplyr)
library(ggsci)  # 用于Nature风格配色
library(ggpubr) # 用于美化图形

# ============================================================
# Step 2: 读取并处理数据
# ============================================================
dir()
# 1. 读取物种丰度表 (请修改路径)
# 注意：假设你的第一列是ASV ID，设为行名
# 正确读取制表符分隔的文件
# read.delim 默认就是制表符分隔
# 虽然名字是csv，但实际上是制表符分隔
otu_raw <- read.csv("archaea_data_rarefied_2913.csv", header = TRUE,row.names = 1,stringsAsFactors = FALSE,check.names = FALSE)
# 2. 读取分组信息
group_raw <- read.csv("group.csv", header=T, stringsAsFactors = FALSE)

# 3. 数据转置：Vegan包要求 行=样本，列=物种
otu_t <- as.data.frame(t(otu_raw))

# ============================================================
# Step 3: 构建完整的元数据 (Metadata) —— 关键步骤！！！
# ============================================================
# 你提供的数据只有 J17, J43 等编号。
# 你需要在这里根据你的实际采样记录，定义每个编号属于哪个 属(Genus) 和 哪个支系(Clade)。
# 下面我是根据常见珊瑚做的一个【示例】，请务必修改为你的真实信息！

# 定义映射关系 (请根据你的实验记录修改这里！)
# 假设：
# Robusta: J17, J19, J37, J43 (比如 Platygyra, Favites)
# Complexa: J54, J61, J66, J77 (比如 Acropora, Porites)

metadata <- group_raw %>%
  mutate(
    # 1. 定义属 (Genus) - 用于形状区分
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
    
    # 2. 定义进化支系 (Clade) - 用于颜色和椭圆区分 (核心假设)
    Clade = case_when(
      group %in% c("J54", "J61", "J66") ~ "Robusta",
      group %in% c("J17", "J19", "J37", "J43") ~ "Complexa",
      TRUE ~ "Unknown"
    )
  )

# 确保样品顺序一致
# 找到 otu_t 中存在的样品，并对齐 metadata
samples_intersect <- intersect(rownames(otu_t), metadata$sample)
otu_t <- otu_t[samples_intersect, ]
metadata <- metadata[match(samples_intersect, metadata$sample), ]

# ============================================================
# Step 4: 计算距离与 PCoA
# ============================================================

# 1. 计算 Bray-Curtis 距离
dist_mat <- vegdist(otu_t, method = "bray")

# 2. 进行 PCoA (PCoA = MDS with metric distance)
pcoa <- cmdscale(dist_mat, k = 2, eig = TRUE)

# 3. 提取坐标用于绘图
pcoa_points <- as.data.frame(pcoa$points)
colnames(pcoa_points) <- c("PCoA1", "PCoA2")
pcoa_data <- cbind(pcoa_points, metadata) # 合并分组信息

# 4. 计算轴解释率 (用于坐标轴标签)
eig_percent <- round(pcoa$eig / sum(pcoa$eig) * 100, 1)

# ============================================================
# Step 5: 统计检验 (PERMANOVA) - 这一步决定你的图是否具有统计意义
# ============================================================

# 检验 Clade (Robusta vs Complexa) 的差异
adonis_res <- adonis2(dist_mat ~ Clade, data = metadata, permutations = 9999)
R2_val <- round(adonis_res$R2[1], 3)
P_val <- adonis_res$`Pr(>F)`[1]

# 生成图上的标注文字
stat_label <- paste0("PERMANOVA\nRobusta vs Complexa\nR2 = ", R2_val, "\nP = ", P_val)

# ============================================================

  # ============================================================
# Step 6: 绘制顶刊风格 PCoA 图 (NC / iMeta Style)
# ============================================================

# 定义8种不同的形状
# 形状编号: 0=正方形, 1=圆形, 2=三角形, 3=十字, 4=X形, 5=菱形, 6=倒三角, 7=方加圆, 8=星形, 
#9=菱形加十字, 10=圆加十字, 11=三角加十字, 12=方加十字, 13=圆加点, 14=方加点, 15=实心方, 16=实心圆, 17=实心三角, 18=实心菱形
custom_shapes <- c(13, 8, 1, 2, 4, 0, 7)  # 7种不同的形状

p <- ggplot(pcoa_data, aes(x = PCoA1, y = PCoA2)) +
  
  # 1. 添加置信椭圆 (背景层)
  stat_ellipse(aes(fill = Clade, color = Clade), 
               type = "t", level = 0.95, geom = "polygon", alpha = 0.15, show.legend = FALSE) +
  stat_ellipse(aes(color = Clade), 
               type = "t", level = 0.95, geom = "path", linetype = 2, alpha = 0.8, show.legend = FALSE) +
  
  # 2. 添加散点：颜色代表Clade，形状代表Genus
  geom_point(aes(color = Clade, shape = Genus), size = 4, alpha = 0.9) +
  
  # 3. 配色和形状方案
  scale_color_npg() + 
  scale_fill_npg() +
  scale_shape_manual(values = custom_shapes) +  # 使用自定义形状
  
  # 4. 坐标轴标签
  labs(x = paste0("PCoA1 (", eig_percent[1], "%)"),
       y = paste0("PCoA2 (", eig_percent[2], "%)"),
       title = "Structure of Endolithic Archaeal Communities",
       subtitle = "Differentiation by Host Evolutionary Lineage") +
  
  # 5. 添加统计检验结果注释
  annotate("text", x = max(pcoa_data$PCoA1), y = max(pcoa_data$PCoA2), 
           label = stat_label, hjust = 1, vjust = 1, size = 4, fontface = "italic") +
  
  # 6. 主题美化
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )


# ============================================================
# Step 7: 输出结果
# ============================================================
print(p)

# 保存为高清PDF (推荐)
ggsave("PCoA_Robusta_vs_Complexa.pdf", p, width = 8, height = 5)
ggsave("PCoA_Robusta_vs_Complexa.png", p, width = 8, height = 5, dpi = 300)


