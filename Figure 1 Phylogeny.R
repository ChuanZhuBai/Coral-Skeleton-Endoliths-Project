if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")
BiocManager::install("treeio")
install.packages("ggplot2")
install.packages("ape")

# ============================================================
# Step 1: 加载必要的包
# ============================================================
library(ggplot2)
library(ggtree)
library(treeio)
library(ape)

# ============================================================
# Step 2: 读取并处理树文件
# ============================================================

# 1. 定义你的 Newick 字符串
nwk_text <- "(((Acropora_sp.:0.05349219,Astreopora_sp.:0.00495695)0.7810:0.02332794,(Goniopora_sp.:0.02048619,Porites_sp.:0.02027785)0.9190:0.03457251)0.9980:0.08608893,(Dipsastraea_sp.:0.00152518,(Favites_sp.:0.00170417,Platygyra_sp.:0.00170300)0.6010:0.00188739)1.0000:0.11920196,Nematostella_vectensis:0.17800233);"

# 2. 读取树
tree <- read.tree(text = nwk_text)

# 3. 指定外群进行定根 (Rooting)
# 这一步非常重要，确保进化方向正确
tree_rooted <- root(tree, outgroup = "Nematostella_vectensis", resolve.root = TRUE)

# ============================================================
# Step 3: 定义分组信息 (Robusta vs Complexa)
# ============================================================

# 创建一个分组列表，用于给树枝上色
groupInfo <- list(
  Complexa = c("Acropora_sp.", "Astreopora_sp.", "Goniopora_sp.", "Porites_sp."),
  Robusta = c("Dipsastraea_sp.", "Favites_sp.", "Platygyra_sp."),
  Outgroup = c("Nematostella_vectensis")
)

# 将分组信息整合到树对象中
tree_grouped <- groupOTU(tree_rooted, groupInfo)

# ============================================================
# Step 4: 顶刊风格绘图 (NC / iMeta Style)
# ============================================================

# 设置配色 (保持与 PCoA 一致: Complexa=红, Robusta=蓝)
# 0(默认黑线), 1(Complexa), 2(Outgroup), 3(Robusta) - 顺序可能根据字母排序变动，下面手动指定
cols <- c("black", "#E64B35", "grey50", "#4DBBD5") 
# 注意：#E64B35 是 Nature红色，#4DBBD5 是 Nature蓝色

p <- ggtree(tree_grouped, aes(color = group), size = 1) + # size调整线条粗细
  
  # 1. 设置树枝颜色
  scale_color_manual(values = c("Complexa" = "#E64B35", # 红色
                                "Robusta" = "#4DBBD5",  # 蓝色
                                "Outgroup" = "grey60",  # 灰色
                                "0" = "black")) +       # 骨架颜色
  
  # 2. 处理末端物种名 (Tip Labels)
  geom_tiplab(
    aes(label = gsub("_", " ", label)), # 把下划线替换为空格
    size = 4.5,                         # 字体大小
    fontface = "italic",                # 斜体
    offset = 0.005,                     # 标签与树枝的距离
    color = "black"                     # 文字颜色始终为黑色
  ) +
  
  # 3. 标记节点置信度 (Bootstrap values)
  # 只显示大于 0.7 (即70%) 的值，避免图太乱
  geom_nodepoint(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) > 0.7), 
                 size = 3, shape = 21, fill = "white", color = "black") +
  
  # 4. 添加右侧的分类色条 (Clade Label Bars)
  # 需要找到每个分支的最近共同祖先(MRCA)节点号
  # 这里使用 ggtree 的自动查找功能，或者手动指定
  # Complexa 条
  geom_cladelabel(node = MRCA(tree_rooted, c("Acropora_sp.", "Porites_sp.")), 
                  label = "Complexa", 
                  align = TRUE, 
                  offset = 0.06,    # 色条离树的距离
                  offset.text = 0.01, # 文字离色条的距离
                  barsize = 2,      # 色条宽度
                  fontsize = 5, 
                  color = "#E64B35") +
  
  # Robusta 条
  geom_cladelabel(node = MRCA(tree_rooted, c("Dipsastraea_sp.", "Platygyra_sp.")), 
                  label = "Robusta", 
                  align = TRUE, 
                  offset = 0.06, 
                  offset.text = 0.01,
                  barsize = 2, 
                  fontsize = 5, 
                  color = "#4DBBD5") +
  
  # 5. 添加比例尺 (Scale Bar)
  geom_treescale(x = 0, y = 1, width = 0.05, fontsize = 4, linesize = 1) +
  
  # 6. 调整显示范围，防止标签文字溢出
  xlim(0, 0.35) +
  
  # 7. 去除图例 (如果不想要group的图例)
  theme(legend.position = "none")

# ============================================================
# Step 5: 保存结果
# ============================================================
print(p)

# 保存为 PDF (矢量图，适合投稿)
ggsave("Figure1_Host_Phylogeny.pdf", p, width = 8, height = 4)
# 保存为 PNG (高清位图)
ggsave("Figure1_Host_Phylogeny.png", p, width = 8, height = 6, dpi = 300)

