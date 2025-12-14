# ============================================================
# Step 1: 加载包
# ============================================================
library(igraph)
library(psych) # 用于计算相关性矩阵和P值
library(dplyr)
library(tibble)

# ============================================================
# Step 2: 数据读取与预处理
# ============================================================

# 1. 读取数据 (使用去除全零后的抽平数据)
bac_otu <- read.csv("bacteria_data_rarefied_9480_no_zeros.csv", row.names=1, check.names=F)
arc_otu <- read.csv("archaea_data_rarefied_2913_no_zeros.csv", row.names=1, check.names=F)

bac_tax <- read.csv("bacteria_tax_table.csv", row.names=1, stringsAsFactors=F)
arc_tax <- read.csv("archaea_tax_table.csv", row.names=1, stringsAsFactors=F)

group <- read.csv("group.csv", stringsAsFactors=F)

# 2. 转换为相对丰度 (Relative Abundance)
# 这一步非常重要，让细菌和古菌在同一量级上比较
bac_rel <- sweep(bac_otu, 2, colSums(bac_otu), "/") 
arc_rel <- sweep(arc_otu, 2, colSums(arc_otu), "/") 

# 3. 补充 Taxonomy 信息 (添加 Domain 列，用于画图形状)
# 确保 Tax 表包含 Phylum
bac_tax$Domain <- "Bacteria"
arc_tax$Domain <- "Archaea"

# 清洗 Phylum 名字 (去掉 p__ 前缀)
bac_tax$Phylum <- gsub("^p__", "", bac_tax$Phylum)
arc_tax$Phylum <- gsub("^p__", "", arc_tax$Phylum)
# 处理空值
bac_tax$Phylum[bac_tax$Phylum == "" | is.na(bac_tax$Phylum)] <- "Unassigned"
arc_tax$Phylum[arc_tax$Phylum == "" | is.na(arc_tax$Phylum)] <- "Unassigned"

# 4. 合并数据
# (1) 合并 OTU 表 (行合并，前提是列名Sample一致)
# 先检查列名是否对齐
common_samples <- intersect(colnames(bac_rel), colnames(arc_rel))
if(length(common_samples) == 0) stop("细菌和古菌没有匹配的样品名！请检查 csv 文件。")

# 只保留公共样品并合并
total_otu <- rbind(bac_rel[, common_samples], arc_rel[, common_samples])

# (2) 合并 Tax 表
# 选取需要的列，保证列名一致
cols_keep <- c("Domain", "Phylum") # 你可以加 Class, Order 等
total_tax <- rbind(bac_tax[, cols_keep], arc_tax[, cols_keep])

# 确保 OTU 和 Tax 行名匹配
common_asvs <- intersect(rownames(total_otu), rownames(total_tax))
total_otu <- total_otu[common_asvs, ]
total_tax <- total_tax[common_asvs, ]

stopifnot( all(rownames(total_otu) == rownames(total_tax)) )

# ============================================================
# Step 3: 定义网络构建与分析函数
# ============================================================
# 1. 拆分数据
# 根据 metadata 找到 Robusta 和 Complexa 的样品名
meta_sub <- group %>% 
  mutate(Clade = case_when(
    group %in% c("J54", "J61", "J66") ~ "Robusta",
    group %in% c("J17", "J19", "J37", "J43") ~ "Complexa",
    TRUE ~ "Unknown"
  ))

samples_robusta <- meta_sub$sample[meta_sub$Clade == "Robusta"]
samples_complexa <- meta_sub$sample[meta_sub$Clade == "Complexa"]

# 提取子集 OTU 表
# 注意：取子集后，有些 ASV 可能全为 0，需要在函数内部过滤
otu_robusta <- total_otu[, intersect(colnames(total_otu), samples_robusta)]
otu_complexa <- total_otu[, intersect(colnames(total_otu), samples_complexa)]

# ============================================================
# 优化后的网络构建函数
# ============================================================
build_network_optimized <- function(otu_table, tax_table, group_name, r_threshold=0.6, p_threshold=0.05) {
  
  print(paste(">>> 正在处理分组:", group_name))
  
  # --- 1. 优化过滤策略 (Prevalence + Abundance) ---
  
  # 计算平均丰度
  mean_abund <- rowMeans(otu_table)
  # 计算出现率 (Prevalence): 在多少比例的样品中 > 0
  prevalence <- rowSums(otu_table > 0) / ncol(otu_table)
  
  # 设定阈值：
  # 策略：保留 (平均丰度 > 0.1%) OR (出现率 > 30%) 的 ASV
  # 这样既保留了高丰度物种，也保留了低丰度但稳定的核心物种
  keep_idx <- (mean_abund > 0.001) | (prevalence > 0.3)
  
  otu_filt <- otu_table[keep_idx, ]
  
  # --- 2. 检查过滤后的 细菌/古菌 比例 (关键步骤) ---
  # 从 tax_table 中提取留下的 ASV 的 Domain 信息
  kept_tax <- tax_table[rownames(otu_filt), ]
  domain_counts <- table(kept_tax$Domain)
  
  print("  [过滤报告]")
  print(paste("  - 原始 ASV:", nrow(otu_table)))
  print(paste("  - 过滤后 ASV:", nrow(otu_filt)))
  print("  - 过滤后各界节点数量:")
  print(domain_counts)
  
  # 如果古菌少于 3 个，可能无法形成跨界网络，给出警告
  if(is.na(domain_counts["Archaea"]) || domain_counts["Archaea"] < 3) {
    print("  警告: 过滤后古菌节点过少，跨界分析可能受限！")
  }
  
  if(nrow(otu_filt) < 5) {
    print("  错误: 符合条件的 ASV 太少，无法构建网络")
    return(NULL)
  }
  
  # --- 3. 计算相关性 (Spearman) ---
  otu_t <- t(otu_filt)
  cor_res <- corr.test(otu_t, use="pairwise", method="spearman", adjust="fdr", ci=FALSE)
  r_matrix <- cor_res$r
  p_matrix <- cor_res$p
  
  diag(r_matrix) <- 0
  diag(p_matrix) <- 1
  
  # --- 4. 筛选边 ---
  adj_matrix <- r_matrix
  adj_matrix[abs(r_matrix) < r_threshold | p_matrix > p_threshold] <- 0
  
  # --- 5. 构建图对象 ---
  g <- graph_from_adjacency_matrix(adj_matrix, mode="undirected", weighted=TRUE, diag=FALSE)
  g <- delete.vertices(g, which(degree(g) == 0)) # 删除孤立点
  
  print(paste("  - 网络构建完成：节点数 =", vcount(g), ", 边数 =", ecount(g)))
  
  if(vcount(g) == 0) return(NULL)
  
  # --- 6. 添加节点属性 ---
  node_ids <- V(g)$name
  node_tax <- tax_table[node_ids, ]
  V(g)$Domain <- node_tax$Domain
  V(g)$Phylum <- node_tax$Phylum
  V(g)$Degree <- degree(g)
  
  # --- 7. 添加边属性 ---
  E(g)$Interaction <- ifelse(E(g)$weight > 0, "Positive", "Negative")
  
  # --- 8. 计算拓扑属性 (无权模式) ---
  g_topo <- delete_edge_attr(g, "weight")
  mod_cluster <- cluster_fast_greedy(g_topo)
  
  topo_stats <- data.frame(
    Group = group_name,
    Nodes = vcount(g),
    # 统计最终网络里还有多少个细菌和古菌
    Nodes_Bac = sum(V(g)$Domain == "Bacteria"),
    Nodes_Arc = sum(V(g)$Domain == "Archaea"),
    Edges = ecount(g),
    Avg_Degree = mean(degree(g)),
    Avg_Path_Length = mean_distance(g_topo, directed=FALSE),
    Clustering_Coeff = transitivity(g_topo, type="global"),
    Modularity = modularity(mod_cluster),
    Graph_Density = edge_density(g),
    Pos_Edges = sum(E(g)$weight > 0),
    Neg_Edges = sum(E(g)$weight < 0),
    Pos_Neg_Ratio = sum(E(g)$weight > 0) / (sum(E(g)$weight < 0) + 0.0001)
  )
  
  return(list(graph = g, stats = topo_stats))
}

# ============================================================
# 执行 (请确保之前的 otu_robusta 等数据已准备好)
# ============================================================
net_rob <- build_network_optimized(otu_robusta, total_tax, "Robusta", r_threshold=0.6, p_threshold=0.05)
net_com <- build_network_optimized(otu_complexa, total_tax, "Complexa", r_threshold=0.6, p_threshold=0.05)

# 导出 (同前)
if(!is.null(net_rob) & !is.null(net_com)) {
  topology_table <- rbind(net_rob$stats, net_com$stats)
  print(topology_table)
  write.csv(topology_table, "TableS4_Network_Topology.csv", row.names = FALSE)
  
  write_graph(net_rob$graph, "Network_Robusta.gml", format = "gml")
  write_graph(net_com$graph, "Network_Complexa.gml", format = "gml")
}


################################################################################
################################################################################
# ============================================================
# Step 1: 加载包
# ============================================================
library(igraph)
library(ggplot2)
library(dplyr)
library(ggrepel) # 用于添加不重叠的标签

# ============================================================
# Step 2: 定义计算 Zi 和 Pi 的函数
# ============================================================
# 这是一个核心算法函数
calc_zi_pi <- function(igraph_obj, group_name) {
  
  # 1. 再次确保计算模块化 (基于无权图以防报错)
  g_topo <- delete_edge_attr(igraph_obj, "weight")
  # 使用 fast_greedy 算法划分模块 (与你Table S4一致)
  wtc <- cluster_fast_greedy(g_topo)
  
  # 获取每个节点的模块归属
  membership <- membership(wtc)
  module_list <- unique(membership)
  
  # 获取邻接矩阵
  adj <- as.matrix(as_adjacency_matrix(g_topo))
  
  # 初始化结果表
  res <- data.frame(
    ID = V(igraph_obj)$name,
    Domain = V(igraph_obj)$Domain,
    Phylum = V(igraph_obj)$Phylum,
    Module = as.numeric(membership),
    Group = group_name,
    Zi = NA,
    Pi = NA
  )
  
  # 计算 degree (k)
  ki <- degree(g_topo)
  
  # --- 开始循环计算 ---
  for (i in 1:nrow(res)) {
    node_id <- res$ID[i]
    node_mod <- res$Module[i]
    
    # 1. 计算 Zi (Within-module connectivity)
    # 找到同模块的所有节点
    nodes_in_mod <- which(membership == node_mod)
    # 节点 i 在模块内的度 (kis)
    kis <- sum(adj[node_id, names(membership)[nodes_in_mod]])
    
    # 计算该模块内所有节点的平均度和标准差
    ks_all <- degree(g_topo, v = nodes_in_mod)
    ksi_bar <- mean(ks_all)
    ksi_sigma <- sd(ks_all)
    
    # 防止除以0 (如果模块太小)
    if (ksi_sigma == 0) {
      res$Zi[i] <- 0
    } else {
      res$Zi[i] <- (kis - ksi_bar) / ksi_sigma
    }
    
    # 2. 计算 Pi (Among-module connectivity)
    # 节点 i 的总度 (ki)
    ki_total <- ki[node_id]
    
    # 计算节点 i 与每个模块的连接数 (kit)
    sum_kit_sq <- 0
    for (m in module_list) {
      nodes_in_m <- which(membership == m)
      kit <- sum(adj[node_id, names(membership)[nodes_in_m]])
      sum_kit_sq <- sum_kit_sq + (kit / ki_total)^2
    }
    
    res$Pi[i] <- 1 - sum_kit_sq
  }
  
  # 定义角色 (Roles)
  res$Role <- case_when(
    res$Zi <= 2.5 & res$Pi <= 0.62 ~ "Peripherals",    # 边缘节点 (大多数)
    res$Zi <= 2.5 & res$Pi > 0.62  ~ "Connectors",     # 模块间连接者
    res$Zi > 2.5 & res$Pi <= 0.62  ~ "Module Hubs",    # 模块内枢纽
    res$Zi > 2.5 & res$Pi > 0.62  ~ "Network Hubs"     # 全网枢纽
  )
  
  return(res)
}

# ============================================================
# Step 3: 执行计算
# ============================================================
# 假设 net_rob 和 net_com 是你上一步生成的对象
# 如果没有对象，请先用 read_graph("Network_Robusta.gml", format="gml") 读取

# 1. 计算 Robusta
df_rob <- calc_zi_pi(net_rob$graph, "Robusta")

# 2. 计算 Complexa
df_com <- calc_zi_pi(net_com$graph, "Complexa")

# 3. 合并数据
plot_data <- rbind(df_rob, df_com)

# 设置因子水平，保证 Robusta 和 Complexa 顺序
plot_data$Group <- factor(plot_data$Group, levels = c("Complexa", "Robusta"))

# ============================================================
# Step 4: 绘制顶刊风格 Zipper Plot
# ============================================================

# 自定义形状: 细菌=圆(16), 古菌=三角(17)
# 自定义颜色: 细菌=灰色, 古菌=红色 (突出古菌)
# 或者按 Role 上色，按 Domain 定形状

p_zipper <- ggplot(plot_data, aes(x = Pi, y = Zi)) +
  
  # 1. 添加阈值线
  geom_hline(yintercept = 2.5, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0.62, linetype = "dashed", color = "grey50") +
  
  # 2. 散点图
  # 形状区分细菌古菌，颜色区分角色(或者直接把古菌标红)
  geom_point(aes(shape = Domain, fill = Domain, size = Zi), alpha = 0.8) +
  
  # 3. 样式设置
  scale_shape_manual(values = c("Bacteria" = 21, "Archaea" = 24)) + # 21圆, 24三角(可填充)
  scale_fill_manual(values = c("Bacteria" = "#8491B4", "Archaea" = "#E64B35")) + # 细菌灰蓝，古菌红
  scale_size_continuous(range = c(2, 6), guide = "none") + # 根据Zi大小调整点大小
  
  # 4. 分面展示
  facet_grid(. ~ Group) +
  
  # 5. 添加区域文字标签
  annotate("text", x = 0.1, y = 4, label = "Module Hubs", size = 3.5, fontface="bold", color="grey40") +
  annotate("text", x = 0.8, y = 4, label = "Network Hubs", size = 3.5, fontface="bold", color="grey40") +
  annotate("text", x = 0.8, y = 1, label = "Connectors", size = 3.5, fontface="bold", color="grey40") +
  annotate("text", x = 0.1, y = 1, label = "Peripherals", size = 3.5, fontface="bold", color="grey40") +
  
  # 6. 标记关键物种 (古菌 或 任何 Hubs)
  # 只标记属于 Module Hubs, Network Hubs, Connectors 的点，或者所有古菌
  geom_text_repel(data = subset(plot_data, Domain == "Archaea" | Role != "Peripherals"),
                  aes(label = ifelse(Domain == "Archaea", paste0(Phylum, "(Arc)"), Phylum)),
                  size = 3, max.overlaps = 20, box.padding = 0.5) +
  
  # 7. 坐标轴与主题
  labs(x = "Among-module connectivity (Pi)", 
       y = "Within-module connectivity (Zi)",
       title = "Keystone Species Analysis (Zi-Pi Plot)") +
  
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 12, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"), # 分面标题大小
    strip.background = element_rect(fill = "white"),
    legend.position = "bottom",
    legend.title = element_blank()
  )

# ============================================================
# Step 5: 保存与导出
# ============================================================
print(p_zipper)

ggsave("Figure4_PanelC_ZipperPlot.pdf", p_zipper, width = 12, height = 6)

# 导出关键物种名单 (放入 Table S5)
keystone_species <- plot_data %>% 
  filter(Role != "Peripherals") %>%
  select(Group, ID, Domain, Phylum, Role, Zi, Pi) %>%
  arrange(Group, desc(Zi))

write.csv(keystone_species, "TableS5_Keystone_Species.csv", row.names = FALSE)
print("关键物种名单已生成: TableS5_Keystone_Species.csv")











#####################################################################################
