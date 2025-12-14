install.packages(c("ggplot2", "ggalluvial", "dplyr", "reshape2", "RColorBrewer", "cowplot"))

# ============================================================
# Step 1: 加载必要的包
# ============================================================
library(ggplot2)
library(ggalluvial) # 核心包：用于画冲击图(连线堆叠图)
library(dplyr)
library(reshape2)
library(RColorBrewer) # 配色
library(cowplot)      # 拼图

# ============================================================
# Step 2: 定义通用的数据处理函数
# ============================================================
# 这个函数负责：读取数据 -> 清洗物种名 -> 按门合并 -> 筛选Top10 -> 按属取平均 -> 整理成绘图格式

process_abundance_data <- function(otu_file, tax_file, group_file, domain_name) {
  
  # 1. 读取数据
  otu <- read.csv(otu_file, header=T, row.names=1, check.names=F)
  tax <- read.csv(tax_file, header=T, row.names=1, check.names=F, stringsAsFactors=F)
  group <- read.csv(group_file, header=T, stringsAsFactors=F)
  
  # 2. 构建 Metadata (映射属和Clade)
  metadata <- group %>%
    mutate(
      Genus = case_when(
        group %in% c("J17") ~ "Astreopora",  
        group %in% c("J19") ~ "Goniopora", 
        group %in% c("J37") ~ "Porites",    
        group %in% c("J43") ~ "Acropora",     
        group %in% c("J54") ~ "Platygyra",   
        group %in% c("J61") ~ "Favites", 
        group %in% c("J66") ~ "Dipsastraea",    
        TRUE ~ "Unknown"
      ),
      Clade = case_when(
        group %in% c("J54", "J61", "J66") ~ "Robusta",
        group %in% c("J17", "J19", "J37", "J43") ~ "Complexa",
        TRUE ~ "Unknown"
      )
    ) %>%
    filter(Clade != "Unknown")
  
  # 3. 整理物种注释 (清洗 p__ 前缀)
  # 确保 Tax 表包含 Phylum 列
  if(!"Phylum" %in% colnames(tax)) stop("Taxonomy file must contain 'Phylum' column")
  
  # 提取 Phylum 并清洗
  # 处理 NA 和空值
  tax$Phylum[is.na(tax$Phylum) | tax$Phylum == ""] <- "Unassigned"
  # 去掉 "p__" 前缀 (如果有)
  clean_phylum <- gsub("^p__", "", tax$Phylum)
  # 处理一些特殊的未分类 (如 "uncultured")
  clean_phylum[clean_phylum == "uncultured"] <- "Unclassified"
  
  # 4. 合并 OTU 和 Tax
  # 转置 OTU 表 (行=ASV) 用于合并
  common_ids <- intersect(rownames(otu), rownames(tax))
  otu_sub <- otu[common_ids, ]
  phylum_vec <- clean_phylum[match(common_ids, rownames(tax))]
  
  # 按门汇总丰度 (aggregate)
  otu_phylum <- aggregate(otu_sub, by=list(Phylum=phylum_vec), FUN=sum)
  rownames(otu_phylum) <- otu_phylum$Phylum
  otu_phylum <- otu_phylum[, -1]
  
  # 5. 转化为相对丰度 (如果输入不是相对丰度的话)
  # 你的文件名叫 relative_abundance... csv，假设已经是相对丰度了
  # 为了保险，再归一化一次
  otu_rel <- sweep(otu_phylum, 2, colSums(otu_phylum), "/")
  
  # 6. 筛选 Top 11 门，其余归为 Others
  mean_abund <- rowMeans(otu_rel)
  top_n <- 11
  top_taxa <- names(sort(mean_abund, decreasing=TRUE)[1:top_n])
  
  # 构建最终表
  final_tab <- otu_rel[top_taxa, ]
  others <- 1 - colSums(final_tab)
  final_tab <- rbind(final_tab, Others=others)
  
  # 7. 按 "珊瑚属 (Genus)" 计算平均丰度
  # 转置为 (行=Sample, 列=Phylum)
  tab_t <- as.data.frame(t(final_tab))
  tab_t$Sample <- rownames(tab_t)
  
  # 合并 Metadata
  tab_merged <- merge(tab_t, metadata, by.x="Sample", by.y="sample")
  
  # 按 Genus 和 Clade 分组求均值
  # 使用 reshape2::melt 先变长，再聚合
  tab_long <- melt(tab_merged, id.vars=c("Sample", "group", "Genus", "Clade"), 
                   variable.name="Phylum", value.name="Abundance")
  
  # 计算均值 (每个属展示一个柱子)
  plot_data <- tab_long %>%
    group_by(Clade, Genus, Phylum) %>%
    summarise(MeanAbundance = mean(Abundance), .groups="drop")
  
  # 8. 设置因子水平 (排序)
  # Phylum 排序: Others 在最下面，其他的按丰度排序
  plot_data$Phylum <- factor(plot_data$Phylum, levels = c("Others", rev(top_taxa)))
  
  # Genus 排序: 
  # Complexa 在左 (Acropora, Astreopora, Goniopora, Porites)
  # Robusta 在右 (Dipsastraea, Favites, Platygyra)
  genus_order <- c("Acropora", "Astreopora", "Goniopora", "Porites", 
                   "Dipsastraea", "Favites", "Platygyra")
  plot_data$Genus <- factor(plot_data$Genus, levels = genus_order)
  
  return(plot_data)
}

# ============================================================
# Step 3: 定义绘图函数 (Alluvial Plot)
# ============================================================

plot_alluvial <- function(data, title_text) {
  
  # 定义颜色板 (Set3 比较柔和，适合物种堆叠)
  # 扩展颜色板以防物种超过12个
  my_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique(data$Phylum)))
  # 强制让 "Others" 变成灰色
  names(my_colors) <- levels(data$Phylum)
  my_colors["Others"] <- "grey90"
  
  p <- ggplot(data, aes(x = Genus, y = MeanAbundance, 
                        fill = Phylum, stratum = Phylum, alluvium = Phylum)) +
    
    # 1. 绘制连线 (Flow)
    geom_flow(alpha = 0.3, curve_type = "linear", width = 0.7) +
    
    # 2. 绘制柱子 (Stratum)
    geom_stratum(width = 0.7, color = "white", size = 0.1) +
    
    # 3. 分面：将 Robusta 和 Complexa 分开，但保持连贯性
    facet_grid(~Clade, scales = "free_x", space = "free_x") +
    
    # 4. 颜色与标签
    scale_fill_manual(values = my_colors) +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    labs(x = NULL, y = "Relative Abundance", title = title_text) +
    
    # 5. 主题美化
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, face = "italic"),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey95"),
      strip.text = element_text(face = "bold", size = 12),
      legend.title = element_blank(),
      legend.position = "right"
    )
  
  return(p)
}

# ============================================================
# Step 4: 执行分析与绘图
# ============================================================

# --- 1. 处理细菌数据 ---
bacteria_data <- process_abundance_data(
  otu_file = "relative_abundance_bacteria_rarefied_9480.csv",
  tax_file = "bacteria_tax_table.csv",
  group_file = "group.csv",
  domain_name = "Bacteria"
)

# 绘制细菌图
p_bac <- plot_alluvial(bacteria_data, "Bacterial Composition (Phylum Level)")

print(p_bac)
# --- 2. 处理古菌数据 ---
archaea_data <- process_abundance_data(
  otu_file = "relative_abundance_archaea_rarefied_2913.csv",
  tax_file = "archaea_tax_table.csv",
  group_file = "group.csv",
  domain_name = "Archaea"
)

# 绘制古菌图
p_arc <- plot_alluvial(archaea_data, "Archaeal Composition (Phylum Level)")

print(p_arc)
# ============================================================
# Step 5: 组合与保存
# ============================================================

# 组合图片 (上下排列)
final_fig3a <- plot_grid(p_bac, p_arc, ncol = 2, align = "v", labels = c("A1", "A2"))

print(final_fig3a)

# 保存
ggsave("Figure3A_Alluvial_Composition.pdf", final_fig3a, width = 10, height = 5)
ggsave("Figure3A_Alluvial_Composition.png", final_fig3a, width = 10, height = 5, dpi = 300)



#########################################################################################
#Complexa和Robusta大类上的连线堆叠图
###########################################################################################
# ============================================================
# Step 1: 加载必要的包
# ============================================================
library(ggplot2)
library(ggalluvial)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(cowplot)

# ============================================================
# Step 2: 定义按 Clade 汇总数据的函数
# ============================================================
process_data_by_clade <- function(otu_file, tax_file, group_file) {
  
  # 1. 读取数据
  otu <- read.csv(otu_file, header=T, row.names=1, check.names=F)
  tax <- read.csv(tax_file, header=T, row.names=1, check.names=F, stringsAsFactors=F)
  group <- read.csv(group_file, header=T, stringsAsFactors=F)
  
  # 2. 构建 Metadata (只关注 Clade)
  metadata <- group %>%
    mutate(Clade = case_when(
      group %in% c("J54", "J61", "J66") ~ "Robusta",
      group %in% c("J17", "J19", "J37", "J43") ~ "Complexa",
      TRUE ~ "Unknown"
    )) %>%
    filter(Clade != "Unknown")
  
  # 3. 整理物种注释 (清洗 Phylum)
  if(!"Phylum" %in% colnames(tax)) stop("Taxonomy file must contain 'Phylum' column")
  tax$Phylum[is.na(tax$Phylum) | tax$Phylum == ""] <- "Unassigned"
  clean_phylum <- gsub("^p__", "", tax$Phylum)
  clean_phylum[clean_phylum == "uncultured"] <- "Unclassified"
  
  # 4. 合并 OTU 和 Tax，并按门汇总
  common_ids <- intersect(rownames(otu), rownames(tax))
  otu_sub <- otu[common_ids, ]
  phylum_vec <- clean_phylum[match(common_ids, rownames(tax))]
  
  otu_phylum <- aggregate(otu_sub, by=list(Phylum=phylum_vec), FUN=sum)
  rownames(otu_phylum) <- otu_phylum$Phylum
  otu_phylum <- otu_phylum[, -1]
  
  # 5. 归一化为相对丰度
  otu_rel <- sweep(otu_phylum, 2, colSums(otu_phylum), "/")
  
  # 6. 筛选 Top 10 门
  mean_abund <- rowMeans(otu_rel)
  top_n <- 11
  top_taxa <- names(sort(mean_abund, decreasing=TRUE)[1:top_n])
  
  final_tab <- otu_rel[top_taxa, ]
  others <- 1 - colSums(final_tab)
  final_tab <- rbind(final_tab, Others=others)
  
  # 7. --- 关键修改：按 Clade 分组取平均 ---
  tab_t <- as.data.frame(t(final_tab))
  tab_t$Sample <- rownames(tab_t)
  
  tab_merged <- merge(tab_t, metadata, by.x="Sample", by.y="sample")
  
  tab_long <- melt(tab_merged, id.vars=c("Sample", "group", "Clade"), 
                   variable.name="Phylum", value.name="Abundance")
  
  # 聚合计算均值 (Group by Clade and Phylum)
  plot_data <- tab_long %>%
    group_by(Clade, Phylum) %>%
    summarise(MeanAbundance = mean(Abundance), .groups="drop")
  
  # 8. 设置因子水平 (排序)
  # Phylum: Others 在最下面
  plot_data$Phylum <- factor(plot_data$Phylum, levels = c("Others", rev(top_taxa)))
  
  # Clade: Complexa 在左，Robusta 在右 (符合进化顺序)
  plot_data$Clade <- factor(plot_data$Clade, levels = c("Complexa", "Robusta"))
  
  return(plot_data)
}

# ============================================================
# Step 3: 定义绘图函数 (Clade Comparision)
# ============================================================
plot_clade_alluvial <- function(data, title_text) {
  
  # 动态生成颜色
  n_colors <- length(unique(data$Phylum))
  my_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_colors)
  names(my_colors) <- levels(data$Phylum)
  my_colors["Others"] <- "grey90"
  
  p <- ggplot(data, aes(x = Clade, y = MeanAbundance, 
                        fill = Phylum, stratum = Phylum, alluvium = Phylum)) +
    
    # 1. 绘制连线 (Flow)
    # alpha=0.6 让连线稍微明显一点
    geom_flow(alpha = 0.4, curve_type = "sigmoid", width = 0.7) +
    
    # 2. 绘制柱子 (Stratum)
    geom_stratum(width = 0.7, color = "white", size = 0.1) +
    
    # 3. 颜色与标签
    scale_fill_manual(values = my_colors) +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    labs(x = NULL, y = "Mean Relative Abundance", title = title_text) +
    
    # 4. 主题
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 14, face = "bold", color = "black"), # 加大字体
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
  
  return(p)
}

# ============================================================
# Step 4: 执行与组合
# ============================================================

# 1. 细菌
df_bac <- process_data_by_clade(
  otu_file = "relative_abundance_bacteria_rarefied_9480.csv",
  tax_file = "bacteria_tax_table.csv",
  group_file = "group.csv"
)
p_bac_clade <- plot_clade_alluvial(df_bac, "Bacterial Community (Phylum)")

print(p_bac_clade)
# 2. 古菌
df_arc <- process_data_by_clade(
  otu_file = "relative_abundance_archaea_rarefied_2913.csv",
  tax_file = "archaea_tax_table.csv",
  group_file = "group.csv"
)
p_arc_clade <- plot_clade_alluvial(df_arc, "Archaeal Community (Phylum)")

print(p_arc_clade)
# 3. 拼图 (左右排列，直观对比)
final_clade_comp <- plot_grid(p_bac_clade, p_arc_clade, ncol = 2, align = "h", labels = c("A", "B"))

print(final_clade_comp)

# 保存
ggsave("Figure3_Clade_Comparison_Alluvial.pdf", final_clade_comp, width = 10, height = 6)
ggsave("Figure3_Clade_Comparison_Alluvial.png", final_clade_comp, width = 12, height = 6, dpi = 300)