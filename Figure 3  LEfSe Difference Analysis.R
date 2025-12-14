# ============================================================
# 定义清洗函数：删除行和为 0 的 ASV
# ============================================================
remove_zero_asvs <- function(input_file, output_file) {
  
  # 1. 读取数据
  # check.names=FALSE 保持列名原样 (防止 J17-A 变成 J17.A)
  otu <- read.csv(input_file, header=TRUE, row.names=1, check.names=FALSE)
  
  # 记录原始数量
  original_rows <- nrow(otu)
  
  # 2. 计算每一行的和 (rowSums)
  # 筛选出和大于 0 的行
  otu_clean <- otu[rowSums(otu) > 0, ]
  
  # 记录清洗后数量
  clean_rows <- nrow(otu_clean)
  removed_count <- original_rows - clean_rows
  
  # 3. 打印报告
  print(paste("正在处理文件:", input_file))
  print(paste("  - 原始 ASV 数量:", original_rows))
  print(paste("  - 删除的全零 ASV:", removed_count))
  print(paste("  - 剩余 ASV 数量:", clean_rows))
  
  # 4. 保存文件
  write.csv(otu_clean, output_file, quote=FALSE)
  print(paste("  - 已保存清洗后的文件至:", output_file))
  print("------------------------------------------------------")
}

# ============================================================
# 执行清洗
# ============================================================

# 1. 清洗细菌数据
remove_zero_asvs(
  input_file = "bacteria_data_rarefied_9480.csv", 
  output_file = "bacteria_data_rarefied_9480_no_zeros.csv"
)

# 2. 清洗古菌数据
remove_zero_asvs(
  input_file = "archaea_data_rarefied_2913.csv", 
  output_file = "archaea_data_rarefied_2913_no_zeros.csv"
)

########################################
#######################################
# ============================================================
# Step 1: 加载包
# ============================================================
library(dplyr)
library(ggplot2)
library(MASS) # 用于 LDA
library(ggpubr)
library(cowplot)
# ============================================================
# 辅助函数：去除前后空白并确保名字一致
# ============================================================
trim_names <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x
}

# ============================================================
# Step 2: 定义通用的 LEfSe 分析函数（终极修正版）
# ============================================================
run_lefse_analysis <- function(otu_file, tax_file, group_file, domain_name, lda_threshold = 1.0, verbose = TRUE) {
  
  # --- 1. 数据读取与预处理 ---
  otu <- read.csv(otu_file, header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  tax <- read.csv(tax_file, header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  group <- read.csv(group_file, header = TRUE, stringsAsFactors = FALSE)
  
  # trim names
  rownames(otu) <- trim_names(rownames(otu))
  rownames(tax) <- trim_names(rownames(tax))
  group$sample <- trim_names(group$sample)
  group$group <- trim_names(group$group)
  
  # Metadata
  metadata <- group %>%
    mutate(Clade = case_when(
      group %in% c("J54", "J61", "J66") ~ "Robusta",
      group %in% c("J17", "J19", "J37", "J43") ~ "Complexa",
      TRUE ~ "Unknown"
    )) %>%
    filter(Clade != "Unknown")
  
  # --- 2. 聚合到 Genus 水平 ---
  if(!"Genus" %in% colnames(tax)) stop("Tax file must contain a column named 'Genus'.")
  tax$Genus[is.na(tax$Genus) | tax$Genus == ""] <- "Unassigned"
  clean_genus <- gsub("^g__", "", tax$Genus)
  clean_genus <- trim_names(clean_genus)
  clean_genus[clean_genus == "" | clean_genus == "uncultured" | clean_genus == "Unassigned"] <- "Unclassified"
  
  common_ids <- intersect(rownames(otu), rownames(tax))
  otu_sub <- otu[common_ids, , drop = FALSE]
  genus_vec <- clean_genus[match(common_ids, rownames(tax))]
  
  # aggregate
  agg_df <- aggregate(otu_sub, by = list(Genus = genus_vec), FUN = sum)
  rownames(agg_df) <- agg_df$Genus
  otu_genus <- agg_df[, -1, drop = FALSE]
  
  # 强制数值化
  otu_genus <- as.data.frame(lapply(otu_genus, function(x) as.numeric(as.character(x))))
  rownames(otu_genus) <- agg_df$Genus # 补回行名
  otu_genus[is.na(otu_genus)] <- 0
  
  # 转化为相对丰度 (百分比 0-100)
  otu_rel <- sweep(otu_genus, 2, colSums(otu_genus), "/") * 100
  otu_rel[is.na(otu_rel)] <- 0
  
  # --- 3. 差异分析 ---
  data_t <- as.data.frame(t(otu_rel))
  data_t$Sample <- rownames(data_t)
  data_merged <- merge(data_t, metadata, by.x = "Sample", by.y = "sample")
  
  results <- data.frame()
  all_genera <- rownames(otu_genus)
  
  for (gen in all_genera) {
    if (!gen %in% colnames(data_merged)) next
    
    vec <- data_merged[[gen]]
    if (sum(vec > 0) < 3) next # 只有不到3个样本有该菌，跳过
    
    # 过滤低丰度 (平均相对丰度 < 0.01%)
    if (mean(vec, na.rm=TRUE) < 0.01) {
      # if (verbose) message(sprintf("Skipping '%s' (low abundance)", gen))
      next
    }
    
    # Wilcoxon
    wt_res <- tryCatch(wilcox.test(vec ~ data_merged$Clade), error = function(e) NULL)
    
    if (!is.null(wt_res) && !is.na(wt_res$p.value) && wt_res$p.value < 0.05) {
      
      means <- tapply(vec, data_merged$Clade, mean, na.rm = TRUE)
      
      # LDA-like Effect Size (Log10 mean ratio)
      # 加一个极小值防止除以0
      effect_size <- log10((means["Complexa"] + 1e-6) / (means["Robusta"] + 1e-6))
      enriched <- ifelse(means["Complexa"] > means["Robusta"], "Complexa", "Robusta")
      
      results <- rbind(results, data.frame(
        Genus = gen,
        P_value = wt_res$p.value,
        LDA = as.numeric(effect_size),
        Enriched = enriched,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # --- 4. 结果筛选 ---
  # 如果没有显著结果，直接返回空表，不报错
  if (nrow(results) == 0) {
    message(paste("No significant taxa found for", domain_name))
    return(data.frame())
  }
  
  # 筛选阈值
  final_res <- results %>%
    filter(abs(LDA) > lda_threshold) %>%
    arrange(desc(abs(LDA)))
  
  if (nrow(final_res) == 0) {
    message(paste("No taxa passed LDA threshold for", domain_name))
    return(data.frame())
  }
  
  # 添加 Domain (确保有行才添加)
  final_res$Domain <- domain_name
  final_res$LDA_plot <- abs(final_res$LDA) # 用于绘图
  
  return(final_res)
}

# 绘图函数（与之前版本兼容）
plot_lefse_bar <- function(res_data, title_text) {
  if (is.null(res_data) || nrow(res_data) == 0) {
    return(ggplot() + annotate("text", x = 1, y = 1, label = "No significant taxa found") + theme_void())
  }
  
  plot_data <- res_data %>%
    mutate(Plot_LDA = ifelse(Enriched == "Robusta", -LDA_plot, LDA_plot))
  
  # 为保证顺序稳定，按 Plot_LDA 排序
  plot_data$Genus <- factor(plot_data$Genus, levels = plot_data$Genus[order(plot_data$Plot_LDA)])
  
  p <- ggplot(plot_data, aes(x = Genus, y = Plot_LDA, fill = Enriched)) +
    geom_bar(stat = "identity", width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = c("Complexa" = "#E64B35", "Robusta" = "#4DBBD5")) +
    labs(x = NULL, y = "LDA-like Score (log10 ratio)", title = title_text) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 10, face = "italic", color = "black"),
      axis.text.x = element_text(size = 10),
      panel.grid = element_blank(),
      legend.position = "top",
      legend.title = element_blank()
    )
  return(p)
}


# 1. 细菌分析 (阈值 1.0)
lefse_bac <- run_lefse_analysis(
  "bacteria_data_rarefied_9480_no_zeros.csv", 
  "bacteria_tax_table.csv", 
  "group.csv", 
  "Bacteria", 
  lda_threshold = 4.0 # 建议先用 1.0，如果太多再加到 2.0
)
p_lefse_bac <- plot_lefse_bar(lefse_bac, "Bacterial Biomarkers")

print(p_lefse_bac)

# 2. 古菌分析 (降低阈值试一试)
# 注意：如果古菌真的没差异，这里可能还是空表
lefse_arc <- run_lefse_analysis(
  "archaea_data_rarefied_2913_no_zeros.csv", 
  "archaea_tax_table.csv", 
  "group.csv", 
  "Archaea", 
  lda_threshold = 2 # 这里的阈值是 Log10 倍数，0.5 约等于 3 倍差异
)
p_lefse_arc <- plot_lefse_bar(lefse_arc, "Archaeal Biomarkers")

print(p_lefse_arc)
# 3. 拼图
fig3b <- plot_grid(p_lefse_bac, p_lefse_arc, ncol = 2, align = "h", labels = c("B1", "B2"))
print(fig3b)
ggsave("Figure3B_LEfSe_Barplot.pdf", fig3b, width = 12, height = 6)


###########################################################
library(dplyr)
library(ggplot2)

# ============================================================
# 改进版绘图函数：支持筛选 Top N
# ============================================================
plot_lefse_bar_optimized <- function(res_data, title_text, top_n = 10) {
  
  if (is.null(res_data) || nrow(res_data) == 0) {
    return(ggplot() + annotate("text", x = 1, y = 1, label = "No significant taxa found") + theme_void())
  }
  
  # --- 1. 筛选 Top N ---
  # 分组筛选：每个 Clade 只保留 LDA 绝对值最大的前 top_n 个
  # 这样防止 Complexa 太多把 Robusta 挤没了
  plot_data_filtered <- res_data %>%
    group_by(Enriched) %>%
    slice_max(order_by = abs(LDA), n = top_n) %>%
    ungroup()
  
  # --- 2. 数据处理用于绘图 ---
  # 为了让 Robusta 朝左，Complexa 朝右
  plot_data <- plot_data_filtered %>%
    mutate(Plot_LDA = ifelse(Enriched == "Robusta", -abs(LDA), abs(LDA)))
  
  # 排序：让柱子按长短排列，看起来更整齐
  plot_data$Genus <- factor(plot_data$Genus, levels = plot_data$Genus[order(plot_data$Plot_LDA)])
  
  # --- 3. 绘图 ---
  p <- ggplot(plot_data, aes(x = Genus, y = Plot_LDA, fill = Enriched)) +
    geom_bar(stat = "identity", width = 0.8, color = "black", size = 0.2) + # 加个细黑边更有质感
    coord_flip() + # 保持横向条形图（这是展示长物种名的最佳方式）
    
    # 配色
    scale_fill_manual(values = c("Complexa" = "#E64B35", "Robusta" = "#4DBBD5")) +
    
    # 调整 Y 轴范围和标签 (去掉负号)
    scale_y_continuous(labels = function(x) abs(x)) + 
    
    # 标题和标签
    labs(x = NULL, y = "LDA Score (log10)", title = title_text,
         subtitle = paste0("Top ", top_n, " Enriched Genera per Clade")) +
    
    # 主题美化
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 11, face = "italic", color = "black"), # 物种名斜体
      axis.text.x = element_text(size = 10, color = "black"),
      axis.title.x = element_text(size = 12, face = "bold"),
      panel.grid.major.y = element_blank(), # 去掉横向网格
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40")
    )
  
  return(p)
}

# ============================================================
# 执行绘图 (使用新的函数)
# ============================================================

# 1. 细菌 (假设 lefse_bac 已经算好了)
# 设置 top_n = 15，只显示 Complexa 前15个和 Robusta 的所有(因为它不到15个)
p_lefse_bac_opt <- plot_lefse_bar_optimized(lefse_bac, "Bacterial Biomarkers", top_n = 10)

# 2. 古菌 (假设 lefse_arc 已经算好了)
p_lefse_arc_opt <- plot_lefse_bar_optimized(lefse_arc, "Archaeal Biomarkers", top_n = 10)

print(p_lefse_bac_opt)
# 3. 组合
library(cowplot)
fig3b_final <- plot_grid(p_lefse_bac_opt, p_lefse_arc_opt, ncol = 2, align = "h", labels = c("B1", "B2"))

print(fig3b_final)

# 保存时调整长宽比
ggsave("Figure3B_LEfSe_Top15.pdf", fig3b_final, width = 12, height = 6)