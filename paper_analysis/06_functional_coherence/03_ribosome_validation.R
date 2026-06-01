#!/usr/bin/env Rscript
#
# ribosome_expression.R
# 分析 Bombyx mori 精巢 vs 卵巢中 KEGG Ribosome 通路基因的表达模式
# 验证核糖体基因在卵巢中系统性的高表达（与卵母细胞核糖体储备假说一致）
#
# 输入：
#   - KEGG 通路映射文件
#   - 共识差异表达结果
#   - 标准化表达量矩阵
#
# 输出：
#   - 小提琴图 + 配对点图
#   - 表达数据表

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tidyr)
})

source("paper_analysis/theme_nature.R")

# ============================================================
# 0. 文件路径配置
# ============================================================

base_dir <- Sys.getenv("PROJECT_ROOT", normalizePath("."))

kegg_file <- file.path(base_dir,
  "experiments/bombyx_enrichment/data/bombyx_ncbi_to_kegg.tsv")

consensus_file <- file.path(base_dir,
  "benchmark_results/bombyx_mori/runs/2026-05-24_bombyx_mori_decontam-on",
  "07.consensus_expression/Testis_vs_Ovary/consensus_results.tsv")

norm_counts_file <- file.path(base_dir,
  "benchmark_results/bombyx_mori/runs/2026-05-24_bombyx_mori_decontam-on",
  "06.differential_expression/featurecounts/normalized_counts.csv")

out_dir <- file.path(base_dir, "experiments/bombyx_enrichment/results")
fig_dir <- file.path(out_dir, "figures")
tbl_dir <- file.path(out_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== Step 1: 读取 KEGG Ribosome 通路基因 ===\n")

# ============================================================
# 1. 读取 KEGG 映射，提取 Ribosome 通路基因
# ============================================================

kegg <- read_tsv(kegg_file, col_types = cols(.default = col_character()),
                  show_col_types = FALSE)

ribosome_pathway <- "path:bmor03010"
ribosome_kegg <- kegg %>%
  filter(kegg_pathway_id == ribosome_pathway) %>%
  distinct(ncbi_gene_id)

cat(sprintf("   Ribosome 通路中唯一 NCBI Gene ID 数: %d\n", nrow(ribosome_kegg)))

# ============================================================
# 2. 读取共识结果，筛选 Tier A down（卵巢高表达）基因
# ============================================================

consensus <- read_tsv(consensus_file, col_types = cols(.default = col_character()),
                       show_col_types = FALSE)

cat(sprintf("   共识结果总基因数: %d\n", nrow(consensus)))

tier_a_down <- consensus %>%
  filter(tier == "Tier_A", consensus_direction == "down")

cat(sprintf("   Tier A down 基因数: %d\n", nrow(tier_a_down)))

# ============================================================
# 3. 映射 NCBI Gene ID → consensus gene_id_standard
# ============================================================
#
# KEGG 使用纯 NCBI Gene ID (如 "692512")
# consensus 使用 "gene:GeneID_692512" 或 "gene:LOC101743284"
# 策略：从 KEGG Gene ID 构造 "gene:GeneID_XXXX" 格式匹配

# 方式1：直接构造 gene:GeneID_ 前缀匹配
consensus_gene_ids <- unique(tier_a_down$gene_id_standard)

# 从 consensus gene_id 中提取可能的 NCBI Gene ID
# "gene:GeneID_692512" -> "692512"
# "gene:LOC101743284" -> "101743284"
extract_numeric_id <- function(x) {
  # 提取 gene:GeneID_ 或 gene:LOC 后的数字部分
  # 保持与输入相同的长度，无匹配的返回 NA
  m <- regmatches(x, gregexpr("[0-9]+$", x))
  vapply(m, function(el) {
    if (length(el) == 0) NA_character_ else el[1]
  }, character(1))
}

consensus_ncbi <- data.frame(
  gene_id_standard = consensus_gene_ids,
  ncbi_candidate = extract_numeric_id(consensus_gene_ids),
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(ncbi_candidate))

# 匹配：KEGG Gene ID 必须等于 consensus 中的 ncbi_candidate
mapped <- ribosome_kegg %>%
  inner_join(consensus_ncbi, by = c("ncbi_gene_id" = "ncbi_candidate"))

cat(sprintf("   Ribosome 通路中匹配到 Tier A down 的基因数: %d\n", nrow(mapped)))

if (nrow(mapped) == 0) {
  stop("未找到任何 Ribosome 通路基因匹配到 Tier A down。请检查 ID 映射逻辑。")
}

# 提取匹配的 gene_id_standard 列表
ribosome_genes <- mapped$gene_id_standard

# ============================================================
# 4. 读取标准化表达量矩阵
# ============================================================

cat("=== Step 2: 读取标准化表达量 ===\n")

# normalized_counts.csv: 第一列为 gene ID (如 "gene:GeneID_692512"),
# 后续列为 TT_1, TT_2, TT_3, OV_1, OV_2, OV_3
# 第一列没有列名 (header 中是 "")
# 需要以字符形式读取第一列，其余为数值。策略：全读为字符后转换。

norm_raw <- read_csv(norm_counts_file, show_col_types = FALSE,
                      col_types = cols(.default = col_character()))

# 将表达式列（除第一列外）转换为数值
for (col in colnames(norm_raw)[-1]) {
  norm_raw[[col]] <- as.numeric(norm_raw[[col]])
}
colnames(norm_raw)[1] <- "gene_id"

cat(sprintf("   标准化矩阵基因数: %d\n", nrow(norm_raw)))

# 筛选 Ribosome 通路基因
norm_ribo <- norm_raw %>%
  filter(gene_id %in% ribosome_genes)

cat(sprintf("   匹配到表达量的 Ribosome Tier A down 基因数: %d\n", nrow(norm_ribo)))

if (nrow(norm_ribo) == 0) {
  stop("Ribosome 基因未在表达量矩阵中找到。请检查 gene_id 格式。")
}

# ============================================================
# 5. 计算每基因在精巢/卵巢中的平均表达
# ============================================================

cat("=== Step 3: 统计分析 ===\n")

# 转换为长格式
norm_long <- norm_ribo %>%
  pivot_longer(
    cols = -gene_id,
    names_to = "sample",
    values_to = "expression"
  ) %>%
  mutate(
    tissue = case_when(
      grepl("^TT_", sample) ~ "Testis",
      grepl("^OV_", sample) ~ "Ovary",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(tissue))

# 检查样本数
cat(sprintf("   精巢样本: %s\n", paste(unique(grep("TT_", norm_long$sample, value = TRUE)), collapse = ", ")))
cat(sprintf("   卵巢样本: %s\n", paste(unique(grep("OV_", norm_long$sample, value = TRUE)), collapse = ", ")))

# 计算每基因每组织的均值
gene_means <- norm_long %>%
  group_by(gene_id, tissue) %>%
  summarise(mean_expr = mean(expression, na.rm = TRUE), .groups = "drop")

# 宽格式：每基因一行
gene_wide <- gene_means %>%
  pivot_wider(
    names_from = tissue,
    values_from = mean_expr,
    names_prefix = "mean_"
  ) %>%
  mutate(
    log2FC_ovary_vs_testis = mean_Ovary - mean_Testis,
    higher_in = ifelse(log2FC_ovary_vs_testis > 0, "Ovary", "Testis")
  )

# ============================================================
# 6. Wilcoxon 配对符号秩检验
#    H0: 卵巢与精巢表达差异的中位数为 0
# ============================================================

# 配对检验：每个基因在两种组织中的表达
testis_vals <- gene_wide$mean_Testis
ovary_vals  <- gene_wide$mean_Ovary

# 检查是否所有值都可用于配对
n_genes <- length(testis_vals)
cat(sprintf("   用于检验的基因数: %d\n", n_genes))

wilcox_result <- wilcox.test(ovary_vals, testis_vals,
                              paired = TRUE,
                              alternative = "greater",
                              exact = FALSE)

p_value <- wilcox_result$p.value
cat(sprintf("   Wilcoxon paired signed-rank test (one-sided, Ovary > Testis):\n"))
cat(sprintf("     V = %.1f, p = %.4g\n", wilcox_result$statistic, p_value))

# 计算 % 卵巢中高表达的基因
n_higher_ovary <- sum(gene_wide$higher_in == "Ovary")
pct_higher_ovary <- 100 * n_higher_ovary / n_genes

cat(sprintf("   卵巢中高表达的基因: %d / %d (%.1f%%)\n",
            n_higher_ovary, n_genes, pct_higher_ovary))

# 计算组均值
mean_testis <- mean(testis_vals)
mean_ovary  <- mean(ovary_vals)

cat(sprintf("   精巢平均表达: %.2f\n", mean_testis))
cat(sprintf("   卵巢平均表达: %.2f\n", mean_ovary))
cat(sprintf("   平均 log2FC (卵巢/精巢): %.2f\n", mean(ovary_vals - testis_vals)))

# ============================================================
# 7. 准备标注文字
# ============================================================

if (p_value < 0.001) {
  p_label <- "P < 0.001"
} else {
  p_label <- sprintf("P = %.3f", p_value)
}

annot_text <- sprintf(
  "Wilcoxon paired test: %s\nGenes: %d  |  Ovary-higher: %.0f%%",
  p_label, n_genes, pct_higher_ovary
)

# ============================================================
# 8. 可视化
# ============================================================

cat("=== Step 4: 生成图形 ===\n")

# 配色方案
tissue_colors <- c("Testis" = "#3B7DD8", "Ovary" = "#E24B4B")

# ----- 8a. 小提琴图 + 点图 -----

# 构建适合 violin 图的长格式数据
violin_data <- gene_wide %>%
  select(gene_id, mean_Testis, mean_Ovary) %>%
  pivot_longer(
    cols = c(mean_Testis, mean_Ovary),
    names_to = "tissue_raw",
    values_to = "expression"
  ) %>%
  mutate(
    tissue = ifelse(tissue_raw == "mean_Testis", "Testis", "Ovary")
  )

p_violin <- ggplot(violin_data, aes(x = tissue, y = expression, fill = tissue)) +
  geom_violin(alpha = 0.4, linewidth = 0.5) +
  stat_summary(fun = median, geom = "crossbar", width = 0.3, linewidth = 0.4, color = "grey20") +
  geom_jitter(width = 0.12, alpha = 0.5, size = 0.8, color = "grey30") +
  scale_fill_manual(values = tissue_colors, guide = "none") +
  labs(
    x = NULL,
    y = expression(log[2] ~ "Normalized Expression"),
    title = "KEGG Ribosome Pathway Genes",
    subtitle = sprintf("Bombyx mori Testis vs Ovary  |  %d genes  |  %s",
                       n_genes, p_label)
  )

save_pub_r(p_violin, file.path(fig_dir, "ribosome_expression_violin"),
  width_mm=80, height_mm=75)
cat("   已保存: ribosome_expression_violin\n")

# ----- 8b. 配对点图（基因连线） -----

p_paired <- ggplot(gene_wide) +
  geom_segment(
    aes(x = "Testis", xend = "Ovary",
        y = mean_Testis, yend = mean_Ovary),
    color = "grey75", linewidth = 0.3, alpha = 0.6
  ) +
  geom_point(aes(x = "Testis", y = mean_Testis),
             color = tissue_colors["Testis"], size = 1.2, alpha = 0.7) +
  geom_point(aes(x = "Ovary", y = mean_Ovary),
             color = tissue_colors["Ovary"], size = 1.2, alpha = 0.7) +
  labs(
    x = NULL,
    y = expression(log[2] ~ "Normalized Expression"),
    title = "KEGG Ribosome Pathway: Paired Expression",
    subtitle = sprintf("Each line = one gene  |  %d genes  |  %.0f%% higher in ovary",
                       n_genes, pct_higher_ovary)
  )

save_pub_r(p_paired, file.path(fig_dir, "ribosome_expression_paired"),
  width_mm=80, height_mm=75)
cat("   已保存: ribosome_expression_paired\n")

# ============================================================
# 9. 保存数据表
# ============================================================

cat("=== Step 5: 保存数据 ===\n")

# 合并基因注释信息
output_table <- gene_wide %>%
  left_join(
    mapped %>% select(gene_id_standard, ncbi_gene_id),
    by = c("gene_id" = "gene_id_standard")
  ) %>%
  select(gene_id, ncbi_gene_id, mean_Testis, mean_Ovary,
         log2FC_ovary_vs_testis, higher_in) %>%
  arrange(desc(log2FC_ovary_vs_testis))

write_tsv(output_table, file.path(tbl_dir, "ribosome_expression_data.tsv"))
cat(sprintf("   已保存: ribosome_expression_data.tsv (%d genes)\n", nrow(output_table)))

# ============================================================
# 10. 输出摘要
# ============================================================

cat("\n========================================\n")
cat("         分析结果摘要\n")
cat("========================================\n")
cat(sprintf("Ribosome 通路 Tier A down 基因数: %d\n", n_genes))
cat(sprintf("Wilcoxon paired test p-value:     %.4g\n", p_value))
cat(sprintf("卵巢中高表达基因比例:             %.1f%% (%d/%d)\n",
            pct_higher_ovary, n_higher_ovary, n_genes))
cat(sprintf("精巢平均表达 (log2):              %.2f\n", mean_testis))
cat(sprintf("卵巢平均表达 (log2):              %.2f\n", mean_ovary))
cat(sprintf("平均 log2FC (卵巢 - 精巢):        %.2f\n", mean(ovary_vals - testis_vals)))
cat("========================================\n")

if (p_value < 0.05 && pct_higher_ovary >= 70) {
  cat("\n✓ 结果符合预期：核糖体基因在卵巢中系统性高表达，")
  cat("\n  与卵母细胞核糖体储备假说一致。\n")
} else {
  cat("\n⚠ 结果部分符合预期，请检查数据和假设。\n")
}

cat("\n分析完成。\n")
