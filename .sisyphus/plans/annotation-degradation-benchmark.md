# Annotation Degradation Benchmark (#1)

## TL;DR

> 对果蝇 ON run 的 DESeq2 结果做 GFF 退化模拟，重新跑共识引擎，计算 precision/recall/F1，输出到 `benchmark_results/drosophila_wolbachia/analysis/annotation_degradation/`

## 前置知识

### GFF 退化原理

- **不重跑 Snakemake**，只拿已有的 4 个 DESeq2 CSV，在 R 内模拟"注释不完整"
- 退化后的基因：`P.Value ← 1`、`adj.P.Val ← 1`、`logFC ← NA`（模拟"注释里没有这个基因所以测不到"）
- 然后重新跑 consensus engine，用完整注释的 Tier_A 集做 gold standard 算 precision/recall/F1

### 5 种退化模式

| 模式 | 操作 | 现实含义 |
|------|------|---------|
| Random drop | 随机删除 X% gene | 基因组碎片化 |
| Length-biased | 最短 X% gene 优先删 | 小基因注释缺失 |
| Expression-biased | 最低表达 X% gene 优先删 | 非模式生物缺组织特异基因 |
| Transcript-level | 删转录本但保留 gene | isoform 不全（最常见） |
| ID corruption | 改 X% gene_id | 跨来源 ID 不兼容 |

### 执行参数

- 每种模式 × 4 级强度 (0%, 25%, 50%, 75%) × 5 seeds = **100 次模拟**
- 纯 R 计算，几分钟出结果

---

## TODOs

- [x] Step 1: 修复 `workflow/rules/benchmark.smk`（5 条过时路径已更新）
- [x] Step 2: 恢复果蝇 ON 结果到 `results/`（symlink 06.differential_expression + 05.quantification）
- [x] Step 3: 切换 config 到果蝇 ON
- [x] Step 4: 干跑检查 `snakemake --dry-run`
- [x] Step 5: 运行 benchmark（`--use-conda --cores 4`）
- [x] Step 6: 归档结果到 `benchmark_results/drosophila_wolbachia/analysis/annotation_degradation/`
- [x] Step 7: 清理临时文件
- [x] Step 8: 更新论文（annotation degradation §3.2 + Table 2 + Abstract）
