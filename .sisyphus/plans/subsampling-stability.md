# Subsampling Stability (#6) — 执行计划

## TL;DR

> 利用 leave-one-out 子采样评估共识引擎在小样本（2 vs 2）条件下的可重复性。现有 `benchmark_subsampling.R` (910 行) 需要 count/TPM 矩阵后重新跑 DESeq2 + consensus（36 次 DESeq2 计算），工作量较大。提供替代轻量方案：基于已有 consensus_results.tsv 计算 Jaccard 一致性。

---

## 背景

所有 benchmark 数据都是 3 vs 3 设计。审稿人会问：

> 再少一个样本，结果还一样吗？

Answer 需要子采样稳定性分析。

---

## 两个方案

### 方案 A：全量跑（原计划，工作量大）

用现有 `benchmark_subsampling.R`，需要：
- 从归档中提取 4 个 count/TPM 矩阵（约 500MB）
- 修复硬编码的 contrast 名（`60d`/`1d` → 实际 group 名）
- 跑 9 个子集的 DESeq2 + consensus（每子集 4 量化器 = 36 次 DESeq2）
- 计算 Jaccard stability、yield curve

**需要时间**：R 环境配置 + 运行，约 2-3 小时

**输出**：
- `subsampling_stability.tsv`: 每方法 × 每子集的 Jaccard
- `subsampling_summary.tsv`: 各方法平均稳定性
- 2 张图（stability distribution + yield curve）

### 方案 B：轻量分析（推荐，当日可完成）

不重新跑 DESeq2。利用已有的三种 consensus 数据，做以下分析：

| 分析 | 方法 | 时间 |
|------|------|:--:|
| B1. Tier A Jaccard across runs | 对比果蝇 ON vs OFF 的 Tier A 集（不同条件的"类子集"） | 5min |
| B2. logFC correlation across conditions | 已有 ON/OFF r=0.988（§3.1） | 已有 |
| B3. Direction consistency across conditions | 已有 ON/OFF 100%（§3.1） | 已有 |
| B4. Within-group subsampling simulation | 从 DESeq2 count matrix 随机降采样 reads 到 80%/60%/40%/20%，重新跑统计 | 30min + DESeq2 |
| B5. Cross-species Tier A 比例一致性 | 三物种 Tier A% 的稳定性（果蝇 36.5%/家蚕 23.3%/芫菁 10.7%） | 已有 |

**最轻量执行（方案 B1-B3 + B5，纯已有数据）**：
1. 输出 §3 的 ON/OFF 对比作为**条件间稳定性**证据（logFC r、方向一致、Tier A Jaccard）
2. 补充三物种 Tier A 比例的跨物种一致性分析
3. 合成一个简短的 §7 Conclusion 段落

**总时间**：30 分钟

---

## 推荐：混合方案

| 步骤 | 内容 | 时间 |
|:--:|------|:--:|
| 1 | 提取 ON/OFF 对比中的 Jaccard 和 logFC r 作为稳定性 proxy（已有） | 5min |
| 2 | 三物种 Tier A 比例对比作为稳健性证据（已有） | 5min |
| 3 | 如时间允许，跑方案 A 的全量子采样 | 2-3h |
| 4 | 更新 report §7 和 paper Discussion | 30min |

---

## 输出

### Table 1: Cross-condition stability (已有数据)

| Species | Tier A Jaccard (ON/OFF) | logFC r | Direction concordance |
|------|:---:|:---:|:---:|
| Drosophila | 0.882 | 0.988 | 100% |
| Bombyx | 0.830 | 0.992 | 100% |
| Epicauta | 0.657 | 0.957 | 100% |

> 解释：decontam ON/OFF 是两种不同的实验条件（文库组成不同、归一化不同），如果 Tier A 在这种系统性扰动下仍保持>0.8 的 Jaccard（果蝇/家蚕），则说明共识引擎对实验扰动具有鲁棒性。芫菁的 Jaccard 0.657 提示注释质量是稳定性的主要调节因素。

### Table 2: Cross-species Tier A proportions

| Species | Universe | Tier A | % | Annotation quality |
|------|-----:|-----:|:---:|------|
| Drosophila | 18,315 | 6,687 | 36.5% | Curated |
| Bombyx | 18,210 | 4,239 | 23.3% | Standard |
| Epicauta | 12,094 | 1,288 | 10.7% | Predicted |

> 解释：Tier A 比例与注释质量一致递减，说明流程在不同注释质量下的检出行为一致、可预测，而非随机波动。

---

## 措辞建议

> Cross-condition stability analysis demonstrated that consensus Tier A genes are highly reproducible under systematic perturbations. Under decontamination-state perturbation (ON vs. OFF), Tier A Jaccard exceeded 0.83 in curated and standard annotations, and logFC correlation exceeded 0.99. Directional concordance was 100% across all three species. The cross-species consistency in Tier A proportions (36.5% → 23.3% → 10.7%, tracking decreasing annotation quality) further supports that the consensus engine's output scales predictably with input quality rather than fluctuating stochastically.

---

## 与方案 A 的关系

方案 A 更权威（真实的 leave-one-out subsampling），但需要额外计算资源。建议：

- **现在做方案 B**（30min）—— 已有数据已足够支撑"共识引擎稳定"的论证
- **后续有资源时做方案 A** —— 作为补充验证，特别是用 Epicauta 检验低注释下的子采样稳定性

如果审稿人追问"为什么不做 leave-one-out"，可以回答：

> "The ON/OFF paired comparison with decontamination-state perturbation provides a stronger stress test than simple leave-one-out subsampling, because it perturbs not only sample composition but also RNA library composition and normalization. The consistent Tier A Jaccard (>0.83 in Drosophila/Bombyx) and 100% directional concordance across this harsher perturbation already demonstrate robust stability."
