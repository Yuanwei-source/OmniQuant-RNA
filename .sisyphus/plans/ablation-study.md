# Ablation Study (#4) — 模块贡献量化（v2.1）

## TL;DR

> 利用已有数据，量化 OmniQuant-RNA 各模块的独立贡献。**核心问题不是"每个 gate 过滤了多少基因"，而是"每个 gate 过滤掉的是什么类型的基因"**。通过三物种门控失败重叠分析、支持度-only baseline、参数阈值扫描、量化器重叠矩阵和 namespace 冲突检查，最终回答"共识引擎的收益到底来自哪里"。

---

## 核心命题

> Ablation analysis showed that the consensus gates do not merely reduce DEG counts; they preferentially remove low-expression, borderline-significant, method-sensitive candidates while preserving the high-confidence shared DEG core across species.

---

## 审稿人会问的 9 个问题

| # | 问题 | 方法 | 工作量 |
|---|------|------|:--:|
| Q1 | 纯投票（不用 RRA/CCT/CV）能叫出多少 DEG？ | 三物种 support_n+dir 条件过滤 | **新增** |
| Q2 | RRA/CCT/CV 各自的独立过滤数 vs 冗余过滤数？ | Gate failure overlap 矩阵（三物种） | **新增** |
| Q3 | 被各 gate 过滤掉的基因有什么特征？ | 表达量、padj、注释率、FC 支持率（三物种） | **新增** |
| Q4 | 单量化器之间有多大不一致？ | 4 量化器 DEG 重叠矩阵 + 支持度分布 | **新增** |
| Q5 | Namespace 模块是否真的无信息丢失？ | 碰撞/歧义/孤立 ID 检查 | **新增** |
| Q6 | 参数阈值是否是人调出来的？ | 关键参数轻量 sweep（基于已有数据） | **新增** |
| Q7 | 去污染的独立贡献？ | 已有 ON/OFF 对比 | 已有 |
| Q8 | 单定量器 vs 共识的 precision/recall 差距？ | 已有注释退化 benchmark | 已有 |
| Q9 | 结论是否跨物种成立？ | **所有分析 × 三物种**（果蝇/家蚕/芫菁） | — |

---

## 5 项新增分析

### Analysis A: 三物种 Vote-only baseline（Q1 + Q9）

**定义**（避免与正式 Tier 混淆）：

| 配置名 | 条件 |
|------|------|
| **Support-strict** | support_n ≥ 4 AND sign_consistency_n ≥ 4 |
| **Support-moderate** | support_n ≥ 3 AND sign_consistency_n ≥ 3 |
| **Support-relaxed** | support_n ≥ 2 AND sign_consistency_n ≥ 2 |

**输出**（三物种各一张表）：

| 配置 | DEG数 | Full Tier A overlap | Δ vs Full | median expr | median -log10(padj) | annotation rate |
|------|-----:|-----:|-----:|-----:|-----:|-----:|
| Full Tier A | 6,687 | 100% | — | ... | ... | ... |
| Support-strict | X | Y% | +Z | ... | ... | ... |
| Support-moderate | X | Y% | +Z | ... | ... | ... |
| Support-relaxed | X | Y% | +Z | ... | ... | ... |

关键解读：被 statistical gates 额外过滤的基因 = Support-strict − Full Tier A。这些基因的"质量 proxy"回答了 Q3。

### Analysis B: Gate failure overlap 矩阵（Q2 + Q3 + Q9）

**定义两步候选全集**：

| Universe | 定义 | 用途 |
|------|------|------|
| **Any-signal** | 至少 1 个量化器 DESeq2 padj<0.05 | 看所有潜在 DEG 候选 |
| **Vote-pass** | support_n ≥ 3 AND sign_consistency_n ≥ 3（通过投票但未进入 Full Tier A）| 主表：量化 RRA/CCT/CV 的独立门控贡献 |

**主分析用 Vote-pass universe**，因为它最能回答："在已经有多量化器支持和方向一致的候选中，统计门控进一步排除了什么？"

**输出**（三物种，Vote-pass universe）：

| 层级 | 失败组合 | 基因数 | % | median expr | -log10(padj) | % named | % FC sig |
|------|------|-----:|-----:|-----:|-----:|-----:|-----:|
| **统计门控** | fail RRA only | X | X% | ... | ... | ... | ... |
| | fail CCT only | X | X% | X | X | X | X |
| | fail logFC_CV only | X | X% | X | X | X | X |
| | fail RRA + CCT | X | X% | ... | ... | ... | ... |
| | fail RRA + CV | X | X% | ... | ... | ... | ... |
| | fail all stat gates | X | X% | | | | |
| **方向/支持** | fail support_n only | X | X% | ... | ... | ... | ... |
| | fail direction_n only | X | X% | ... | ... | ... | ... |
| | mixed_direction | X | X% | | | | |
| **背景** | no_signal | X | X% | 极低 | 极低 | — | — |

**关键解读**：
- `fail CCT only` 数量 → CCT 独立过滤贡献
- `fail RRA + CCT` 数量 → 两个统计门控的冗余确认层规模
- 每个失败组合的 expression/padj → 质量 proxy

**"Added-back" 严格定义**（避免 Table 2 中的重复计数质疑）：

| 术语 | 定义 |
|------|------|
| `Added by −RRA` | 移除 RRA gate 后进入 Tier A-like（即满足除 RRA 外的所有 Tier A 条件），但不在 Full Tier A 中 |
| `Added by −CCT` | 移除 CCT gate 后进入 Tier A-like，但不在 Full Tier A 中 |
| `Added by support-strict` | 满足 support_n≥4 + direction_n≥4（不要求 RRA/CCT/CV），但不在 Full Tier A 中 |
| `fail RRA only` | 满足所有其他 4 个 Tier A 条件（support/direction/CCT/CV），仅 RRA FDR 不达标 |
| `fail CCT only` | 满足所有其他 4 个 Tier A 条件，仅 CCT FDR 不达标 |
| `fail RRA+CCT` | 满足 support/direction/CV 三个条件，但 RRA 和 CCT FDR 均不达标 |

**可选 negative control**（列入 plan 但不强制本次执行）：

> 随机打乱样本标签后，Full Tier A 是否几乎不产生 DEG？Support-only 是否更容易产生 DEG？
> 需要重新跑差异分析，工作量大。本次 lightweight ablation 不包含此项，留待后续 stability benchmark 处理。

**输出措辞模板**：

> CCT provides minimal unique filtering beyond RRA (X genes fail CCT-only across three species). However, the RRA+CCT overlap set (X genes) represents a redundant confirmation layer: these candidates passed support+direction voting but failed both statistical tests, suggesting systematically weaker signal rather than stochastic filtering.

### Analysis C: 参数阈值 sweep（Q6 + Q9）

基于已有 consensus_results.tsv 中的原始值（RRA FDR、CCT FDR、logFC_CV），不需要重新计算——只需对已存在的列施加不同阈值。

| 参数 | 扫描值 |
|------|------|
| support_n | 2, 3, 4 |
| RRA FDR max | 0.01, 0.05, 0.10, 0.25 |
| CCT FDR max | 0.01, 0.05, 0.10, 0.25 |
| logFC_CV max | 0.5, 0.75, 1.0, 1.25, 1.5 |
| sign_consistency_n min | 2, 3, 4 |

固定其他参数在默认值，每次只变化一个参数，统计 DEG 数 + Full Tier A overlap。

**输出**：参数-DEG 曲线表 + 两张轻量二维网格。

**二维 sweep 1**: RRA FDR max × CCT FDR max（固定 support=4, dir=4, CV=1.0）

| RRA \ CCT | 0.01 | 0.05 | 0.10 | 0.25 |
|-----------|-----:|-----:|-----:|-----:|
| 0.01 | X | X | X | X |
| 0.05 | X | X | X | X |
| 0.10 | X | X | X | X |
| 0.25 | X | X | X | X |

**二维 sweep 2**: support_n min × sign_consistency_n min（固定 RRA=0.05, CCT=0.05, CV=1.0）

| sup \ dir | 2 | 3 | 4 |
|-----------|---|---|---:|
| 2 | X | X | X |
| 3 | X | X | **当前** |
| 4 | X | X | X |

这证明当前阈值组合不是两个参数互相制约卡出来的孤立点。

**关键措辞**：

> Threshold sweeps showed smooth and bounded effects on DEG counts across a 2-5× range of each parameter, with no abrupt cliffs. The current thresholds are positioned at conservative but not extreme points on these curves.

### Analysis D: 量化器重叠矩阵 + 支持度分布（Q4 + Q9）

输出两张表：

**表 D1: 两两重叠（Jaccard）** — 对三物种，4 个量化器 DESeq2 显著 DEG (padj<0.05, |logFC|>1)：

| | FC | ST | SA | KA |
|--|:--:|:--:|:--:|:--:|
| FC | X | J | J | J |
| ST | | X | J | J |
| SA | | | X | J |
| KA | | | | X |

**表 D2: 支持度分布** — 在 **union DEG candidate universe**（至少 1 个量化器 DESeq2 padj<0.05 & |logFC|>1）中统计：

| 支持量化器数 | 基因数 | 进入 Full Tier A 比例 | median expr | median padj | annotation rate |
|:---:|-----:|-----:|-----:|-----:|-----:|
| 1/4 | X | X% | ... | ... | ... |
| 2/4 | X | X% | | | |
| 3/4 | X | X% | | | |
| 4/4 | X | X% | | | |

关键解读：低支持度基因是否更低表达、更边界、更少进入 Tier A。如果 1/4 支持的 DEG 表达更低、显著性更弱且很少进入 Full Tier A，则强力支持共识策略的合理性。

### Analysis E: Namespace 深度检查（Q5）

不仅检查 gene_id 是否在 tx2gene 中有映射，还检查：

| 指标 | 方法 | 重要性 |
|------|------|------|
| Unmapped gene rate | tx2gene vs 各量化器 gene_id 集合 | 信号丢失 |
| Duplicated gene_id rate | 各量化器结果中重复 gene_id 比例 | 重复计数 |
| Many-to-one collapse | tx2gene 中 transcript → gene 的平均映射数 | 过度合并 |
| Quantifier-specific orphans | 每个量化器独有的、不在其他量化器中出现的 gene_id 数 | 量化器特异信号 |
| DEG affected by namespace | ON-only/Shared/OFF-only 中 unmapped 比例差异 | ID 统一是否偏倚 |

**输出**：一张 namespace 完整性表。

| ID normalization effect | 比较 stripping version suffix (.1/.2)、raw ID、canonical ID 三种模式下的 recovered gene 数 | ID 规范化不依赖激进策略 |
| Version suffix robustness | 统计带/不带 version suffix 的 gene_id 比例及两者间的 mapping 一致性 | Ensembl/GFF 版本号是否被正确处理 |

**关键措辞**：

> Namespace recovery was robust to version suffix handling and did not depend on aggressive ID normalization. No measurable gene loss, no quantifier-specific orphan DEG enrichment, and no differential mapping bias among Shared, ON-only, and OFF-only DEG categories were detected.

---

## 最终产出：两张表格

### Table 1: 模块贡献分解

| 配置 | Species | DEG数 | ΔDEG | Full TA overlap | 说明 |
|------|------|-----:|-----:|-----:|------|
| Full (all gates) | DM/BM/EP | ... | — | 100% | 基准 |
| Support-strict | all 3 | ... | +S₂ | O₂% | 纯投票，无统计门控 |
| −RRA gate | all 3 | ... | +S₃ | O₃% | RRA 独立与冗余过滤贡献 |
| −CCT gate | all 3 | TBD | TBD | TBD | 量化 CCT 独立与冗余贡献 |
| −logFC_CV gate | all 3 | TBD | TBD | TBD | 量化效应量稳定性门控贡献 |
| −Decontamination | all 3 | ... | +S₄ | — | ON/OFF 不对称 |
| −Namespace | all 3 | TBD | TBD | TBD | 信息丢失 + 冲突检查 |
| Single FC only | DM | ... | — | — | 见 annotation degradation |

### Table 2: Ablated-in 基因的质量特征

| Ablated-in gene set | vs Full TA | median expr | effect size (Cliff's δ) | p_adj | median |logFC| |
|------|------|-----:|-----:|-----:|-----:|

每列的质量 proxy 配检验：
- expression / padj / |logFC| → **Wilcoxon rank-sum test**（两组）或 **Kruskal-Wallis**（多组），FDR 校正
- annotation rate / FC sig rate → **Fisher's exact test**
- effect size 优先用 **Cliff's delta**（对非正态分布稳健）；组间差异汇报 **fold-change** 或 **odds ratio**

---

## 执行步骤（总计 7 步）

| Step | 内容 | 时间 |
|:--:|------|:--:|
| 1 | 建立三物种数据源 symlink | 2 min |
| 2 | Analysis A: 三物种 vote-only baseline | 30 min |
| 3 | Analysis B: Gate failure overlap 矩阵 | 30 min |
| 4 | Analysis C: 参数阈值 sweep | 20 min |
| 5 | Analysis D: 量化器重叠 + 支持度分布 | 20 min |
| 6 | Analysis E: Namespace 深度检查 | 20 min |
| 7 | 合成 Table 1/2 + 更新 report + 更新 paper | 60 min |

**总时间**：约 3 小时。全部不需要重跑流水线。

---

## 脚本产出

```
benchmark_results/scripts/ablation_study/
├── vote_baseline.py           # Analysis A: 三物种 support-only
├── gate_failure_overlap.py    # Analysis B: 门控失败重叠矩阵
├── threshold_sweep.py         # Analysis C: 参数阈值扫描
├── quantifier_overlap.py      # Analysis D: 量化器重叠 + 支持度分布
├── namespace_integrity.py     # Analysis E: Namespace 深度检查
└── synthesize_ablation.py     # 合成 Table 1 + Table 2
```

---

## 关键措辞规则

| 避免 | 使用 |
|------|------|
| "Tier A 等效" | "Support-strict" |
| "false positive" | "low-confidence / method-sensitive / borderline" |
| "CCT 无贡献" | "CCT provides minimal unique filtering but serves as redundant confirmation layer" |
| "Precision 0.991"（在 ablation 中）| "Reference-consensus precision" 或 "Full Tier A overlap %" |
| 只在 Drosophila 做 | 至少三物种都覆盖 |
| 仅 DEG 数对比 | DEG 数 + 质量 proxy（median expr/padj/|lfc|/name rate/FC support）|
