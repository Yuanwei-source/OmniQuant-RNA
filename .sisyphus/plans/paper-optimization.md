# Paper Optimization Plan (v2)

## TL;DR

> 6 个 phase：冻结事实 → 润色英语 → 提取引用声明 → 补引用 → reviewer readiness audit → 最终一致性检查。关键：先锁数字再改语言；引用优先方法原始论文而非 Nature；最后做 claims overreach 审计。

## 执行流程

```
Phase 0: factual lock table (冻结所有数字)
    ↓
Phase 0.5: figure/table architecture (图表策略)
    ↓
Phase 1: nature-polishing (语言润色，围绕 frozen values + 图表)
    ↓
Phase 1.5: citation-claim extraction (提取需要引用的声明)
    ↓
Phase 2: nature-citation (补引用，优先原始方法论文)
    ↓
Phase 3: reviewer-readiness audit (论证完整性 + claims overreach)
    ↓
Phase 4: final consistency polish (数字/引用/结论一致性)
```

---

## Phase 0 — factual lock table（5 min）

在改任何文字前，锁定所有关键数字。Phase 1 只允许围绕这些值改语言，不允许改数字或扩大解释。

```
FACTUAL LOCK TABLE
═══════════════════════════════════════════════════════════════
Drosophila ON/OFF Tier A Jaccard      0.882
Bombyx ON/OFF Tier A Jaccard          0.830
Epicauta ON/OFF Tier A Jaccard        0.657
Drosophila logFC r                    0.988
Bombyx logFC r                        0.992
Epicauta logFC r                      0.957
Direction concordance                 100% across three species
Drosophila LOO Consensus Jaccard      0.856
Best single-quantifier LOO Jaccard    Salmon 0.851
Tier A proportions                    36.5%, 23.3%, 10.7%
Top-100 consensus limitation          Consensus not best (0.573 vs FC 0.607)
Expression-biased consensus F1 (50%)   0.825
Expression-biased FC F1 (50%)          0.746
Decontamination OFF-only/ON-only asymmetry  3.0× / 4.0× / 2.3×
CCT unique filtering                  DM:0 BM:0 EP:12
RRA unique filtering                  DM:887 BM:299 EP:0
Namespace recovery                    100%
```

### Phase 0 产出
- 事实锁定表（供后续所有 phase 对照）

---

## Phase 0.5 — Figure/Table architecture（15 min）

### 目的
在 Phase 1 润色之前确定图表结构，因为图表会反过来影响 Results 怎么写。

### 主文配置：5 图 + 2 表

#### Figure 1: Workflow schematic
概念图，一眼说明方法创新。
- 双轨去污染（Kaiju + Bowtie2）
- 四定量器并行（FC/ST/SA/KA）
- RRA + CCT + logFC_CV + support/direction gates
- 输出 Tier A/B/C

#### Figure 2: Annotation degradation benchmark
证明共识引擎在注释退化下 precision-oriented。
- A: 退化模式示意
- B: Consensus vs FC precision/recall
- C: F1 across degradation modes
- D: Detectable recall

已有数据：`benchmark_results/drosophila_wolbachia/analysis/annotation_degradation/figures/`

#### Figure 3: Decontam perturbation + P1 OFF-only characterization
证明 OFF-only 是低表达、边界显著、方法敏感候选。
- A: ON/OFF shared vs OFF-only vs ON-only 定义
- B: 四类基因表达量分布（violin/boxplot）
- C: -log10(padj) 分布
- D: featureCounts support rate
- E: annotation status（ON-only 未注释基因富集 4.6×）

#### Figure 4: Ablation study
回答"模块到底贡献什么"。
- A: Gate failure overlap（三物种）
- B: Added-back genes quality profile
- C: RRA×CCT threshold heatmap
- D: support×direction sweep

#### Figure 5: Stability analyses
证明不是偶然结果。
- A: ON/OFF Tier A Jaccard（三物种）
- B: ON/OFF logFC correlation
- C: Drosophila leave-one-out Jaccard（Consensus vs single）
- D: Top-100/Top-1000/full DEG set stability

已有数据：`benchmark_results/drosophila_wolbachia/analysis/subsampling_stability/figures/`

#### Table 1: Datasets and benchmark design

| Species | Annotation | Contrast | Replicates | Genome | Benchmark role |
|------|------|:--:|:--:|------|------|
| Drosophila melanogaster | Curated (BDGP6.54) | Wolb_inf vs Free | 3 vs 3 | 18,315 genes | Model / high-quality |
| Bombyx mori | Standard (Ensembl Metazoa) | Testis vs Ovary | 3 vs 3 | 18,210 genes | Intermediate |
| Epicauta impressicornis | Predicted (AUGUSTUS) | Diapause vs Non | 3 vs 3 | 12,094 genes | Non-model / low-quality |

#### Table 2: Key benchmark summary

| Benchmark | Main metric | Key result | Interpretation |
|------|------|------|------|
| Annotation degradation | Reference-consensus precision | ≥0.985 | High precision under degradation |
| Decontam perturbation | ON/OFF Tier A Jaccard | 0.882/0.830/0.657 | Robust, modulated by annotation |
| P1 OFF-only | Expression/padj/FC support | Lower, weaker | Enriched for low-confidence |
| Ablation | Added-back quality | Weaker than Tier A | Gates refine candidate quality |
| Leave-one-out | Jaccard | 0.856 vs Salmon 0.851 | Comparable stability, conservative |

### 补充材料（6-8 图 + 5-6 表）

| Figure | 内容 | 数据状态 |
|------|------|:--:|
| S1 | 所有 degradation 模式 precision/recall | ✅ 已有 |
| S2 | ON/OFF DEG category counts | ✅ 可生成 |
| S3 | P1 四分位数分布 | ✅ 可生成 |
| S4 | 三物种 gate failure overlap | ✅ 已有 |
| S5 | 完整 threshold sweep | ✅ 已有 |
| S6 | Namespace recovery | ✅ 已有 |
| S7 | Top-k stability curves | ✅ 已有 |
| S8 | Quantifier overlap matrix | ✅ 已有 |

### Phase 0.5 产出
- 主图/主表清单 + panel 级信息标注
- 补充图/表清单
- Figure-to-claim 映射表（每张图支撑论文的哪个 claim）
- 缺图清单（哪些需要新做，哪些已有）

---

### 输入
`experiments/drafts/omniquant-paper-v1.md` + Phase 0 事实表

### 润色策略

| Section | 重点 |
|------|------|
| **Title** | curiosity + credibility，不夸大 |
| **Abstract** | 只放 3 类核心证据：annotation degradation + decontam perturbation + ablation/stability。不放所有数字 |
| **Introduction** | gap-driven：三个挑战 → 现有方案不足 → 本方案填补什么 |
| **Methods** | 可复现性描述，术语一致性（RRA/CCT/Tier/gate 全文统一） |
| **Results** | Results tone（was detected）≠ Discussion tone（may reflect） |
| **Discussion** | 强化 hedging，避免 "proves/demonstrates"，保持 reviewer-safe |
| 全文 | 句子 ≤30 words，paragraph 单一思想，中文痕迹消除 |

### 不能碰
- 数字不改（对照 Phase 0 表）
- 结论强度不改（已经过专家审核）
- 不发明数据/引用

### 产出
- 润色后论文 + 3-5 条修订说明

---

## Phase 1.5 — citation-claim extraction（5 min）

Phase 1 后、Phase 2 前，提取所有需要外部引用的声明。避免乱补引用。

### 提取格式

```
C001: "Wolbachia is present in 40-66% of insect species"
  → needs citation (meta-analysis/review)
C002: "RRA integrates ranked lists across independent studies"
  → needs citation (original RRA paper)
C003: "CCT provides valid combined p-values under arbitrary correlation"
  → needs citation (original CCT paper)
C004: "nf-core/rnaseq provides multi-quantifier outputs without integration"
  → needs citation (nf-core paper)
C005: "OmniQuant-RNA achieved reference-consensus precision 0.991"
  → INTERNAL result, no external citation needed
```

### 产出
- 引用声明列表（10-15 条），标记 citation needed vs internal
- 每条声明标注引用优先级（1=方法原始论文, 2=review, 3=高影响期刊补强）

---

## Phase 2 — nature-citation（15-20 min）

### 引用优先级（非 Nature-first）

| 优先级 | 引什么 | 例子 |
|:--:|------|------|
| **1st** | 原始方法/工具论文 | DESeq2, RRA, CCT, Salmon, Kallisto, featureCounts, StringTie, Bowtie2, Kaiju, Snakemake |
| **2nd** | 领域综述/benchmark | non-model insect transcriptomics, Wolbachia prevalence, RNA-seq best practices |
| **3rd** | 竞争 pipeline 论文 | nf-core/rnaseq, other consensus DEG tools |
| **4th** | 高影响期刊补强 | Nature Methods, Genome Biology, Bioinformatics |

### 不搜
- 数据库 accession（DRA008737）
- 软件版本引用（pipeline 依赖文档）
- 内部结果（Tier A Jaccard, OFF-only feature, ablation, namespace recovery, LOO Jaccard 等）

### 产出
- 分段引用对应表（S001 → 候选 + 支撑等级）
- `.enw` 参考文献管理器文件
- 缺失引用缺口清单

---

## Phase 3 — reviewer-readiness audit（10-15 min）

### 论证完整性检查

| 检查项 | 说明 |
|------|------|
| Abstract 覆盖全部核心结果？ | deco+degradation+ablation/stability，不放所有数字 |
| 每个 claim 有对应 benchmark/table？ | Results 无 unreferenced claims |
| Discussion 完整？ | ablation ✅, stability ✅ |
| Limitations 与最新结果一致？ | subsampling 表述匹配 §8.2 |
| 数量跨章节一致？ | Abstract/Results/Discussion 数字对照 Phase 0 |
| 引用密度？ | Introduction 无连续 3+ 句无引用 |

### Claim strength audit

| 风险措辞 | 改法 |
|------|------|
| proves / demonstrates | 是否真有因果证据？若无 → suggest / indicate / is consistent with |
| false positive | → low-confidence / method-sensitive / borderline-significant |
| caused by microbial reads | → consistent with decontamination-sensitive normalization effects |
| superior to | → comparable or slightly higher / more conservative |
| robust across species | → robust under tested perturbations |

### 关键审稿问题预检
1. 方法创新点是否区别于 nf-core/rnaseq、MultiQC、普通 consensus DEG？
2. 是否有 self-referential evaluation 的透明声明？
3. 负面结果（Top-100 consensus 不如单定量器、Epicauta stability 低）是否正面解释？
4. Limitations 是否主动承认 annotation degradation 可能不代表真实非模式场景？

### 产出
- 论证缺口清单（如有）
- Claim strength 审计结果（哪些句子需要弱化）

---

## Phase 4 — final consistency polish（10 min）

Phase 2 补引用后，句子可能变长变乱。最后过一遍：

1. 数字一致性（对照 Phase 0 表逐项检查）
2. 引用是否真正支撑对应句子（每个 [N] 点开检查）
3. Abstract 是否仍不过度承诺
4. Limitations 是否与最新 benchmark 结果一致

### 产出
- 定稿论文

---

## 预计时间

| Phase | 内容 | 时间 |
|:--:|------|:--:|
| 0 | Factual lock table | 5 min |
| 0.5 | Figure/Table architecture | 15 min |
| 1 | nature-polishing | 20-30 min |
| 1.5 | Claim extraction | 5 min |
| 2 | nature-citation | 15-20 min |
| 3 | Reviewer-readiness audit | 10-15 min |
| 4 | Final polish | 10 min |
| **Total** | | **~90 min** |
