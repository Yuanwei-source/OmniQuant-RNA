# OmniQuant-RNA 方法与分析边界

## 文档目的

这份文档只记录当前仍然有效、会影响代码修改和结果解释的方法设计。

不在这里保留以下内容：

- 一次性排错记录
- AI 过程草稿
- 阶段性审计全文
- 重复的 benchmark 叙事版本

## 主目标

当前仓库的主目标是宿主转录组定量、差异表达分析与多定量器共识整合。

因此，所有默认设计都优先服务以下主链路：

1. 更干净的宿主 reads 输入
2. 更稳定的宿主比对与定量
3. 更可信的 DEA 结果
4. 更稳健的共识差异表达输出

## Decontam 的边界

Decontam 的首要职责不是“尽量保留所有生态信号”，而是“保护宿主表达分析不被非宿主 reads 干扰”。

当前默认原则如下：

1. `clean` FASTQ 默认只服务宿主主链路。
2. `Uncertain`、非宿主、污染相关 reads 不默认回灌 clean。
3. 这部分 reads 作为审计证据和独立微生物线索侧支输入保留。

对应地，宿主分析主链路为：

1. fastp
2. decontam
3. host alignment
4. quantification
5. single-method DEA
6. consensus aggregation

## 共识层统计原则

当前共识层以多定量器 DESeq2 结果为输入，服务目标是获得更稳健的 gene-level DEG 结论，而不是把多个 quantifier 视为相互独立的重复实验。

当前正式保留的统计原则如下：

### 双重共识引擎（Dual-Consensus Framework）

采用 RRA 与 CCT 作为**并列的双重共识引擎**。任何基因必须同时通过两者的阈值检验才能进入 Tier A/B/C：

1. **Cauchy Combination Test (CCT) — 参数层共识引擎（co-primary）。** 对同方向原始 P 值进行等权 Cauchy 组合。其核心数学优势在于对任意相关结构的稳健性（robust to arbitrary dependency structures）——当 Salmon 与 Kallisto 因算法相似性产生非生物学的排序相关性时，CCT 不会像 RRA 那样面临独立性假设被违反的潜在膨胀风险。因此 CCT 作为并列主引擎而非辅助证据。

2. **Robust Rank Aggregation (RRA) — 非参数层共识引擎（co-primary）。** 评估基因在多定量器排序中的稳定靠前程度。未进入单向 rank list 的基因被保守赋予底层秩次（rank = N），系统性地惩罚缺乏跨定量器方向一致性的基因。该惩罚机制使零假设下 P 值分布在高端向 1.0 偏移，有效压制假阳性。

3. **双重守门（Dual Gatekeeping）。** Tier A 要求 `RRA FDR < 0.05` 且 `CCT FDR < 0.05`。任何因 RRA 相关性膨胀产生的假阳性都会被 CCT 拦截。Tier A 基因经受住了参数和非参数双重考验。

### 方向冲突与胜者诅咒的消除

- 若上调与下调支持数相等 → 直接判 `mixed`，不使用 P 值断案
- 若双方均有信号但一方未达 2 倍优势（如 2 up vs 1 down）→ 判 `mixed`
- `mixed` 与 `none` 方向基因的共识 P 值设为 `NA`，彻底排除出显著层级
- 仅单一工具支持的基因（support_n < 2）自动归入 `unclassified`

### 置换检验验证（100 次，8595 共享基因，4 定量器）

- **关键区间的完美校准**：α = 0.05 处 RRA 经验假阳性率为 4.98%（期望 5.00%）
- **高端的保守性偏移**：零假设下 P 值分布在 >0.5 区间显著向 1.0 偏移，系统惩罚缺乏方向一致性的基因
- **真实信号的强富集**：真实数据在 α = 0.05 处富集约 23%（4.7×），信号由真实生物学差异驱动

该置换检验使用基因标签独立打乱（破坏跨定量器相关性），主要验证 RRA 算法实现正确性。对相关性的稳健性由 CCT 引擎提供。

### Tier 阈值配置（config.yaml）

| Tier | min_support | min_sign_consistency | max_rra_fdr | max_cct_fdr | max_logFC_CV |
|------|-------------|----------------------|-------------|-------------|--------------|
| A | 4 | 4* | 0.05 | 0.05 | 1.00 |
| B | 3 | 3 | 0.10 | 0.10 | 1.25 |
| C | 2 | 2 | 0.25 | 0.25 | 1.50 |

\* Tier A 额外要求 sign_consistency_n == support_n（所有显著工具方向完全一致）

### 不使用动态 CCT 权重

避免基于相关性阈值的动态权重方案，防止方法描述与实现不一致，也避免缺乏稳定方法学依据的启发式惩罚。CCT 的等权设计已有完备数学保证。

### RRA 方法敏感性

为验证 Tier 分级结果对 rank aggregation 方法选择的稳健性，我们并行计算了基于 Borda 均值聚合（`method = "mean"`）的替代 RRA 得分，并以相同的双引擎守门条件重新判定 Tier。输出表 `tier_mean` 列和 `best_mean_rra_fdr` 列记录了基于均值聚合的分级结果。`consensus_summary.tsv` 中的 `rra_mean_tier_a_concordance` 等指标报告了两种聚合方法下维持相同 Tier 分类的基因数。高一致性表明结果对聚合方法选择不敏感。

可直接用于论文或报告的方法学表述如下：

> To rigorously synthesize the results from four correlated quantifiers, we implemented a dual-consensus framework avoiding traditional heuristic weights.
>
> First, to address the issue of non-independence among quantifiers—which can inflate standard rank aggregation methods—we employed the **Cauchy Combination Test (CCT)** as a co-primary parametric engine. CCT is mathematically proven to be robust to arbitrary dependency structures under equal weights, providing reliable P-value combination even when quantifier outputs are correlated due to algorithmic similarities.
>
> Second, as an orthogonal non-parametric validation, we utilized **Robust Rank Aggregation (RRA)** to capture rank consistency across quantifiers. Genes not ranked in a given directional list were conservatively assigned the lowest rank, systematically penalizing those lacking cross-quantifier concordance.
>
> Crucially, to eliminate the "winner's curse," any gene exhibiting conflicting directional signals across quantifiers was strictly labeled as *mixed*—its consensus P-value was set to NA, effectively excluding it from all significance tiers. The consensus direction was determined by a strict majority rule: when support counts were tied or when neither direction held a ≥2× advantage, the direction was deemed *mixed*. Genes supported by fewer than two quantifiers were classified as *unclassified*.
>
> Only genes passing **both** the parametric (CCT FDR < 0.05) and non-parametric (RRA FDR < 0.05) thresholds, with consistent directional support and adequate effect size stability, were assigned to high-confidence Tier A. This dual-gatekeeping strategy ensures that false positives arising from correlation-induced RRA inflation are intercepted by the correlation-robust CCT.

## 微生物线索侧支的定位

微生物线索分析是独立侧支，不是主结果支柱。

它的目标是：

1. 解释样本异质性
2. 标记高污染或异常样本
3. 做重点目标的存在性筛查
4. 为后续课题提供入口

它不承担以下任务：

- 替代正式宏转录组研究
- 替代微生物表达定量研究
- 作为宿主 DEA 主输入

## 侧支分析的最小范围

为了保证有意义但不过度扩张，微生物线索侧支只做三件事：

1. 样本级微生物负担概览
2. 重点目标存在性筛查
3. 与宿主结果做轻量关联

对应只需要四类输出：

1. 样本级汇总表：每个样本的非宿主比例、细菌比例、真菌比例、病毒比例、Top taxa
2. 重点目标面板：Wolbachia、病毒、真菌、环境细菌等存在性与相对强弱
3. 微生物组成图：堆叠条形图或热图
4. 宿主结果解释图：在宿主 PCA 或样本聚类图上叠加微生物负担标记

这些侧支结果统一挂载到 `results/03.decontam/clues/`，与证据源保持邻近，但不回灌宿主主链路。

侧支中的重点目标 taxid 和 Top taxa 保留数量通过 `config.yaml` 中的 `decontam.clues.priority_targets` 与 `decontam.clues.top_taxa_n` 显式控制，而不是依赖脚本硬编码。

## 文档保留原则

远程仓库只保留两类文档：

1. 总入口文档
2. 当前有效的方法与决策文档

以下内容默认不推远程仓库：

1. `experiments/audits/` 下的一次性审计和排错记录
2. `experiments/drafts/` 下的草稿文案
3. 未定稿的 benchmark 方案
4. AI 生成但未被采纳的过程性总结

如果某份本地文档里有长期有效的结论，应将结论提炼并合并进本文件或 [README.md](README.md)，而不是继续堆积版本。

## 变更原则

以后只在以下情况下新增正式文档：

1. 新增了一个会长期保留的分析模块
2. 现有 README 无法承载核心方法边界
3. 某个决策会持续影响后续代码实现与结果解释

否则优先更新现有正式文档，不新增同主题重复文件。
