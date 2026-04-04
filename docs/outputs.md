# OmniQuant-RNA 正式输出说明

## 文档目的

这份文档只解释正式结果目录的语义，以及每类输出该如何被后续分析使用。

它不记录一次性 benchmark、临时排错截图或过程性草稿。

## 输出分层总览

主结果目录分为四层语义：

1. 参考与命名空间输出
2. 质量控制与去污染输出
3. 宿主定量与差异表达输出
4. 共识整合与报告输出

## 1. 参考与命名空间

目录：`results/00.reference`

用于统一各定量器后续导入与聚合的参考层输出，包括：

- `tx2gene_reference.tsv`
- `stringtie_tx2gene_bridge.tsv`
- `tx2gene_master.tsv`
- `gene_namespace.tsv`
- `import_manifests/`

这些文件属于主流程正式输入，不应被草率手工覆盖。

## 2. 质量控制与去污染

### 原始与修剪后 QC

- `results/01.raw_qc`：原始 reads 的 FastQC 结果
- `results/02.trimmed_data`：fastp 输出、修剪后 reads、修剪后 FastQC 结果

### Decontam 主链路输出

目录：`results/03.decontam`

这一层同时包含主链路输出和旁路证据，但两者用途不同。

#### 主链路输出

- `clean/`：默认供宿主比对、定量、DEA、共识使用的 clean FASTQ
- `stats/`：样本级和项目级 decontam 决策统计
- `qc/`：clean FASTQ 的 FastQC 结果

默认语义：

1. `clean/` 只服务宿主主链路
2. 不应把微生物侧支证据重新混入主分析输入

#### 旁路证据输出

- `audit/uncertain/`：未决 reads
- `audit/removed/`：被剔除的非宿主或污染相关 reads
- `reference/`：参考污染审计与 scaffold 黑名单等证据
- `clues/tables/`：微生物线索侧支的正式汇总表
- `clues/plots/`：微生物线索侧支的正式图形输出

这些输出默认用于：

1. 审计
2. 异常样本解释
3. 微生物线索侧支分析

它们默认不作为宿主 DEA 主输入。

## 3. 宿主比对、定量与 DEA

### 比对

目录：`results/04.alignment`

- HISAT2 或 STAR 的 BAM 与索引

这是宿主主链路核心中间结果，直接依赖 `results/03.decontam/clean`。

### 定量

目录：`results/05.quantification`

- `native/`：各定量器原生输出
- `matrices/`：下游正式使用的 counts 或 abundance 矩阵
- `audit/`：轻量审计表

默认应优先把 `matrices/` 视为正式分析输入，而不是零散使用各工具原生文件。

### 单定量器差异表达

目录：`results/06.differential_expression`

每个 quantifier 子目录通常包含：

- 原始 DEA 结果表
- `dea_session.rds`
- 归一化计数矩阵

这是共识层的直接上游输入。

## 4. 共识整合与报告

### 共识输出

目录：`results/07.consensus_expression`

这是仓库的主结果层，汇总多定量器 DEA 证据后得到正式的共识差异表达结果。

### 报告输出

目录：`results/08.reports`

- MultiQC
- 其他汇总性可视化与报告材料

## 微生物线索侧支的最小输出

如果启用独立微生物线索分析，正式建议只保留以下四类输出：

1. 样本级微生物负担汇总表
2. 重点目标存在性面板
3. 微生物组成图
4. 与宿主结果并排的解释图

当前正式文件名约定为：

1. `results/03.decontam/clues/tables/sample_microbial_burden.tsv`
2. `results/03.decontam/clues/tables/priority_targets.tsv`
3. `results/03.decontam/clues/plots/microbial_composition_stacked_bar.pdf`
4. `results/03.decontam/clues/plots/host_context_overlay.pdf`

其中 `priority_targets.tsv` 的重点目标 taxid 和 `sample_microbial_burden.tsv` 的 Top taxa 保留数量由 `config.yaml` 中的 `decontam.clues.priority_targets` 与 `decontam.clues.top_taxa_n` 控制。

不要把侧支结果和宿主主结果混放在同一张正式 DEG 表里。

## 使用原则

1. 宿主主链路只吃正式主链路输出。
2. 旁路证据只做审计和侧支分析。
3. 如果某份输出未来会被长期解释或复用，就把其语义补充到本文件，而不是在 `experiments/` 中继续堆新说明。
